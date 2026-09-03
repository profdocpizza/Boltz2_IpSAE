#!/usr/bin/env python3
"""Run a directory of Boltz2 YAML jobs with one persistent model instance.

The regular ``boltz predict`` CLI is a process-level entry point and loads the
checkpoint once per invocation.  This module keeps the model in memory and
creates only the per-job preprocessing, datamodule, writer, and Trainer
objects.  Jobs are deliberately executed serially because they target one
GPU.
"""

from __future__ import annotations

import argparse
import os
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Callable, Iterable


@dataclass(frozen=True)
class WorkerConfig:
    """Boltz settings required to reproduce the generated CLI commands."""

    cache: Path | None = None
    recycling_steps: int = 3
    sampling_steps: int = 200
    diffusion_samples: int = 1
    max_parallel_samples: int | None = 5
    step_scale: float | None = None
    write_full_pae: bool = False
    write_full_pde: bool = False
    output_format: str = "mmcif"
    num_workers: int = 2
    preprocessing_threads: int = 1
    override: bool = False
    use_msa_server: bool = False
    msa_server_url: str = "https://api.colabfold.com"
    msa_pairing_strategy: str = "greedy"
    use_potentials: bool = False
    no_kernels: bool = False
    max_msa_seqs: int = 8192
    subsample_msa: bool = False
    num_subsampled_msa: int = 1024
    write_embeddings: bool = False
    checkpoint: Path | None = None


def _prediction_dir(output_dir: Path, input_yaml: Path) -> Path:
    job_root = output_dir / f"boltz_results_{input_yaml.stem}"
    return job_root / "predictions" / input_yaml.stem


def prediction_exists(output_dir: Path, input_yaml: Path) -> bool:
    """Return whether a job has at least one structure prediction."""
    prediction_root = _prediction_dir(output_dir, input_yaml)
    if not prediction_root.is_dir():
        return False
    return any(
        path.is_file() and path.suffix.lower() in {".cif", ".mmcif", ".pdb"}
        for path in prediction_root.iterdir()
    )


def discover_jobs(output_root: Path) -> list[tuple[Path, Path]]:
    """Discover generated YAML jobs and their existing per-binder output dirs."""
    jobs: list[tuple[Path, Path]] = []
    for input_yaml in sorted(output_root.glob("binder_*/*.yaml")):
        if input_yaml.is_file():
            jobs.append((input_yaml, input_yaml.parent / "outputs"))
    return jobs


class HotBoltzWorker:
    """Persistent Boltz2 inference worker for one-GPU serial execution."""

    def __init__(self, config: WorkerConfig) -> None:
        self.config = config
        self.cache: Path | None = None
        self.model: Any = None
        self._boltz: dict[str, Any] = {}

    def start(self) -> None:
        """Load Boltz data and model weights exactly once."""
        if self.model is not None:
            return

        import torch
        from boltz.main import (
            Boltz2,
            Boltz2DiffusionParams,
            BoltzSteeringParams,
            MSAModuleArgs,
            PairformerArgsV2,
            download_boltz2,
            get_cache_path,
        )
        torch.set_grad_enabled(False)
        torch.set_float32_matmul_precision("highest")

        cache = (self.config.cache or Path(get_cache_path())).expanduser().resolve()
        cache.mkdir(parents=True, exist_ok=True)
        download_boltz2(cache)

        diffusion_params = Boltz2DiffusionParams()
        if self.config.step_scale is not None:
            diffusion_params.step_scale = self.config.step_scale

        msa_args = MSAModuleArgs(
            subsample_msa=self.config.subsample_msa,
            num_subsampled_msa=self.config.num_subsampled_msa,
        )
        steering_args = BoltzSteeringParams()
        steering_args.fk_steering = self.config.use_potentials
        steering_args.physical_guidance_update = self.config.use_potentials

        checkpoint = self.config.checkpoint or cache / "boltz2_conf.ckpt"
        predict_args = {
            "recycling_steps": self.config.recycling_steps,
            "sampling_steps": self.config.sampling_steps,
            "diffusion_samples": self.config.diffusion_samples,
            "max_parallel_samples": self.config.max_parallel_samples,
            "write_confidence_summary": True,
            "write_full_pae": self.config.write_full_pae,
            "write_full_pde": self.config.write_full_pde,
        }
        self.model = Boltz2.load_from_checkpoint(
            checkpoint,
            strict=True,
            predict_args=predict_args,
            map_location="cpu",
            diffusion_process_args=asdict(diffusion_params),
            ema=False,
            use_kernels=not self.config.no_kernels,
            pairformer_args=asdict(PairformerArgsV2()),
            msa_args=asdict(msa_args),
            steering_args=asdict(steering_args),
        )
        self.model.eval()
        self.cache = cache
        self._boltz = {
            "Boltz2InferenceDataModule": self._import("Boltz2InferenceDataModule"),
            "BoltzProcessedInput": self._import("BoltzProcessedInput"),
            "BoltzWriter": self._import("BoltzWriter"),
            "Manifest": self._import("Manifest"),
            "Trainer": self._import("Trainer", "pytorch_lightning"),
            "check_inputs": self._import("check_inputs"),
            "filter_inputs_structure": self._import("filter_inputs_structure"),
            "process_inputs": self._import("process_inputs"),
        }
        print(f"[hot-boltz] ready; checkpoint={checkpoint}", flush=True)

    @staticmethod
    def _import(name: str, module_name: str = "boltz.main") -> Any:
        module = __import__(module_name, fromlist=[name])
        return getattr(module, name)

    def run(self, input_yaml: Path, output_dir: Path) -> dict[str, Any]:
        """Run one YAML job using the already-loaded model."""
        self.start()
        assert self.model is not None
        assert self.cache is not None

        input_yaml = input_yaml.resolve()
        output_dir = output_dir.resolve()
        job_root = output_dir / f"boltz_results_{input_yaml.stem}"
        job_root.mkdir(parents=True, exist_ok=True)

        check_inputs: Callable[[Path], list[Path]] = self._boltz["check_inputs"]
        process_inputs = self._boltz["process_inputs"]
        manifest_cls = self._boltz["Manifest"]
        filter_inputs_structure = self._boltz["filter_inputs_structure"]

        data = check_inputs(input_yaml)
        process_inputs(
            data=data,
            out_dir=job_root,
            ccd_path=self.cache / "ccd.pkl",
            mol_dir=self.cache / "mols",
            use_msa_server=self.config.use_msa_server,
            msa_server_url=self.config.msa_server_url,
            msa_pairing_strategy=self.config.msa_pairing_strategy,
            msa_server_username=os.environ.get("BOLTZ_MSA_USERNAME"),
            msa_server_password=os.environ.get("BOLTZ_MSA_PASSWORD"),
            api_key_header=None,
            api_key_value=os.environ.get("MSA_API_KEY_VALUE"),
            boltz2=True,
            preprocessing_threads=self.config.preprocessing_threads,
            max_msa_seqs=self.config.max_msa_seqs,
        )

        manifest = manifest_cls.load(job_root / "processed" / "manifest.json")
        filtered_manifest = filter_inputs_structure(
            manifest=manifest,
            outdir=job_root,
            override=self.config.override,
        )
        if not filtered_manifest.records:
            return {
                "status": "skipped_existing",
                "return_code": 0,
                "input_yaml": str(input_yaml),
                "output_dir": str(output_dir),
            }

        processed_dir = job_root / "processed"
        processed = self._boltz["BoltzProcessedInput"](
            manifest=filtered_manifest,
            targets_dir=processed_dir / "structures",
            msa_dir=processed_dir / "msa",
            constraints_dir=(processed_dir / "constraints") if (processed_dir / "constraints").exists() else None,
            template_dir=(processed_dir / "templates") if (processed_dir / "templates").exists() else None,
            extra_mols_dir=(processed_dir / "mols") if (processed_dir / "mols").exists() else None,
        )
        writer = self._boltz["BoltzWriter"](
            data_dir=processed.targets_dir,
            output_dir=job_root / "predictions",
            output_format=self.config.output_format,
            boltz2=True,
            write_embeddings=self.config.write_embeddings,
        )
        trainer = self._boltz["Trainer"](
            default_root_dir=job_root,
            strategy="auto",
            callbacks=[writer],
            accelerator="gpu",
            devices=1,
            precision="bf16-mixed",
        )
        data_module = self._boltz["Boltz2InferenceDataModule"](
            manifest=processed.manifest,
            target_dir=processed.targets_dir,
            msa_dir=processed.msa_dir,
            mol_dir=self.cache / "mols",
            num_workers=self.config.num_workers,
            constraints_dir=processed.constraints_dir,
            template_dir=processed.template_dir,
            extra_mols_dir=processed.extra_mols_dir,
            override_method=None,
        )
        trainer.predict(self.model, datamodule=data_module, return_predictions=False)
        if not prediction_exists(output_dir, input_yaml):
            raise RuntimeError(f"Boltz2 produced no structure predictions for {input_yaml.name}")
        return {
            "status": "complete",
            "return_code": 0,
            "input_yaml": str(input_yaml),
            "output_dir": str(output_dir),
        }


def run_batch(
    output_root: Path,
    config: WorkerConfig,
    worker_factory: type[HotBoltzWorker] = HotBoltzWorker,
) -> dict[str, int]:
    """Execute all generated jobs serially, skipping completed jobs."""
    jobs = discover_jobs(output_root.resolve())
    pending = [
        (input_yaml, output_dir)
        for input_yaml, output_dir in jobs
        if config.override or not prediction_exists(output_dir, input_yaml)
    ]
    print(f"[hot-boltz] jobs={len(jobs)} pending={len(pending)} workers=1", flush=True)
    if not pending:
        return {"total": len(jobs), "completed": 0, "skipped": len(jobs)}

    worker = worker_factory(config)
    completed = 0
    skipped = len(jobs) - len(pending)
    for index, (input_yaml, output_dir) in enumerate(jobs, start=1):
        if not config.override and prediction_exists(output_dir, input_yaml):
            print(f"[hot-boltz] [{index}/{len(jobs)}] skip {input_yaml}", flush=True)
            continue
        print(f"[hot-boltz] [{index}/{len(jobs)}] run {input_yaml}", flush=True)
        result = worker.run(input_yaml, output_dir)
        if result.get("status") == "skipped_existing":
            skipped += 1
        else:
            completed += 1
    return {"total": len(jobs), "completed": completed, "skipped": skipped}


def _parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--cache", type=Path, default=None)
    parser.add_argument("--checkpoint", type=Path, default=None)
    parser.add_argument("--recycling-steps", type=int, default=3)
    parser.add_argument("--sampling-steps", type=int, default=200)
    parser.add_argument("--diffusion-samples", type=int, default=1)
    parser.add_argument("--max-parallel-samples", type=int, default=5)
    parser.add_argument("--step-scale", type=float, default=None)
    parser.add_argument("--write-full-pae", action="store_true")
    parser.add_argument("--write-full-pde", action="store_true")
    parser.add_argument("--num-workers", type=int, default=2)
    parser.add_argument("--preprocessing-threads", type=int, default=max(os.cpu_count() or 1, 1))
    parser.add_argument("--override", action="store_true")
    parser.add_argument("--use-msa-server", action="store_true")
    parser.add_argument("--msa-server-url", default="https://api.colabfold.com")
    parser.add_argument("--msa-pairing-strategy", default="greedy")
    parser.add_argument("--use-potentials", action="store_true")
    parser.add_argument("--no-kernels", action="store_true")
    parser.add_argument("--max-msa-seqs", type=int, default=8192)
    parser.add_argument("--subsample-msa", action="store_true")
    parser.add_argument("--num-subsampled-msa", type=int, default=1024)
    parser.add_argument("--write-embeddings", action="store_true")
    return parser.parse_args(list(argv) if argv is not None else None)


def main(argv: Iterable[str] | None = None) -> int:
    args = _parse_args(argv)
    config = WorkerConfig(
        cache=args.cache,
        checkpoint=args.checkpoint,
        recycling_steps=args.recycling_steps,
        sampling_steps=args.sampling_steps,
        diffusion_samples=args.diffusion_samples,
        max_parallel_samples=args.max_parallel_samples,
        step_scale=args.step_scale,
        write_full_pae=args.write_full_pae,
        write_full_pde=args.write_full_pde,
        num_workers=args.num_workers,
        preprocessing_threads=args.preprocessing_threads,
        override=args.override,
        use_msa_server=args.use_msa_server,
        msa_server_url=args.msa_server_url,
        msa_pairing_strategy=args.msa_pairing_strategy,
        use_potentials=args.use_potentials,
        no_kernels=args.no_kernels,
        max_msa_seqs=args.max_msa_seqs,
        subsample_msa=args.subsample_msa,
        num_subsampled_msa=args.num_subsampled_msa,
        write_embeddings=args.write_embeddings,
    )
    try:
        summary = run_batch(args.output_root, config)
    except Exception as exc:  # noqa: BLE001
        print(f"[hot-boltz] failed: {exc}", file=sys.stderr)
        return 1
    print(
        f"[hot-boltz] complete={summary['completed']} skipped={summary['skipped']} total={summary['total']}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
