import tempfile
import unittest
from pathlib import Path

from hot_boltz_worker import HotBoltzWorker, WorkerConfig, discover_jobs, prediction_exists, run_batch


class FakeWorker(HotBoltzWorker):
    instances = []

    def __init__(self, config):
        super().__init__(config)
        self.runs = []
        self.instances.append(self)

    def run(self, input_yaml, output_dir):
        self.runs.append((input_yaml, output_dir))
        prediction_dir = output_dir / f"boltz_results_{input_yaml.stem}" / "predictions" / input_yaml.stem
        prediction_dir.mkdir(parents=True, exist_ok=True)
        (prediction_dir / f"{input_yaml.stem}_model_0.cif").write_text("data_test\n")
        return {"status": "complete", "return_code": 0}


class HotBoltzWorkerTests(unittest.TestCase):
    def test_batch_runs_pending_jobs_once_and_skips_completed_jobs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            binder_a = root / "binder_a"
            binder_b = root / "binder_b"
            binder_a.mkdir()
            binder_b.mkdir()
            yaml_a = binder_a / "binder_a_vs_target.yaml"
            yaml_b = binder_b / "binder_b_vs_target.yaml"
            yaml_a.write_text("version: 1\n")
            yaml_b.write_text("version: 1\n")
            output_b = binder_b / "outputs"
            prediction_dir = output_b / "boltz_results_binder_b_vs_target" / "predictions" / "binder_b_vs_target"
            prediction_dir.mkdir(parents=True)
            (prediction_dir / "binder_b_vs_target_model_0.cif").write_text("data_test\n")

            FakeWorker.instances.clear()
            summary = run_batch(root, WorkerConfig(), worker_factory=FakeWorker)

            self.assertEqual(summary, {"total": 2, "completed": 1, "skipped": 1})
            self.assertEqual(len(FakeWorker.instances), 1)
            self.assertEqual(FakeWorker.instances[0].runs, [(yaml_a, binder_a / "outputs")])

    def test_all_completed_jobs_do_not_start_worker(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            binder = root / "binder_a"
            binder.mkdir()
            input_yaml = binder / "binder_a_vs_target.yaml"
            input_yaml.write_text("version: 1\n")
            output_dir = binder / "outputs"
            prediction_dir = output_dir / "boltz_results_binder_a_vs_target" / "predictions" / "binder_a_vs_target"
            prediction_dir.mkdir(parents=True)
            (prediction_dir / "model_0.cif").write_text("data_test\n")

            FakeWorker.instances.clear()
            summary = run_batch(root, WorkerConfig(), worker_factory=FakeWorker)

            self.assertEqual(summary, {"total": 1, "completed": 0, "skipped": 1})
            self.assertEqual(FakeWorker.instances, [])

    def test_job_discovery_and_prediction_detection(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            binder = root / "binder_a"
            binder.mkdir()
            input_yaml = binder / "binder_a_vs_target.yaml"
            input_yaml.write_text("version: 1\n")
            output_dir = binder / "outputs"

            self.assertEqual(discover_jobs(root), [(input_yaml, output_dir)])
            self.assertFalse(prediction_exists(output_dir, input_yaml))


if __name__ == "__main__":
    unittest.main()
