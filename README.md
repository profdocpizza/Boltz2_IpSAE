# Boltz2_IpSAE

This repository contains scripts for co-folding using Boltz2, calculating and visualising [ipSAE](https://github.com/DunbrackLab/IPSAE) metric. With a single .yaml file you can cofold your binder to multiple targets and itself to monitor specificity of binding.

## How to Run

1.  **Clone the repository and set up the [`Boltz2`](https://github.com/jwohlwend/boltz) environment.**
2.  **Navigate to the example directory:**

    ```bash
    cd ./Boltz2_IpSAE/example_run
    ```

3.  **Run the main script to generate the validation scripts:**

    ```bash
    python ../make_binder_validation_scripts.py --config config.yaml
    ```

4.  **Run all co-folding validation scripts with a single script**

    ```bash
    ./boltz_validation/run_all_cofolding.sh
    ```

5.  **When co-folding is done, visualize the results:**

    ```bash
    ./boltz_validation/visualise_cofolding_results.sh
    ```

This will generate a plot of the results, which will be saved in the `boltz_validation/summary` directory.
