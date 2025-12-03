# Boltz2_IpSAE

This repository contains scripts for calculating and visualizing protein-protein interaction scores using `ipsae.py`. The `ipsae.py` script calculates several scores, including ipSAE, pDockQ, and LIS, which are used to evaluate the quality of predicted protein structures from tools like AlphaFold2, AlphaFold3, and Boltz1.

## Example Output

![Example Plot](placeholder.png)

*(This is a placeholder image. To generate the actual plot, please follow the instructions below.)*

## How to Run

1.  **Clone the repository and set up the `boltz2` environment.**
2.  **Navigate to the example directory:**

    ```bash
    cd ./Boltz2_IpSAE/example_run
    ```

3.  **Run the main script to generate the validation scripts:**

    ```bash
    ./run.sh
    ```

4.  **Run the co-folding validation scripts:**

    ```bash
    ./boltz_validation/run_all_cofolding.sh
    ```

5.  **Visualize the results:**

    ```bash
    ./boltz_validation/visualise_cofolding_results.sh
    ```

This will generate a plot of the results, which will be saved in the `boltz_validation` directory.
