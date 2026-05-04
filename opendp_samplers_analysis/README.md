# CKS20 sampler study

Compares three variants of the discrete Laplace sampler (original CKS20, Taylor, fixed) on sampling quality and runtime

Run things in this order:

1. `cargo run --release` runs the Rust sampler and dumps CSVs under `./csvs/` 
Default is 1 run, change source code for more. Requires OpenSSL. For windows users, might want to change the rng module and use rand instead
2. `csv_to_mat.m` (in MATLAB) converts the CSVs into `.mat` files under `./data_mats/`
3. Once the `.mats` exist, run `statistical_analysis.m` for the per-layer tests and/or `plrv_analysis.m` for the PLRV analysis
