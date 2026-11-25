# FDBAH

R code for *Initial Subset-Size Selection for Minimum Covariance Determinant Estimation*, by Chikun Li and Baisuo Jin. Figure and table numbers follow the English manuscript of 25 November 2026.

## Install

From the repository root (the directory containing `DESCRIPTION`), run:

```sh
Rscript inst/reproduce/setup.R
R CMD INSTALL .
```

Alternatively, install the source archive with `R CMD INSTALL FDBAH_0.1.0.tar.gz` after installing the dependencies above.

For Tables 4 and S6, also run `Rscript inst/reproduce/setup.R --comparators`. This installs the included authors' StableMCD 1.0 source and builds fsdac 0.6-1 from [FSDA commit 85ede5f](https://github.com/UniprJRC/FSDA-MATLAB_Coder/tree/85ede5feeb450ac0329ae0c1e420d9430aee8980). A C/C++ toolchain for R is required (Rtools on Windows). The build omits two unused inverse-envelope interfaces whose wrapper headers are absent from that commit. To supply the pinned upstream ZIP locally, set `FDBAH_FSDA_ZIP`.

The three real-data scripts use Python 3 with NumPy and pandas for reading the original data files. Install those dependencies with `python -m pip install numpy pandas`. Set `FDBAH_PYTHON` to the interpreter path if needed. The original computations used R 4.5.2, ddalpha 1.3.16, robustbase 0.99-6 and MASS 7.3-65.

## Main function

```r
fit <- FDBAH::fdb_fit(x, seed = 1L)
fit$h         # initial subset size
fit$initial   # row indices before reweighting
fit$retained  # row indices after reweighting
```

`x` is a numeric observation-by-variable matrix. `fit$center` and `fit$cov` are the final estimates. `FDBAH::fdb_four(x, seed = 1L)` fits FDB-AH, FDB-FS, FDB-HF and FDB-75 using a common projection-depth ordering.

## Reproduce one item

For example, run `Rscript inst/reproduce/Table1.R`. Each script can run independently. To find the scripts after installation, use `system.file("reproduce", package = "FDBAH")` in R. Output goes to `output/` in the working directory. Set `FDBAH_OUTPUT` to choose a different directory and `FDBAH_WORKERS` to set the worker count (default 1).

| Script | Design | Output |
|---|---|---|
| Figure1.R | SKAB normal and rotor records, three FDB methods, 100 projection runs. The upper panels use run 1. | Figure1.pdf |
| TableS6.R | Same two SKAB records, six methods, 100 runs each. | TableS6.csv |
| Table4.R | 2,500 Skin pixels, six methods, 100 runs each. | Table4.csv |
| Table1.R | Uncontaminated Gaussian samples, settings A, B and C, 1,000 replicates each. | Table1.csv |
| Table2.R | Four contamination designs, settings A, B and C, contamination 0.10, signal 15, 1,000 replicates per configuration. | Table2.csv |
| Figure3.R | Mean initial subset size in setting B, four designs, two contamination levels and nine signal strengths, 1,000 replicates per configuration. | Figure3.pdf and Figure3.csv |

Simulation scripts generate their own observations. Real-data scripts obtain the original SKAB and ADBench files and verify SHA-256 before selecting rows. For offline use, set `FDBAH_DATA_DIR` to a directory containing the required originals: `anomaly-free.csv` and `7.csv` for SKAB, or `33_skin.npz` for Skin. The two figures use a Cairo PDF device and Times New Roman when available.

Running times in Tables 4 and S6 are measured on the current machine. The StableMCD projection code uses C++ `std::random_device` in addition to the R bootstrap seed, so its results may vary across fresh runs. Each script starts directly from the specified input, without archived fits or timing records.
