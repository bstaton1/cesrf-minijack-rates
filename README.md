# cesrf-minijack-rates
 
 This repository houses the code and data for conducting the minijack rate analyses found in the manuscript titled: _Precocious maturation of hatchery-raised spring Chinook Salmon as age-2 minijacks is not detectably affected by sire age_ by P.F. Galbreath, C.A. Stockton, C.M. Knudsen, L.R. Medeiros, I.J. Koch, B.A. Staton, W.J. Bosch, H. Nuetzel, and A.L. Pierce.

## File structure

* `00-packages.R`: loads all R packages used by the code in this repo
* `01-data-prep.R`: processes the data file `minijack-data.csv` for analysis
* `02-functions.R`: contains a handful of custom functions to facilitate parametric bootstrapping
* `03-main-analysis.R`: fits models and performs bootstrapping. Saves `.rds` files that are read in by subsequent scripts -- these output files are not tracked by this Git repo. To execute other files in this repo, you will need to run this file -- it takes several hours under the pre-populated settings.
* `04-ms-output.R`: creates files for the table and figure output found in the main text of the manuscript
* `minijack-data.csv`: contains the raw experimental data
* `power-analysis/power-simulation.R`: conducts the power analysis calculations and saves output files that are processed by code in the supplemental Rmd file
* `supplemental-material.Rmd`: source code to create the supplemental material document.

## Dependencies

All analyses were conducted in R version 4.0.2 with these package versions:

| Package    | Version | Used For                                   |
| ---------- | ------- | ------------------------------------------ |
| `glmmTMB`  | 1.0.2.1 | Primary fitting of GLMMs                   |
| `DHARMa`   | 0.3.3.0 | Residual diagnostics for GLMMs             |
| `stringr`  | 1.4.0   | Misc. string manipulations                 |
| `lme4`     | 1.1.23  | Methods for parametric bootstrap           |
| `boot`     | 1.3.25  | Methods for summarizing bootstrap          |
| `scales`   | 1.1.1   | Transparent colors on graphics             |
| `snow`     | 0.4.3   | Parallel computing                         |
| `rlecuyer` | 0.3.5   | Setting random seeds in parallel computing |
