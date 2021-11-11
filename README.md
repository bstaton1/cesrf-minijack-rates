# cesrf-minijack-rates

This repository houses the code and data for conducting the minijack rate analyses found in the manuscript titled: _Precocious maturation of hatchery-raised spring Chinook Salmon as age-2 minijacks is not detectably affected by sire age_ by P. F. Galbreath, H. M. Nuetzel, B. A. Staton, C. A. Stockton, C. M. Knudsen, L. R. Medeiros, I. J. Koch, W. J. Bosch, and A. L. Pierce published in _Transactions of the American Fisheries Society_.

[![ArticleDOI](https://img.shields.io/badge/Article%20DOI-10.1002%2Ftafs.10343-blue)](https://doi.org/10.1002/tafs.10343)

[![GitHub Repo Archive DOI](https://img.shields.io/badge/GitHub%20Repo%20Archive%20DOI-10.5281%2Fzenodo.4730682-blue)](https://doi.org/10.5281/zenodo.4730682)

## File structure

* `00-packages.R`: loads all R packages used by the code in this repo
* `01-data-prep.R`: processes the data file `minijack-data.csv` for analysis
* `02-functions.R`: contains a handful of custom functions to facilitate parametric bootstrapping
* `03-main-analysis.R`: fits models and performs bootstrapping. Saves `.rds` files that are read in by subsequent scripts -- these output files are not tracked by this Git repo. To execute other files in this repo, you will need to run this file -- it takes several hours under the pre-populated settings.
* `04-ms-output.R`: creates files for the table and figure output found in the main text of the manuscript
* `minijack-data.csv`: contains the raw experimental data (variables described in table below)
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

## Variables in Data

The file `minijack-data.csv` contains all of the raw data for the analysis. Each row is an individual progeny with variables shown in the table below:

| Variable             | Description                                                  |
| -------------------- | ------------------------------------------------------------ |
| `year`               | The year the cross was made                                  |
| `dam_id`             | A unique identifier for the dam (female parent)              |
| `dam_age`            | The total age of the dam (winters after the dam was spawned) |
| `dam_FL`             | The fork length of the dam (cm)                              |
| `dam_POH`            | The dam post-orbital hypural length (cm)                     |
| `dam_wt`             | The weight of the dam (kg)                                   |
| `egg_wt`             | The average weight of each of the dam's eggs (g)             |
| `sire_id`            | A unique identifier for the sire (male parent)               |
| `sire_age`           | The total age of the sire (winters after the sire was spawned) |
| `sire_FL`            | The sire fork length (cm)                                    |
| `sire_POH`           | The sire post-orbital hypural length (cm)                    |
| `sire_wt`            | The weight of the sire (kg)                                  |
| `spawn_date`         | The date of spawning (M/D/YYYY)                              |
| `spawn_doy`          | The day of the year of spawning                              |
| `cross_type`         | `dam_age` x `sire_age`                                       |
| `progeny_mat_status` | The maturation status of the progeny; "Minijack" or "Non-maturing" |
| `progeny_length`     | The total length of the progeny (mm)                         |
| `progeny_wt`         | The weight of the progeny (g)                                |