# Code and data for R-wiring

This repository contains the data and analysis code supporting the manuscript entitled “R-wiring: An R package for quantifying interaction rewiring and species turnover in ecological networks”.

The datasets used in the three case studies are located in the `DATA` folder and were obtained from the source studies cited in the manuscript.

The analyses were implemented in R and Python. Instructions for reproducing the figures and tables are provided below.

## R analyses

Case study 1 and the combined Figure 4 were implemented in R. All R scripts are located in the `R` folder and should be run with the repository root as the working directory.

### Execution order

To generate Figure 2 and Figures S1–S4, run the following scripts in order:

1. `01_prepare_pollobs_data.R`
2. `02_analyse_pollobs_rwiring.R`
3. `03_plot_pollobs_results.R`

The first script also calculates the network summaries and sample coverage statistics used in Tables S1 and S2.

Additional summaries for case study 1 can be obtained by running `04_analyse_pollobs_summary_results.R` after `02_analyse_pollobs_rwiring.R`.

To generate Figure S5, run `05a_analyse_temporal_lags_pollobs.R` after completing the data preparation and analysis steps above.

To generate Figure 4, run the following scripts in order:

1. `05a_analyse_temporal_lags_pollobs.R`
2. `05b_analyse_temporal_lags_beefun.R`
3. `05c_analyse_spatial_lags_clough.R`
4. `06_plot_combined_reorganization_cases.R`

If `05a_analyse_temporal_lags_pollobs.R` has already been run to generate Figure S5, proceed with the remaining scripts.

### Input and output files

The scripts read their input data from `DATA` and save processed data and analysis results within `DATA`.

Figures are saved in the `PollObs`, `BeeFun`, and `Clough` subfolders of `FIGURES`. Figure 4 is saved in `FIGURES/Combined`, and supplementary figures for case study 1 are saved in `FIGURES/PollObs/Supplementary`.

### R dependencies and installation

The analyses require the following R packages:

- Rwiring
- dplyr
- tidyr
- purrr
- readr
- lubridate
- iNEXT
- ggplot2
- patchwork
- geosphere

Install the dependencies using:

```r
install.packages(c(
  "remotes", "dplyr", "tidyr", "purrr", "readr", "lubridate",
  "iNEXT", "ggplot2", "patchwork", "geosphere"
))
remotes::install_github("AlfonsoAllen/Rwiring")
```

Further documentation and examples for Rwiring are available in the [GitHub repository](https://github.com/AlfonsoAllen/Rwiring).

### R environment

The R version, system information, and versions of the packages loaded for the analyses are provided in `sessionInfo.txt`.

This file documents the environment used in the study. The installation commands above do not automatically reproduce the recorded package versions.

## Python analyses

Case studies 2 and 3 were implemented in Python.

### Notebooks

To generate Figure 3, Figures S6–S9, and Tables S10 and S11, run `Rewiring_temporal.ipynb`.

To generate Figure S10 and Tables 1 and S12, run `Rewiring_spatial.ipynb`.

### Input and output files

Both notebooks read the input datasets from `DATA` and save the generated figures and tables in `FIGURES`.

### Python environment

The specifications of the Conda environment used for the Python analyses are provided in `spec-file.txt`.
