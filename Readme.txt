The code and data provided here are a suppot for 
R-wiring: An R package for quantifying interaction rewiring and species turnover in ecological networks

The databases for the three case studies are located in the DATA folder. 
They have been downloaded from the cited studies.

The code is written in R and Python. 

### R ###########################################

Case study 1 and the combined Figure 4 were implemented in R. All scripts are located in the R folder and should be run with the repository root as the working directory.

To generate Figure 2 and Figures S1–S4, run the following scripts in order:
- `01_prepare_pollobs_data.R`
- `02_analyse_pollobs_rwiring.R`
- `03_plot_pollobs_results.R`

The first script also calculates the network summaries and sample coverage statistics used in Tables S1 and S2. Additional summaries for case study 1 can be obtained by running `04_analyse_pollobs_summary_results.R` after `02_analyse_pollobs_rwiring.R`.

To generate Figure S5, run `05a_analyse_temporal_lags_pollobs.R` after completing the data preparation and analysis steps above.

To generate Figure 4, run `05a_analyse_temporal_lags_pollobs.R`, `05b_analyse_temporal_lags_beefun.R`, and `05c_analyse_spatial_lags_clough.R`, followed by `06_plot_combined_reorganization_cases.R`.

The scripts read their input data from `DATA` and save processed data and analysis results within `DATA`. Figures are saved in the `PollObs`, `BeeFun`, and `Clough` subfolders of `FIGURES`, while Figure 4 is saved in `FIGURES/Combined`. Supplementary figures for case study 1 are saved in `FIGURES/PollObs/Supplementary`.

R configuration: The required packages are Rwiring, dplyr, tidyr, purrr, readr, lubridate, iNEXT, ggplot2, patchwork, and geosphere. Install the dependencies and Rwiring using:

```r
install.packages(c(
  "remotes", "dplyr", "tidyr", "purrr", "readr", "lubridate",
  "iNEXT", "ggplot2", "patchwork", "geosphere"
))

remotes::install_github("AlfonsoAllen/Rwiring")
```

Further documentation and examples for Rwiring are available in the [GitHub repository](https://github.com/AlfonsoAllen/Rwiring).

The R version, system information, and versions of the packages loaded for the analyses are provided in `sessionInfo.txt`.

###################################################

### PYTHON ######################################3

Case studies 2 and 3 have been coded in python. 

To generate Figure 3, Figures S6-S9, and tables S10 and S11 run the notebook Rewiring_temporal.ipynb 

To generate Figure S10 and Tables 1 and S12, run the notebook Rewiring_spatial.ipynb. 

Both notebooks take the databases from folder DATA, and creates the figures and tables in the folder FIGURES. 

Python configuration: A file with the specifications of the anaconda environment used to generate the figures is also provided, named spec-file.txt. If needed it can be used to re-create the computation environment.

#####################################################3
