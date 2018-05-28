# duckConfit
Crispy estimated tracks rendered from their own error-prone data. This package prepares RAATD standardised track data, fits a state-space model to the prepared data and produces quality-control plots.

This package contains the functions and workflow vignette that reproduce the RAATD data filtering process. `duckConfit` imports additional R packages: `dplyr`, `tibble`, `readr`, `lubridate`, `ggplot2`, `gridExtra`, `diveMove`, `TMB`. All of these are available on CRAN. `duckConfit` also requires the package `ssmTMB`, which is available from [GitHub](https://github.com/ianjonsen/ssmTMB) and can be installed in R via: `devtools::install_github("ianjonsen/ssmTMB")`

