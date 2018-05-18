## Trim individual tracks
## Ryan Reisinger
## Last modified 2018-05-18

## ############################################

setwd("D:/RAATD/Data/GitHubClones/raatd_data")
source("D:/RAATD/Data paper/trimR.R")

## Make sure that the development version of ggmap is installed,
## the CRAN version will fail
# library(devtools)
# devtools::install_github("dkahle/ggmap")

## Run trimR function

sp = "SOES" # RAATD abbreviation for the species of interest
id = "ct16_241_06" # Individual of interest

rm(dat.new)

dat.new <- trimR(sp = sp, id = id)

## Write output back to repo

write.csv(dat.new,
          paste0("D:/RAATD/Data/GitHubClones/raatd_data/data_raw_trimmed/RAATD2017_",
                 sp, ".csv"),
          row.names = F, na = "")