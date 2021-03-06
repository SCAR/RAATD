# This script identifies the start and end dates of user-defined
# breeding stages, based on fitted tracking data

# Mark Hindell & Ryan Reisinger
# Last modified 2017-09-26

#--------------------
## Clear work space 
#rm(list = ls())

#--------------------
## Set species code
sp <- "BBAL"

#--------------------
## Libraries
library(stringr)
library(geosphere)

#--------------------
## Set working directory to GitHub repo location
setwd("D:/RAATD/Data/GitHubClones/raatd_data")

#--------------------
## Source the required functions

# Source the species-specific 'config' file
source(paste0("data_split/config/", sp, "config.r"))

# Source the 'trip_split' function
source("data_split/code/trip_split.R")

# Source the 'failr' function
source("data_split/code/failr.R")

# Source the 'trip_start' function
source("data_split/code/trip_start.R")

# Source the 'stage_start' function
source("data_split/code/stage_start.R")

#--------------------
## Prerequisites

## Metadata
meta <- read.csv("metadata/SCAR_Metadata_2017_forWEBDAV.csv", stringsAsFactors = F)

## Load SSM fits
## This section copes with the old (.Rdata) and new (.RDS) SSM files

if (length(grep(".RDS", filter.rdata)) > 0) {
fits <- readRDS(filter.rdata) # SSM fits new file
} else {
  load(filter.rdata) # SSM fits
  fits <- ssm_by_id #'trip_split' uses object name 'fits'
}


## Discard any fits which failed
fits <- failr(fits)

##-----------------------------------
## Remove dates in the 3 weeks before and after Equinoxes
removeEquinox <- function(fits) {
  #Autumn equinox
  equi.A <- as.integer(format(as.POSIXct("20 March", format = "%d %B"), "%j"))
  equi.A.start <- equi.A - 21
  equi.A.end <- equi.A + 21
  #Spring equinox
  equi.S <- as.integer(format(as.POSIXct("22 September", format = "%d %B"), "%j"))
  equi.S.start <- equi.S - 21
  equi.S.end <- equi.S + 21
  #Remove equinox data
  for (i in 1:length(fits$ssm)) {
    dat <- fits$ssm[[i]]$predicted
    id <- fits$id[i]
    if (substr(id, 1, 14) == "BBAL-Kerguelen") {
    dat$jday <- as.integer(format(dat$date, "%j"))
    dat <- dat[dat$jday < equi.A.start | dat$jday > equi.S.end | (dat$jday > equi.A.end & dat$jday < equi.S.start), ]
    dat$jday <- NULL
    }
    fits$ssm[[i]]$predicted <- dat
  }
  return(fits)
}

fits <- removeEquinox(fits)
##-----------------------------------

## Get the IDs
s.list <- fits$id


## If neccessary, read in 'SOEShauloutSites.csv' for other haulout sites.
#haulouts <- read.csv("data_split/code/SOEShauloutSites.csv", stringsAsFactors = F)

#########################################################################
## Identify trips
## (function 'trip_split')

d <- trip_split(fits, s.list, meta, allow.haulout = FALSE, haulouts = NULL)
tmp <- lapply(d, function(x) x$file)
f <- do.call(rbind, tmp)

#########################################################################
## End of trip and start of next trip, assuming the individual is on shore mid-way between 
## (function 'trip_start')

t.list <- unique(f$ref)
ts <- trip_start(f, t.list)
#tmp <- lapply(ts, function(x) x$file) # unneccesary with my edit of trip_start_v2
t_start <- do.call(rbind, ts)

#########################################################################
## Start time of each stage, based on the start time of the trips 
## (function 'stage_start')

st <- stage_startR(stages, t_start, t.list, force.pb = TRUE, pb.days = 183)
tmp <- lapply(st, function(x) x$file)
sten_stage <- do.call(rbind, tmp)

#--------------------
## Save the output
saveRDS(sten_stage,
        file = paste0("./data_split/stage_dates/",
                      sp,
                      "_stages.RDS"))
