## trimR
## Function to trim a single individual's track in the RAATD data
## Ryan Reisinger
## Last modified 2018-05-18

trimR <- function(sp, id){
  
  ##'sp' = four letter species abbreviation (character)
  ##'id' = individual_id (character)

  ## Returns the dataframe
  
#-----------------------------------

library(geosphere)
library(lubridate)
library(ggmap)
library(ggplot2)

## Read in the metadata
met <- read.csv("D:/RAATD/Data/GitHubClones/raatd_data/metadata/SCAR_Metadata_2017_forWEBDAV.csv", stringsAsFactors = F)

## Define species of interest, using four letter RAATD abbreviations
## for abbreviaitons, see
## unique(met$abbreviated_name)

## Read in data for the species of interest
dat <- read.csv(paste0("D:/RAATD/Data/GitHubClones/raatd_data/data_raw_trimmed/RAATD2017_", sp, ".csv"), stringsAsFactors = F)

#-----------------------------------
## Trim start and end portions of tracks using locator
#-----------------------------------

## This section of code loops through the data by individual.
## It plots displacement v. time,
## and then allows the user to click the locator on the plot twice,
## defining the start and end of the track section TO KEEP.
## This will be reflected in the column 'location_to_keep'.


## Be sure to have your plot window large enough before starting

  hold <- dat[dat$individual_id == id, ]

if (nrow(hold) == 0) {
  stop("No data selected")
}
  
  #Add Distance
  hold$Distance <- distGeo(p1 = met[met$individual_id == id, c("deployment_decimal_longitude", "deployment_decimal_latitude")], p2 = hold[ , c("decimal_longitude", "decimal_latitude")])
  hold$Distance <- hold$Distance/1000
  
  #Add date.time
  hold$date.time <- ymd_hms(paste0(hold$year,"-",hold$month,"-",hold$day," ",hold$time),quiet = TRUE,truncated = 1)
  m <- met[met$individual_id == id, ]
  
  if (m$deployment_time == "") {
    m$deployment_time <- "00:00:00"
  }
  
  date.deploy <- ymd_hms(paste0(m$year,"-",m$month,"-",m$day," ",m$deployment_time),quiet = TRUE,truncated = 1)
  hold$dif <- difftime(hold$date.time, date.deploy, units = "days")
  hold$dif <- as.numeric(hold$dif)
  
  #Now trim
  hold$d <- as.numeric(hold$date.time)
  
  #Plot
  p1 <- ggplot(data = hold, aes(x = hold$dif, y = hold$Distance, color = as.factor(location_to_keep))) +
    geom_point(show.legend = T) +
    theme_bw() +
    geom_vline(xintercept = 0) +
    guides(color = guide_legend(title = "location_to_keep")) +
    geom_line(color = "gray40") + 
    labs(x = "Days since deployment", y = "Distance from deployment")
  print(p1)
  
  #Pick start and end using locator
  t <- gglocator(2, mercator = F)
  t <- t[ , 1]
  
  #Assign 'location_to_keep' based on locator
  for (j in 1:nrow(hold)) {
    if (hold$dif[j] > t[1] && hold$dif[j] < t[2]) {
      hold$location_to_keep[j] <- 1
    } else {
      hold$location_to_keep[j] <- 0
    }
  }
  
  #Remove extraneous columns
  hold$d <- NULL
  hold$dif <- NULL
  hold$date.time <- NULL
  hold$Distance <- NULL

if (nrow(dat[dat$individual_id == id, ]) == nrow(hold)) {
  dat[dat$individual_id == id, "location_to_keep"] <- hold$location_to_keep
} else {
  stop("Length of 'location_to_keep' in new data did not match length of original data")
}

return(dat)
  
}
