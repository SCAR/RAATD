# Function to find the start and end time of each breeding stage

# Mark Hindell & Ryan Reisinger
# Last modified 2017-08-22

stage_startR <- function(stages, t_start, t.list, force.pb = FALSE, pb.days){
  stage_times <- lapply(1:length(t.list), function(i){
    cat(sprintf('\nprocessing %d of %d\n',i,length(t.list))); flush.console()
    
    stg <- t_start[t_start$ref == t.list[i],]
    #f2 <- f[f$ref==t.list[i],]

##set up julian days for the stages that span Jan 1
stages$jday <- stages$start_day
#stages$jday[stages$start_day<=180] <- stages$start_day[stages$start_day<=180]+365
#stages <- stages[order(stages$jday),]                     

##set up a new data frame for the start and end times
sten_stage <- data.frame(id = NULL, stage = NULL, start = NULL, end = NULL)

##get maximum distance and start time for each trip
#t <- aggregate(dnh~ref+trip, data=f, max, na.rm=T) #maximum dfh for each trip
#names(t) <- c("ref", "trip", "maxd")
#s <- aggregate(gmt~ref+trip, data=f, min) #start date of each trip
names(stg) <- c("ref", "trip", "start", "end")

##calculate duration of each trip
stg$dur <- as.numeric(difftime(stg$end, stg$start, units = "hours"))

##Calculate start day
stg  <- stg[order(stg$ref, stg$start),]
stg$doy <- as.POSIXlt(stg$start)$yday

#--------------------------------
## Here, pull in original data and calculate number of
## locations falling into each stage

#calculate day of year
dat <- f[f$ref == t.list[i], ]
dat$doy <- as.POSIXlt(dat$gmt)$yday

#assign a stage to each location
dat$stage <- NA

for (j in 1:dim(stages)[1]) {
  s <- stages[j,]
  if (stages$start_day[j] > stages$end_day[j]) {
    dat$stage[which(dat$doy >= s$start_day | dat$doy <= s$end_day)] <- as.character(s$stage)
  }
  if (stages$start_day[j] < stages$end_day[j]) {
    dat$stage[which(dat$doy >= s$start_day & dat$doy <= s$end_day)] <- as.character(s$stage)
  }
}

#aggregate
dat <- aggregate(ref ~ stage + trip, data = dat, FUN = length)

#find stage with the greatest number of locations for each trip
hld <- list()
for (i in dat$trip) {
  h <- dat[dat$trip == i, ]
  h <- h[h$ref == max(h$ref), c("trip", "stage")]
hld[i] <- list(h)
}

hld <- do.call(rbind, hld)
names(hld) <- c("trip", "tn")

f3 <- merge(stg, hld, by = "trip")

##for tracks that last for than 1 year need to keep track of stages in each year
##do this by adding a unique number to each stage
f3$tn1 <- as.numeric(as.factor(f3$tn))
f3$tn2 <- ceiling(cumsum(abs(c(0, diff(f3$tn1))))) + 1

##for wanderers only
##if trip last for more than year assign status to post-breeding
if (force.pb) {
f3$tn[f3$dur > (pb.days*24)] <- "post-breeding"
}


#sten_trip <- aggregate(start~trip, data=f3, range)## start and end of each trip
sten_stg_start <- aggregate(start~tn+tn2, data = f3, min)##start and end of each stage
sten_stg_end <- aggregate(end~tn+tn2, data = f3, max)
sten_stg <- merge(sten_stg_start, sten_stg_end, by = c("tn", "tn2"))

sten_stg$id <- unique(f3$ref)

##add start and end times for each stage to the dataframe
# start <- ISOdatetime(1970,1,1,0,0,0, tz="GMT")+sten_stg[,3][,1]
# end<- ISOdatetime(1970,1,1,0,0,0, tz="GMT")+sten_stg[,3][,2]
# sten <- data.frame(stage=sten_stg$tn, start=start, end=end)
sten <- sten_stg[ , c("id", "tn", "start", "end")]
names(sten) <- c("id", "stage", "start", "end")
#sten_stage <- rbind(sten_stage, sten)
tmp <- list(file = sten)
return(tmp)
  })
  return(stage_times)
}
