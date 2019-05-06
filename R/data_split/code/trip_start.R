# Function to locate start and end time of each trip
# This is halfway between the end of one trip and start of another
# (defined using 'trip_split' function)

# Mark Hindell & Ryan Reisinger
# Last modified 2017-02-22

trip_start <- function(f, t.list){
  trip_times <- lapply(1:length(t.list), function(i){
    cat(sprintf('\nprocessing %d of %d\n',i,length(t.list))); flush.console()
    
ff <- f[f$ref == t.list[i],]
tsrt <- aggregate(gmt~ref+trip, data = ff, min)
names(tsrt) <- c("ref", "trip", "start")
tsrt$end <- aggregate(gmt~ref+trip, data = ff, max)$gmt

if (nrow(tsrt) > 1) {
tsrt$t2 <- c(tsrt[2:nrow(tsrt), "start"], NA) #get start of the next trip
tsrt$tsrt <- tsrt$end + (as.numeric(tsrt$t2) - as.numeric(tsrt$end))/2 #add the halfway time
tsrt$stfin <- c(tsrt$start[1], tsrt$tsrt[1:(nrow(tsrt) - 1)]) #add the halfway time as the start of the next trip
tsrt$endo <- c(tsrt$stfin[2:(nrow(tsrt))], tsrt$end[nrow(tsrt)]) #define end of trip as stfin of the next trip, except last trip which gets 'end'
stend <- tsrt[,c("ref", "trip", "stfin", "endo")]

}else{
 stend <- data.frame(ref = tsrt$ref, trip = tsrt$trip, stfin = tsrt$start, endo = tsrt$end) 
}


  }
)
  return(trip_times)
}