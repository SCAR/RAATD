# This function does the track splitting

# It assumes an individual is at the colony when it is within threshold distance of
# the deployment site.

# Mark Hindell & Ryan Reisinger
# Last modified 2017-08-22

trip_split <- function(g, s.list, meta, allow.haulout = FALSE, haulouts = NULL){
  out_trip <- lapply(1:length(s.list), function(i){
    cat(sprintf('\nprocessing %d of %d\n',i,length(s.list))); flush.console()
    #Unit type and colony buffer - this depends on the unit type and is defined in the config file
    unit <- as.character(meta$device_type[meta$individual_id == s.list[i]])
    tbuff <- buff$tbuff[buff$type == unit]
    
    f <- fits$ssm[[i]]$predicted
    f <- f[,c("date", "lon", "lat")]
    f$ref <- s.list[i]
    names(f) <- c("gmt", "lon", "lat", "ref")
    f$lon[f$lon > 720] <- f$lon[f$lon > 720] - 720
    f$lon[f$lon > 360] <- f$lon[f$lon > 360] - 360
    f$lon[f$lon > 180] <- f$lon[f$lon > 180] - 360
    
    f$doy <- as.POSIXlt(f$gmt)$yday

    #calculate the distance of all locations from deployment site
    
    dep <- c(meta$deployment_decimal_longitude[meta$individual_id == s.list[i]],
             meta$deployment_decimal_latitude[meta$individual_id == s.list[i]])
    l <- matrix(c(f$lon, f$lat), ncol = 2)
    f$dnh <- distGeo(dep, l)
    
    ### dist to other haulouts
    
    if (allow.haulout) {
      haulouts.temp <- haulouts[haulouts$site == s.list[i], ] #get only specific sites
      if (nrow(haulouts.temp) > 0) {
      dis <- matrix(nrow = nrow(l), ncol = nrow(haulouts.temp) + 1)
      for (k in 1:nrow(haulouts.temp)) {
        dd <- haulouts.temp[k, c("lon", "lat")]
        dis[, k] <- distGeo(dd, l)
      }
      dis[,k+1] <- f$dnh
      f$dnh <- apply(X = dis, MARGIN = 1, FUN = min)
      }
    }
    ### end dist to other haulouts
    
    f$home <- 0
    f$home[f$dnh >= tbuff] <- 1
    f$trip <- ceiling(cumsum(abs(c(0, diff(f$home == 0))))/2) + 1 ##adds a number to each trip
    f$idx <- seq(1, nrow(f))
    #  trips <- aggregate(gmt~trip, f[f$dnh>=tbuff,], range)
    #  trips$start <- ISOdatetime(1970, 1, 1, 0, 0, 0, tz="GMT")+trips$gmt[,1]
    #  trips$end <- ISOdatetime(1970, 1, 1, 0, 0, 0, tz="GMT")+trips$gmt[,2]
    t.list <- names(with(f[f$home > 0,], table(trip))) ##exclude trips that don't go to sea
    trip_times <- function(f, t.list){
      out <- lapply(1:length(t.list), function(j){
        #cat(sprintf('\nprocessing %d of %d\n',j,length(t.list))); flush.console()
        t <- f[f$trip %in% t.list[j],]
        t1 <- t[t$home > 0,]
        st1 <- f$gmt[f$idx == (t1$idx[1] - 1)]
        st2 <- f$gmt[f$idx == (t1$idx[1])]
        
        dnh1 <- f$dnh[f$idx == (t1$idx[1] - 1 )]
        dnh2 <- f$dnh[f$idx == (t1$idx[1])]
        if (length(st1) != 0) {        
        stps <- seq(dnh1, dnh2, length.out = 100)
        tims <- seq(st1, st2, length.out = 100)
        cross.leave <- ISOdatetime(1970, 1, 1, 0, 0, 0, tz = "GMT") + approx(stps, tims,  xout = tbuff)$y
        
        lon1 <- f$lon[f$idx == (t1$idx[1] - 1)]
        lon2 <- f$lon[f$idx == (t1$idx[1])]
        
        lat1 <- f$lat[f$idx == (t1$idx[1] - 1)]
        lat2 <- f$lat[f$idx == (t1$idx[1])]
        
        lons <- seq(lon1, lon2, length.out = 100)
        lats <- seq(lat1, lat2, length.out = 100)
        cross.lon.leave <- approx(tims, lons, xout = cross.leave)$y
        cross.lat.leave <- approx(tims, lats, xout = cross.leave)$y
        }else{
          cross.lon.leave <- f$lon[t$idx[1]]
          cross.lat.leave <- f$lat[t$idx[1]]
          cross.leave <- f$gmt[t$idx[1]]
          
        }
        
#        plot(lons, lats, pch=19)
#        abline(h=cross.lat.leave, col="red")
#        abline(v=cross.lon.leave, col="blue")
        
        en1 <- f$gmt[f$idx == t1$idx[nrow(t1)]]
        en2 <- f$gmt[f$idx == (t1$idx[nrow(t1)] + 1)]
         
        dnh1 <- f$dnh[f$idx == t1$idx[nrow(t1)]]
        dnh2 <- f$dnh[f$idx == (t1$idx[nrow(t1)] + 1)]
        if (length(en2) != 0) {
        stps <- seq(dnh1, dnh2, length.out = 100)
        tims <- seq(en1, en2, length.out = 100)
        cross.return <- ISOdatetime(1970, 1, 1, 0, 0, 0, tz = "GMT") + approx(stps, tims,  xout = tbuff)$y
        
        lon1 <- f$lon[f$idx == t1$idx[nrow(t1)]]
        lon2 <- f$lon[f$idx == (t1$idx[nrow(t1)] + 1)]
        
        lat1 <- f$lat[f$idx == t1$idx[nrow(t1)]]
        lat2 <- f$lat[f$idx == (t1$idx[nrow(t1)] + 1)]
        
        lons <- seq(lon1, lon2, length.out = 100)
        lats <- seq(lat1, lat2, length.out = 100)
        cross.lon.return <- approx(tims, lons, xout = cross.return)$y
        cross.lat.return <- approx(tims, lats, xout = cross.return)$y
        } else {
        cross.lon.return <- f$lon[t$idx[nrow(t1)]]
        cross.lat.return <- f$lat[t$idx[nrow(t1)]]
        cross.return <- f$gmt[t$idx[nrow(t1)]]
        }
#        plot(lons, lats, pch=19)
#        abline(h=cross.lat.return, col="red")
#        abline(v=cross.lon.return, col="blue")
        
        ##data for the trip
        
        f2 <- f[f$idx >= t1$idx[1] & f$idx <= t1$idx[nrow(t1)], c("ref", "trip", "gmt", "lon", "lat","dnh", "idx")]
        lz <- data.frame(ref = f2$ref[1], trip =  f2$trip[1], gmt = cross.leave, lon = cross.lon.leave, lat = cross.lat.leave,  dnh = NA, idx = NA )
        rz <- data.frame(ref = f2$ref[1], trip =  f2$trip[1], gmt = cross.return, lon = cross.lon.return, lat = cross.lat.return, dnh = NA,idx = NA )
        
        f3  <- rbind(lz, f2, rz)        
        tmp <- list(file = f3)
      })
      out
    }
    
    if (length(t.list > 0)) {
      
    d <- trip_times(f, t.list)
    tmp <- lapply(d, function(x) x$file)
    trips <- do.call(rbind, tmp)
    
    
    ##plot the data for visual checking
    nf <- layout(matrix(c(1,1,1,1,1,1,1,1,2,2,2,2), 3, 4, byrow = TRUE))
    #layout.show(nf)
    
    #Plot the data
    print(paste(i, f$ref[1], sep = "_"))
    
    tit <- paste(f$ref[1], f$trip[i], sep = "_")
    plot(f$lon, f$lat, pch = 19, type = "b", main = tit)
    with(f[f$home >= 1,], points(lon, lat, pch = 19, col = "red"))
    with(f[f$home == 0,], points(lon, lat, pch = 19, col = "black"))
    st.mark <- aggregate(gmt~trip, data = trips, min)
    end.mark <- aggregate(gmt~trip, data = trips, max)
    plot(f$gmt, f$dnh, pch = 19, type = "b")
    abline(h = tbuff, col = "red")
    abline(v = st.mark$gmt, col = "blue")
    abline(v = end.mark$gmt, col = "green")
    #scan()
    tmp <- list(file = trips)
}
  })
  out_trip
}