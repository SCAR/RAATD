# Function to rasterize id by id

rasterizeSplit <- function(dat, grd, fun, field, background) {
  
  ids <- unique(dat$id)
  
  these.rasters <- stack()
  
  for (i in 1:length(ids)) {
    print(i)
    this.dat <- dat[dat$id == ids[i], ]
    
    if (is.character(field)) {
      this.field <- this.dat[ , field]
    } else {
      this.field <- rep(1, nrow(this.dat))
    }
    
    this.raster <- rasterize(this.dat[ , c("lon", "lat")],
                           grd,
                           fun = fun,
                           field = this.field,
                           background = background)
    these.rasters <- stack(these.rasters, this.raster)
  }
  calc(these.rasters, sum)
}

# Implementation in binTracks

## Function to bin tracks
binTracks <- function(trks) {
  
  trkcount <- function(x,...) length(unique(na.omit(x)))
  
  ## Remove non-fixed points that do not meet the mask
  trks <- trks[trks$fixed | isSea(trks$date,trks$lon,trks$lat),]
  ## Separate into actual and simulated tracks
  actl <- trks[trks$sim==0,]
  siml <- trks[trks$sim>0,]
  ## Compute mean deployment location
  deploy <- c(mean(wrapLon(trks$deploy_lon,xmin(grd))),mean(trks$deploy_lat))
  
  list(
    avail=rasterizeSplit(dat = trks,grd,
                    fun="count",
                    field = NA,
                    background=0),
    usage=rasterizeSplit(dat = actl,grd,
                    fun="count",
                    field = NA,
                    background=0),
    nactl=rasterizeSplit(dat = actl,grd,
                    field="track",
                    fun=trkcount,
                    background=0),
    nsiml=rasterizeSplit(dat = siml,grd,
                    field="sim",
                    fun=trkcount,
                    background=0),
    dist=distanceFromPoints(grd,deploy)
  )
}
