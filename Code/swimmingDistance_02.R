### ---
## calculate swimming distance using connected grid with >8 neighbours
## similar to raster::gridDistance but the extra neighbours make a big difference to accuracy of results

## BR, March 2018, updated Sep 2018

setwd("~/RAATD_01/RAATD/")

this.species <- "WHCP"

## Threshold colony size?
threshold <- TRUE # yes or no
thr <- 0 # colony size threshold

library(igraph)
library(Matrix)
library(dplyr)
library(raadtools)
library(viridis)

## construct undirected graph representing area
lon_step <- 0.1; lat_step <- 0.1

## earlier BR definition
##lon <- seq(from=-180,to=180,by=lon_step)[-1] ## -1 so that 180 is not repeated at -180
##lat <- rev(seq(from=-80,to=-40,by=lat_step))

## make this lon/lat grid match the one RR is generating, but extending further north to accommodate DMSA colonies
grid_northern_limit <- -36
lon <- seq(from = -180 + lon_step/2, to = 180 - lon_step/2, by = lon_step)
lat <- rev(seq(from = -80 + lat_step/2, to = grid_northern_limit - lat_step/2, by = lat_step))
## matching raster template
raster_templ <- raster(ext = extent(c(-180, 180, -80, grid_northern_limit)), res = c(lon_step, lat_step), crs = CRS("+proj=longlat +ellps=WGS84"))

if (FALSE) {
  ## small test example
  lon_step <- 1
  lat_step <- 1
  lon <- 1:8
  lat <- rev(-65:-60)
}

ngrid <- length(lon)*length(lat)
ll <- as_tibble(expand.grid(lon,lat)) %>% dplyr::rename(lon=Var1,lat=Var2)

## need to remove vertices that fall on land
## figure these out now on the basis of unjittered positions
x0 <- readderivaadc("bathymetry")
## maybe use finer res bathy for this??

## resample to grid aligned with our lon lat grid
x0 <- resample(x0, raster_templ)
## coords in x0 should align with ours, and be in same order
chk <- nrow(coordinates(x0))==nrow(ll) && sum(abs(coordinates(x0)[,1]-ll[,1]))<1e-6 && sum(abs(coordinates(x0)[,2]-ll[,2]))<1e-6
if (!chk) stop("coordinate misalignment")

## discard values over land
to_strip <- values(x0)>0
## check plot
##with(ll[to_strip,],plot(lon,lat))
ll <- ll[!to_strip,]

cached_adjacency_matrix <- file.path("/perm_storage","home","shared","cached_stuff","swimming_distance",paste0("Aw_",lon_step,"_",lat_step,".rds"))
cached_adjacency_graph <- file.path("/perm_storage","home","shared","cached_stuff","swimming_distance",paste0("adjgraph_",lon_step,"_",lat_step,".rds"))

if (!file.exists(cached_adjacency_matrix)) {
  ## horrible code to build connected graph
  ## each cell in our grid gets connected to neighbouring cells

  ## figure out how many connections we expect in our graph
  ## expect that we have full 8-neighbours for all except min and max latitudes, and 5 neighbours for anything at max or min latitude
  expected_sum <- (8*length(lon)*(length(lat)-2)+5*2*length(lon))/2 ## /2 because we will have an undirected graph, so don't need to replicate edges
  ## add SSW, SSE neighbours: 2 neighbours but not for last 2 rows
  expected_sum <- expected_sum+(length(lat)-2)*length(lon)*2
  ## add EES, WWS neighbours: 2 neighbours but not for last row
  expected_sum <- expected_sum+(length(lat)-1)*length(lon)*2
  ## add SSSE, SSSW: 2 neighbours but not for last 3 rows
  expected_sum <- expected_sum+(length(lat)-3)*length(lon)*2
  ## add SSSEE, SSSWW: 2 neighbours but not for last 3 rows
  expected_sum <- expected_sum+(length(lat)-3)*length(lon)*2
  # EEES, WWWS: 2 but not for last row
  expected_sum <- expected_sum+(length(lat)-1)*length(lon)*2
  ## EEESS, WWWSS: 2 but not for last 2 rows
  expected_sum <- expected_sum+(length(lat)-2)*length(lon)*2

  ##edge_list <- matrix(NA_integer_,nrow=expected_sum,ncol=2)
  edge_from <- rep(NA_integer_,expected_sum)
  edge_to <- rep(NA_integer_,expected_sum)
  ptr <- 1

  ## rc to index, row-filled
  rc2ind <- function(r,c,sz=c(length(lat),length(lon))) c+sz[2]*(r-1)
  ind2rc <- function(ind,sz=c(length(lat),length(lon))) {
    temp <- (ind-1)/sz[2]+1
    matrix(c(floor(temp),(temp %% 1)*sz[2]+1),ncol=2,byrow=FALSE)
  }
  ind2c <- function(ind,sz=c(length(lat),length(lon))) ((ind-1) %% sz[2])+1
  ind2r <- function(ind,sz=c(length(lat),length(lon))) floor((ind-1)/sz[2]+1)

  ## eastern neighbour
  nbfrom <- seq_len(ngrid)
  nbto <- nbfrom+1
  chk <- ind2r(nbfrom)==ind2r(nbto) ## only those on the same row (same lat)
  edge_from[ptr:(ptr+sum(chk)-1)] <- nbfrom[chk]
  edge_to[ptr:(ptr+sum(chk)-1)] <- nbto[chk]
  ptr <- ptr+sum(chk)
  ## and wrapped from last col to first col
  nbfrom <- rc2ind(seq_len(length(lat)),length(lon))
  nbto <- rc2ind(seq_len(length(lat)),1)
  edge_from[ptr:(ptr+length(nbfrom)-1)] <- nbfrom
  edge_to[ptr:(ptr+length(nbfrom)-1)] <- nbto
  ptr <- ptr+length(nbfrom)
  ## southern neighbour (exclude last row)
  nbfrom <- seq_len(ngrid-length(lon))
  nbto <- nbfrom+length(lon)
  edge_from[ptr:(ptr+length(nbfrom)-1)] <- nbfrom
  edge_to[ptr:(ptr+length(nbfrom)-1)] <- nbto
  ptr <- ptr+length(nbfrom)
  ## SW neighbour (exclude last row and first column)
  nbfrom <- setdiff(seq_len(ngrid-length(lon)),rc2ind(seq_len(length(lat)-1),1))
  nbto <- nbfrom+length(lon)-1
  edge_from[ptr:(ptr+length(nbfrom)-1)] <- nbfrom
  edge_to[ptr:(ptr+length(nbfrom)-1)] <- nbto
  ptr <- ptr+length(nbfrom)
  ## first col
  nbfrom <- rc2ind(seq_len(length(lat)-1),1)
  nbto <- nbfrom+2*length(lon)-1
  edge_from[ptr:(ptr+length(nbfrom)-1)] <- nbfrom
  edge_to[ptr:(ptr+length(nbfrom)-1)] <- nbto
  ptr <- ptr+length(nbfrom)
  ## SE neighbour (exclude last row and last column)
  nbfrom <- setdiff(seq_len(ngrid-length(lon)),rc2ind(seq_len(length(lat)-1),length(lon)))
  nbto <- nbfrom+length(lon)+1
  edge_from[ptr:(ptr+length(nbfrom)-1)] <- nbfrom
  edge_to[ptr:(ptr+length(nbfrom)-1)] <- nbto
  ptr <- ptr+length(nbfrom)
  ## last col
  nbfrom <- rc2ind(seq_len(length(lat)-1),length(lon))
  nbto <- nbfrom+1
  edge_from[ptr:(ptr+length(nbfrom)-1)] <- nbfrom
  edge_to[ptr:(ptr+length(nbfrom)-1)] <- nbto
  ptr <- ptr+length(nbfrom)

  ## SSE neighbour (exclude last 2 rows and last column)
  nbfrom <- setdiff(seq_len(ngrid-2*length(lon)),rc2ind(seq_len(length(lat)),length(lon)))
  nbto <- nbfrom+2*length(lon)+1
  edge_from[ptr:(ptr+length(nbfrom)-1)] <- nbfrom
  edge_to[ptr:(ptr+length(nbfrom)-1)] <- nbto
  ptr <- ptr+length(nbfrom)
  ## SSE last col
  nbfrom <- rc2ind(seq_len(length(lat)-2),length(lon))
  nbto <- nbfrom+length(lon)+1
  edge_from[ptr:(ptr+length(nbfrom)-1)] <- nbfrom
  edge_to[ptr:(ptr+length(nbfrom)-1)] <- nbto
  ptr <- ptr+length(nbfrom)

  ## SSW neighbour (exclude last 2 rows and first column)
  nbfrom <- setdiff(seq_len(ngrid-2*length(lon)),rc2ind(seq_len(length(lat)),1))
  nbto <- nbfrom+2*length(lon)-1
  edge_from[ptr:(ptr+length(nbfrom)-1)] <- nbfrom
  edge_to[ptr:(ptr+length(nbfrom)-1)] <- nbto
  ptr <- ptr+length(nbfrom)
  ## SSW first col
  nbfrom <- rc2ind(seq_len(length(lat)-2),1)
  nbto <- nbfrom+3*length(lon)-1
  edge_from[ptr:(ptr+length(nbfrom)-1)] <- nbfrom
  edge_to[ptr:(ptr+length(nbfrom)-1)] <- nbto
  ptr <- ptr+length(nbfrom)

  ## EES neighbour (exclude last row and last 2 columns)
  nbfrom <- setdiff(seq_len(ngrid-length(lon)),rc2ind(seq_len(length(lat)),length(lon)))
  nbfrom <- setdiff(nbfrom,rc2ind(seq_len(length(lat)),length(lon)-1))
  nbto <- nbfrom+length(lon)+2
  edge_from[ptr:(ptr+length(nbfrom)-1)] <- nbfrom
  edge_to[ptr:(ptr+length(nbfrom)-1)] <- nbto
  ptr <- ptr+length(nbfrom)
  ## EES last 2 cols excluding last row element
  nbfrom <- c(rc2ind(seq_len(length(lat)-1),length(lon)),rc2ind(seq_len(length(lat)-1),length(lon)-1))
  nbto <- nbfrom+2
  edge_from[ptr:(ptr+length(nbfrom)-1)] <- nbfrom
  edge_to[ptr:(ptr+length(nbfrom)-1)] <- nbto
  ptr <- ptr+length(nbfrom)

  ## WWS neighbour (exclude last row and first 2 columns)
  nbfrom <- setdiff(seq_len(ngrid-length(lon)),rc2ind(seq_len(length(lat)),1))
  nbfrom <- setdiff(nbfrom,rc2ind(seq_len(length(lat)),2))
  nbto <- nbfrom+length(lon)-2
  edge_from[ptr:(ptr+length(nbfrom)-1)] <- nbfrom
  edge_to[ptr:(ptr+length(nbfrom)-1)] <- nbto
  ptr <- ptr+length(nbfrom)
  ## WWS first 2 cols excluding last row element
  nbfrom <- c(rc2ind(seq_len(length(lat)-1),1),rc2ind(seq_len(length(lat)-1),2))
  nbto <- nbfrom+2*length(lon)-2
  edge_from[ptr:(ptr+length(nbfrom)-1)] <- nbfrom
  edge_to[ptr:(ptr+length(nbfrom)-1)] <- nbto
  ptr <- ptr+length(nbfrom)

  ## SSSE neighbour (exclude last 3 rows and last column)
  nbfrom <- setdiff(seq_len(ngrid-3*length(lon)),rc2ind(seq_len(length(lat)),length(lon)))
  nbto <- nbfrom+3*length(lon)+1
  edge_from[ptr:(ptr+length(nbfrom)-1)] <- nbfrom
  edge_to[ptr:(ptr+length(nbfrom)-1)] <- nbto
  ptr <- ptr+length(nbfrom)
  ## SSSE last col
  nbfrom <- rc2ind(seq_len(length(lat)-3),length(lon))
  nbto <- nbfrom+2*length(lon)+1
  edge_from[ptr:(ptr+length(nbfrom)-1)] <- nbfrom
  edge_to[ptr:(ptr+length(nbfrom)-1)] <- nbto
  ptr <- ptr+length(nbfrom)

  ## SSSW neighbour (exclude last 3 rows and first column)
  nbfrom <- setdiff(seq_len(ngrid-3*length(lon)),rc2ind(seq_len(length(lat)),1))
  nbto <- nbfrom+3*length(lon)-1
  edge_from[ptr:(ptr+length(nbfrom)-1)] <- nbfrom
  edge_to[ptr:(ptr+length(nbfrom)-1)] <- nbto
  ptr <- ptr+length(nbfrom)
  ## SSSW first col
  nbfrom <- rc2ind(seq_len(length(lat)-3),1)
  nbto <- nbfrom+4*length(lon)-1
  edge_from[ptr:(ptr+length(nbfrom)-1)] <- nbfrom
  edge_to[ptr:(ptr+length(nbfrom)-1)] <- nbto
  ptr <- ptr+length(nbfrom)

  ## SSSEE neighbour (exclude last 3 rows and last 2 columns)
  nbfrom <- setdiff(seq_len(ngrid-3*length(lon)),rc2ind(seq_len(length(lat)),length(lon)))
  nbfrom <- setdiff(nbfrom,rc2ind(seq_len(length(lat)),length(lon)-1))
  nbto <- nbfrom+3*length(lon)+2
  edge_from[ptr:(ptr+length(nbfrom)-1)] <- nbfrom
  edge_to[ptr:(ptr+length(nbfrom)-1)] <- nbto
  ptr <- ptr+length(nbfrom)
  ## SSSEE last 2 cols
  nbfrom <- rc2ind(seq_len(length(lat)-3),length(lon))
  nbfrom <- c(nbfrom,rc2ind(seq_len(length(lat)-3),length(lon)-1))
  nbto <- nbfrom+2*length(lon)+2
  edge_from[ptr:(ptr+length(nbfrom)-1)] <- nbfrom
  edge_to[ptr:(ptr+length(nbfrom)-1)] <- nbto
  ptr <- ptr+length(nbfrom)

  ## SSSWW neighbour (exclude last 3 rows and first 2 columns)
  nbfrom <- setdiff(seq_len(ngrid-3*length(lon)),rc2ind(seq_len(length(lat)),1))
  nbfrom <- setdiff(nbfrom,rc2ind(seq_len(length(lat)),2))
  nbto <- nbfrom+3*length(lon)-2
  edge_from[ptr:(ptr+length(nbfrom)-1)] <- nbfrom
  edge_to[ptr:(ptr+length(nbfrom)-1)] <- nbto
  ptr <- ptr+length(nbfrom)
  ## SSSWW first 2 cols
  nbfrom <- rc2ind(seq_len(length(lat)-3),1)
  nbfrom <- c(nbfrom,rc2ind(seq_len(length(lat)-3),2))
  nbto <- nbfrom+4*length(lon)-2
  edge_from[ptr:(ptr+length(nbfrom)-1)] <- nbfrom
  edge_to[ptr:(ptr+length(nbfrom)-1)] <- nbto
  ptr <- ptr+length(nbfrom)

  ## EEES neighbour (exclude last row and last 3 columns)
  nbfrom <- setdiff(seq_len(ngrid-length(lon)),rc2ind(seq_len(length(lat)),length(lon)))
  nbfrom <- setdiff(nbfrom,rc2ind(seq_len(length(lat)),length(lon)-1))
  nbfrom <- setdiff(nbfrom,rc2ind(seq_len(length(lat)),length(lon)-2))
  nbto <- nbfrom+length(lon)+3
  edge_from[ptr:(ptr+length(nbfrom)-1)] <- nbfrom
  edge_to[ptr:(ptr+length(nbfrom)-1)] <- nbto
  ptr <- ptr+length(nbfrom)
  ## EEES last 3 cols excluding last row element
  nbfrom <- c(rc2ind(seq_len(length(lat)-1),length(lon)),rc2ind(seq_len(length(lat)-1),length(lon)-1),rc2ind(seq_len(length(lat)-1),length(lon)-2))
  nbto <- nbfrom+3
  edge_from[ptr:(ptr+length(nbfrom)-1)] <- nbfrom
  edge_to[ptr:(ptr+length(nbfrom)-1)] <- nbto
  ptr <- ptr+length(nbfrom)

  ## WWWS neighbour (exclude last row and first 3 columns)
  nbfrom <- setdiff(seq_len(ngrid-length(lon)),rc2ind(seq_len(length(lat)),1))
  nbfrom <- setdiff(nbfrom,rc2ind(seq_len(length(lat)),2))
  nbfrom <- setdiff(nbfrom,rc2ind(seq_len(length(lat)),3))
  nbto <- nbfrom+length(lon)-3
  edge_from[ptr:(ptr+length(nbfrom)-1)] <- nbfrom
  edge_to[ptr:(ptr+length(nbfrom)-1)] <- nbto
  ptr <- ptr+length(nbfrom)
  ## WWWS first 3 cols excluding last row element
  nbfrom <- c(rc2ind(seq_len(length(lat)-1),1),rc2ind(seq_len(length(lat)-1),2),rc2ind(seq_len(length(lat)-1),3))
  nbto <- nbfrom+2*length(lon)-3
  edge_from[ptr:(ptr+length(nbfrom)-1)] <- nbfrom
  edge_to[ptr:(ptr+length(nbfrom)-1)] <- nbto
  ptr <- ptr+length(nbfrom)

  ## EEESS neighbour (exclude last 2 rows and last 3 columns)
  nbfrom <- setdiff(seq_len(ngrid-2*length(lon)),rc2ind(seq_len(length(lat)),length(lon)))
  nbfrom <- setdiff(nbfrom,rc2ind(seq_len(length(lat)),length(lon)-1))
  nbfrom <- setdiff(nbfrom,rc2ind(seq_len(length(lat)),length(lon)-2))
  nbto <- nbfrom+2*length(lon)+3
  edge_from[ptr:(ptr+length(nbfrom)-1)] <- nbfrom
  edge_to[ptr:(ptr+length(nbfrom)-1)] <- nbto
  ptr <- ptr+length(nbfrom)
  ## EEESS last 3 cols excluding last 2 rows
  nbfrom <- c(rc2ind(seq_len(length(lat)-2),length(lon)),rc2ind(seq_len(length(lat)-2),length(lon)-1),rc2ind(seq_len(length(lat)-2),length(lon)-2))
  nbto <- nbfrom+3+length(lon)
  edge_from[ptr:(ptr+length(nbfrom)-1)] <- nbfrom
  edge_to[ptr:(ptr+length(nbfrom)-1)] <- nbto
  ptr <- ptr+length(nbfrom)

  ## WWWSS neighbour (exclude last 2 rows and first 3 columns)
  nbfrom <- setdiff(seq_len(ngrid-2*length(lon)),rc2ind(seq_len(length(lat)),1))
  nbfrom <- setdiff(nbfrom,rc2ind(seq_len(length(lat)),2))
  nbfrom <- setdiff(nbfrom,rc2ind(seq_len(length(lat)),3))
  nbto <- nbfrom+2*length(lon)-3
  edge_from[ptr:(ptr+length(nbfrom)-1)] <- nbfrom
  edge_to[ptr:(ptr+length(nbfrom)-1)] <- nbto
  ptr <- ptr+length(nbfrom)
  ## WWWSS first 3 cols excluding last 2 rows
  nbfrom <- c(rc2ind(seq_len(length(lat)-2),1),rc2ind(seq_len(length(lat)-2),2),rc2ind(seq_len(length(lat)-2),3))
  nbto <- nbfrom+3*length(lon)-3
  edge_from[ptr:(ptr+length(nbfrom)-1)] <- nbfrom
  edge_to[ptr:(ptr+length(nbfrom)-1)] <- nbto
  ptr <- ptr+length(nbfrom)

  A <- sparseMatrix(i=edge_from,j=edge_to,dims=c(ngrid,ngrid))
  if (sum(A)!=expected_sum) stop("A not as expected")

  ## clean up
  rm(nbfrom,nbto,edge_from,edge_to)

  if (FALSE) {
    ## this constructor will be very slow on a large matrix
    ## don't need g yet except for plot check below
    g <- graph_from_adjacency_matrix(as(as(A, "dgCMatrix"),"dgTMatrix"),mode="undirected")

    ## don't do this for anything other than small test grids
    plot(g,layout=layout_on_grid(g,width=length(lon)),edge.curved=FALSE)
    plot(g,layout=layout_on_grid(g,width=length(lon)),edge.curved=TRUE)
  }

  ## discard values over land
  A <- A[!to_strip,]
  A <- A[,!to_strip]
  ngrid <- nrow(ll)

  ## calculate edge weights as great-circle distance between jittered vertices
  ## haversine formula, ok for small distances
  idx <- which(A,arr.ind=TRUE)

  if (max(lon_step,lat_step)<=0.1) {
    ## whole lot at once
    W <- geosphere::distHaversine(ll[idx[,1],],ll[idx[,2],])
  } else {
    chsz <- 1e6
    templl1 <- ll[idx[,1],]
    templl2 <- ll[idx[,2],]
    W <- rep(NA,nrow(templl1))
    for (ci in seq(1,length(W)-chsz,by=chsz)) {
      chidx <- ci:min(ci+chsz-1,length(W))
      W[chidx] <- geosphere::distHaversine(templl1[chidx,],templl2[chidx,])
    }
  }
  Aw <- sparseMatrix(idx[,1],idx[,2],x=W)

  ## need Aw to be symmetric to construct the graph
  ## check that we don't have non-zero entries in both upper and lower triangular parts of matrix
  if (any(Aw>0 & t(Aw)>0)) stop("have conflicting entries in Aw")
  Aw <- Aw+t(Aw)

  ## build graph from Aw
  g <- graph_from_adjacency_matrix(Aw, mode = "undirected", weighted = TRUE)

  ## save these to re-use in later sessions
  saveRDS(Aw, file = cached_adjacency_matrix)
  saveRDS(g, file = cached_adjacency_graph)

} else {
  g <-  readRDS(cached_adjacency_graph)
}


## now calculate shortest paths from each starting point to all cells

## first check the connectedness of g: we should find that we have one very large connected subgraph, plus a tail of very small ones
if (vcount(g) != nrow(ll)) stop("graph size does not match ll grid")
chk <- components(g)
## this was with old BR grid def if (chk$no != 15 || sum(chk$csize < 20) != (chk$no - 1)) {
if (chk$no != 11 || sum(chk$csize < 35) != (chk$no - 1)) {
    stop("graph connectedness is unexpected")
} else {
    ## throw away the small disconnected components
    disconnected_cells <- chk$membership != which.max(chk$csize)
    if (sum(disconnected_cells) != 62) stop("removing different number of disconnected graph nodes than expected")
    ll <- ll[!disconnected_cells, ]
    g <- delete_vertices(g, which(disconnected_cells))
}
## check again
if (vcount(g) != nrow(ll)) stop("graph size does not match ll grid")
if (!is.connected(g)) stop("graph is not fully connected")

## start points to calculate distances from (i.e. colony locations)
cx <- read.table("/perm_storage/home/shared/github/raatd_modelling/dataLayers/distanceColony/colonies.csv",header=TRUE,sep=",",stringsAsFactors=FALSE)
cx <- cx[complete.cases(cx[ , c("lat", "lon")]), ]
cx <- cx %>% dplyr::filter(abbreviation == this.species)

## remove colonies smaller than the threshold defined above, if threshold is TRUE
## will also remove colonies with NA count
if (threshold) {
    chk <- is.na(cx$count) | cx$count <= thr
    if (any(chk)) warning(sprintf("removing %d colonies with count below the threshold of %d", sum(chk), thr))
    cx <- dplyr::filter(cx, count > thr)
}

## get distinct locations
from <- cx %>% dplyr::select(lon, lat) %>% distinct
## indices of these points in our ll grid
from_idx <- from %>% rowwise %>% dplyr::summarize(idx=which.min((ll[,1]-lon)^2+(ll[,2]-lat)^2)) %>% pull(idx) %>% unique
## nb for locations on land, this will choose nearest non-land point

chk <- sapply(adjacent_vertices(g, from_idx), length)
if (any(chk < 1)) stop("unreachable colonies")

if (any(cx$lat > max(lat))) stop("colony outside of grid extent")

D <- distances(g, v = from_idx) ## distance from each colony to each grid point, warning: SLOW!
if (any(is.infinite(D))) stop("infinite distance values")

D <- apply(D, 2, min) ## min distance to colony for each grid point

## cram it back into raster form
if (FALSE) {
    ## old BR raster
    Dr <- raster(ext=extent(c(-180+lon_step/2,180+lon_step/2,min(lat)-lat_step/2,max(lat)+lat_step/2)),res=c(lon_step,lat_step),crs=CRS("+proj=longlat +ellps=WGS84"))
    temp <- rep(NA,prod(dim(Dr)))
    temp[which(!to_strip)] <- D
    values(Dr) <- temp
    ## this won't have the -180 cells
    Dr <- extend(Dr,raster(ext=extent(c(-180-lon_step/2,180+lon_step/2,min(lat)-lat_step/2,max(lat)+lat_step/2)),res=c(lon_step,lat_step),crs=CRS("+proj=longlat +ellps=WGS84")))
    Dr[,1] <- Dr[,ncol(Dr)] ## fill -180 values with +180 values

    ## back to our 0.1-degree grid if necessary
    if (lon_step!=0.1 || lat_step!=0.1)
        Dr <- resample(Dr,raster(ext=extent(c(-180.05,180.05,-80.05,-39.95)),res=c(0.1,0.1),crs=CRS("+proj=longlat +ellps=WGS84")),method="ngb")
} else {
    ## matching RR rasters
    Dr <- raster_templ
    temp <- rep(NA, prod(dim(Dr)))
    temp[which(!to_strip)[!disconnected_cells]] <- D
    values(Dr) <- temp
    ## now crop back to RR extent
    Dr <- crop(Dr, extent(c(-180, 180, -80, -40)))
}

plot(Dr, col = rev(viridis(101)))
points(from$lon, from$lat, col = 1, bg = 1, pch = 21)

if (FALSE) {
  ## checks on method and results
  writeRaster(Dr,file=paste0("/perm_storage/home/shared/distanceColony/distance_to_colony-", this.species, "4.grd"),format="raster",overwrite=TRUE)

  ## compare to the native raster method
  dx <- raster(paste0("/perm_storage/home/shared/distanceColony/distance_to_colony-", this.species, "_1.grd"))
  dx2 <- resample(dx,Dr,method="ngb")
  hist(dx2[dx2>0 & dx2<100e3],1000)
  plot(Dr-dx2)
  points(from$lon,from$lat,col=1,bg=1,pch=21)

  ## HA! raster doesn't wrap properly

  stersouth='+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0' ## define the NSIDC projection and grid size for the Southern Hemisphere
  plot(projectRaster(Dr,crs=stersouth))

  plot(projectRaster(dx2,crs=stersouth))


  ## compare to exact distances for purposes of assessing quantization error
  tev <- geosphere::distHaversine(from[1,],ll)
  for (ci in 2:nrow(from)) tev <- pmin(tev,geosphere::distHaversine(from[ci,],ll))

  De <- raster(ext=extent(c(-180+lon_step/2,180+lon_step/2,min(lat)-lat_step/2,max(lat)+lat_step/2)),res=c(lon_step,lat_step),crs=CRS("+proj=longlat +ellps=WGS84"))
  temp <- rep(NA,prod(dim(De)))
  temp[which(!to_strip)] <- tev
  values(De) <- temp
  ## this won't have the -180 cells
  De <- extend(De,raster(ext=extent(c(-180-lon_step/2,180+lon_step/2,min(lat)-lat_step/2,max(lat)+lat_step/2)),res=c(lon_step,lat_step),crs=CRS("+proj=longlat +ellps=WGS84")))
  De[,1] <- De[,ncol(De)] ## fill -180 values with +180 values

  idx100 <- values(De)<100e3 & values(De)>0 & coordinates(De)[,1]>0
  plot(De[idx100],Dr[idx100]/De[idx100])
  plot(De[idx100],dx2[idx100]/De[idx100])

  plot(De[idx100],Drj[idx100]/De[idx100])

  mean(dx2[idx100]/De[idx100],na.rm=TRUE) ## 1.29
  range(dx2[idx100]/De[idx100],na.rm=TRUE) ## 0 to 60.80

  mean(Dr[idx100]/De[idx100],na.rm=TRUE) ## 1.03 with >8N
  range(Dr[idx100]/De[idx100],na.rm=TRUE) ## 0 to 2.35
}

writeRaster(Dr,file=paste0("/perm_storage/home/shared/distanceColony/distanceColony_", this.species, ".grd"),format="raster",overwrite=TRUE)