### ---
## colony weighting calculations
## uses much the same code as swimming distance
## calculate swimming distance using connected grid with >8 neighbours
## similar to raster::gridDistance but the extra neighbours make a big difference to accuracy of results

## BR, March 2018, Aug 2018

if (FALSE) {
this_species <- "BBAL"
##this_phase <- "chick-rearing"
##this_phase <- "incubation"
this_phase <- "post-breeding"

this_species <- "ANFS"
this_phase <- "breeding"
##this_phase <- "post-moult"
}
## masked out by if (FALSE) temporarily, this_species and this_phase defined in outer script that is sourcing this one

avail_fit_file <- file.path("/perm_storage/home/ryan/RAATD_01/RAATD", this_species, "fittedMods", sprintf("%s_%s_availSCAM.RDS", this_species, this_phase))
if (!file.exists(avail_fit_file)) stop("no fitted availability file appears to be available (looking for ", avail_fit_file, ")")

cat(sprintf("Calculating colony weighting layer for %s/%s\n", this_species, this_phase))

## Threshold colony size?
threshold <- FALSE # yes or no, don't do this for colony weighting
thr <- 49 # colony size threshold

library(scam)
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
ll <- ll_full <- as_tibble(expand.grid(lon,lat)) %>% dplyr::rename(lon=Var1,lat=Var2)

## need to remove vertices that fall on land
## figure these out now on the basis of unjittered positions
x0 <- readderivaadc("bathymetry")
## maybe use finer res bathy for this??

## resample to our lon lat grid
x0 <- resample(x0, raster_templ)
## coords in x0 should align with ours, and be in same order
chk <- nrow(coordinates(x0)) == nrow(ll) && sum(abs(coordinates(x0)[, 1]-ll[, 1])) < 1e-6 && sum(abs(coordinates(x0)[, 2]-ll[, 2])) < 1e-6
if (!chk) stop("coordinate misalignment")

## discard values over land
to_strip <- values(x0) > 0
## check plot
##with(ll[to_strip,],plot(lon,lat))
ll <- ll[!to_strip, ]

cache_dir <- file.path("/perm_storage", "home", "shared", "cached_stuff", "colony_weighting")
if (!dir.exists(cache_dir)) stop("cache directory does not exist")
cached_adjacency_matrix <- file.path(cache_dir, paste0("Aw_", lon_step, "_", lat_step, ".rds"))
cached_adjacency_graph <- file.path(cache_dir, paste0("adjgraph_", lon_step, "_", lat_step, ".rds"))

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

    A <- sparseMatrix(i = edge_from, j = edge_to, dims = c(ngrid, ngrid))
    if (sum(A) != expected_sum) stop("A not as expected")

    ## clean up
    rm(nbfrom, nbto, edge_from, edge_to)

    if (FALSE) {
        ## this constructor will be very slow on a large matrix
        ## don't need g yet except for plot check below
        g <- graph_from_adjacency_matrix(as(as(A, "dgCMatrix"),"dgTMatrix"),mode="undirected")

        ## don't do this for anything other than small test grids
        plot(g,layout=layout_on_grid(g,width=length(lon)),edge.curved=FALSE)
        plot(g,layout=layout_on_grid(g,width=length(lon)),edge.curved=TRUE)
    }


    ## discard values over land
    A <- A[!to_strip, ]
    A <- A[, !to_strip]
    ngrid <- nrow(ll)

    ## calculate edge weights as great-circle distance between jittered vertices
    ## haversine formula, ok for small distances
    idx <- which(A, arr.ind = TRUE)

    if (max(lon_step, lat_step) <= 0.1) {
        ## whole lot at once
        W <- geosphere::distHaversine(ll[idx[, 1], ], ll[idx[, 2], ])
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
    Aw <- sparseMatrix(idx[, 1], idx[, 2], x = W)

    ## need Aw to be symmetric to construct the graph
    ## check that we don't have non-zero entries in both upper and lower triangular parts of matrix
    if (any(Aw > 0 & t(Aw) > 0)) stop("have conflicting entries in Aw")
    Aw <- Aw + t(Aw)

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
cx <- cx %>% dplyr::filter(abbreviation == this_species)

chk <- is.na(cx$count)
if (any(chk)) {
    warning(sprintf("removing %d colony entries with missing count", sum(chk)))
    cx <- cx %>% dplyr::filter(!is.na(count))
}

## remove colonies smaller than the threshold defined above, if threshold is TRUE
if (threshold) {
    chk <- cx$count <= thr
    if (any(chk)) warning(sprintf("removing %d colonies with count below the threshold of %d", sum(chk), thr))
    cx <- dplyr::filter(cx, count > thr)
}

## any duplicates?
chk <- cx %>% count(lon, lat) %>% dplyr::filter(n > 1)
if (nrow(chk) > 1) stop("duplicate colony locations")

from <- cx %>% dplyr::select(lon, lat)
## indices of these points in our ll grid
from_idx <- from %>% rowwise %>% dplyr::summarize(idx = which.min((ll[, 1] - lon)^2 + (ll[, 2] - lat)^2)) %>% pull(idx)
## nb for locations on land, this will choose nearest non-land point

chk <- sapply(adjacent_vertices(g, from_idx), length)
if (any(chk < 1)) stop("unreachable colonies")

if (any(cx$lat > max(lat))) stop("colony outside of grid extent")

D <- distances(g, v = from_idx) ## distance from each colony to each grid point, warning: SLOW!
if (nrow(D) != nrow(cx)) stop("dimension mismatch") ## shouldn't happen, but, you know ...

## scale by all colony sizes, using a 1/distance weighting
scaled_counts <- cx$count/sum(cx$count) ## scale each count to proportion of global population
if (FALSE) {
    invw <- 1/D ## strict 1/distance weighting
} else {
    ## use fitted availability function
    av_fit <- readRDS(avail_fit_file)
    ## predict by rows (colonies)
    invw <- t(apply(D, 1, function(z) predict(av_fit, newdata = data.frame(dist_deploy = z), type = "response", block.size = -1)))
    ## or whole lot at once, but might fail for many colonies
    ##invw <- predict(av_fit, newdata = data.frame(dist_deploy = as.numeric(D)), type = "response", block.size = -1)
    ##invw <- matrix(invw, nrow = nrow(D))
}
llw <- apply(invw, 2, function(z) sum(z * scaled_counts))


## cram it back into raster form
if (FALSE) {
    ## old BR raster
    llwr <- raster(ext = extent(c(-180 + lon_step/2, 180 + lon_step/2, min(lat) - lat_step/2, max(lat) + lat_step/2)), res = c(lon_step, lat_step), crs = CRS("+proj=longlat +ellps=WGS84"))
    temp <- rep(NA, prod(dim(llwr)))
    temp[which(!to_strip)[!disconnected_cells]] <- llw
    values(llwr) <- temp

    ## this won't have the -180 cells
    llwr <- extend(llwr, raster(ext = extent(c(-180-lon_step/2, 180+lon_step/2, min(lat)-lat_step/2, max(lat)+lat_step/2)), res = c(lon_step,lat_step), crs = CRS("+proj=longlat +ellps=WGS84")))
    llwr[, 1] <- llwr[, ncol(llwr)] ## fill -180 values with +180 values

    ## back to our 0.1-degree grid if necessary
    if (lon_step!=0.1 || lat_step!=0.1)
        llwr <- resample(llwr,raster(ext=extent(c(-180.05,180.05,-80.05,-39.95)),res=c(0.1,0.1),crs=CRS("+proj=longlat +ellps=WGS84")),method="ngb")
} else {
    ## matching RR rasters
    llwr <- raster_templ
    temp <- rep(NA, prod(dim(llwr)))
    temp[which(!to_strip)[!disconnected_cells]] <- llw
    values(llwr) <- temp
    ## now crop back to RR extent
    llwr <- crop(llwr, extent(c(-180, 180, -80, -40)))
}

plot(llwr, col = viridis(101), main = sprintf("%s/%s", this_species, this_phase))
points(from$lon, from$lat, col = "white", bg = "red", pch = 21)

save_dir <- "/perm_storage/home/shared/colonyWeighting"
writeRaster(llwr, file = file.path(save_dir, sprintf("colonyWeights_%s_%s.grd", this_species, this_phase)), format = "raster", overwrite = TRUE)

if (FALSE) {
    ## check model predictions with and without this scaling
    ## this is not quite how the final one will look, because we'll apply the colony weighting prior to the percentile transformation
    bbx <- raster(file.path("/perm_storage/home/ryan/RAATD_01/RAATD/allSpRasters", sprintf("%s_%s_gbm_rast_trans.grd", this_species, this_phase)))
    plot(bbx, col = viridis(101))
    plot(bbx * llwr, col = viridis(101))
}

