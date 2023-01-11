# Calculating uncertainty of predictions using bootstrap

# PART B (PART A must have been run already)

# Ryan Reisinger & Ben Raymond

# April 2019

library(raster)
library(scam)

library(gamsim)

library(pals)

# Libraries for parallel processing
library(iterators)
library(foreach)
library(parallel)
library(doParallel)

#--------------------------------------------------------
# Define species and stage

this.species <- "WESE"
this.stage <- "no-stage"

runInPar <- TRUE # Run in parallel?
ncores <- 10 # Number of cores, will use this number only if it is less than available cores

this.species.low <- tolower(this.species)

#--------------------------------------------------------
# Number of bootstraps

K <- 50

#-------------------
# Species and stage specific stuff:

setwd(paste0("~/RAATD_01/RAATD/", this.species))
sp <- this.species

#-------------------
# Data for prediction
# Load the raster values for prediction
load(paste0("./envdata/", this.species.low, "_envdata_", this.stage, ".Rdata"))
stk <- stack(envdata)
rm(envdata)
pred.dat <- rasterToPoints(stk)
rm(stk)

# Load the GBM predictions
out <- readRDS(paste0("../predictionsError/", this.species, "_", this.stage, "_gbm_out.RDS"))

#--------------------------------------------------------
## Availability model
#--------------------------------------------------------

# Load model
scamMod <- readRDS(paste0("fittedMods/", this.species, "_", this.stage, "_availSCAM.RDS"))

#-------------------
# Get K values from the fitted model using Ben's gamsim package

# Get data for prediction
distIce <- raster(paste0("/perm_storage/home/shared/distanceColony/distanceColony_", this.species, ".grd"))
dist <- as.data.frame(pred.dat[ , c("x", "y")])
dist$dist_ice <- raster::extract(distIce, dist)

## get K draws from the fitted scam object
## note that the enforce behaviour is applied to the data that comes back from get_draws, so 
##  newdata needs to be passed as a vector of increasing distances (so that our curve is
##  decreasing wrt distance)
## no need to pass all distances in our dist$dist_ice data.frame, just pass a representative 
##  range (take sample quantiles to get even coverage across our dist values)
qsamp <- seq(from = 0, to = 1, length.out = 1001)
temp <- data.frame(dist_ice = quantile(na.omit(dist$dist_ice), qsamp))
temp.scam <- gamsim::get_draws(object = scamMod, N = K, newdata = temp, enforce = "nonincreasing")
## so each column in temp.scam defines a response curve drawn from the fitted model
## apply each of these curves to our actual distances in dist$dist_ice via approx
## I do not deny that this is inelegant, even by my low, low standards - BR
pred.scam <- do.call(cbind, lapply(seq_len(K), function(z) plogis(approx(x = temp$dist_ice, y = temp.scam[, z], xout = dist$dist_ice)$y)))
## use plogis because get_draws returns values on the linear predictor scale



if (FALSE) {
  ## visual check
  dist$av <- pred.scam[, 40]
  plot(rasterFromXYZ(dist[, c("x", "y", "dist_ice")]), main = "dist_ice")
  plot(rasterFromXYZ(dist[, c("x", "y", "av")]), main = "sampled av")
}

# Save in case next step fails
saveRDS(pred.scam, file = paste0("../predictionsError/", this.species, "_", this.stage, "_scam_out.RDS"))
# pred.scam <- readRDS(paste0("../predictionsError/", this.species, "_", this.stage, "_scam_out.RDS"))

# Multiply the two matrices to get habitat selection given availability for each of K
pred_gbm.av <- out*pred.scam

# Importance conversion
suitability_to_percentiles <- function(x, N = 200) {
    ## expect that x is a raster containing predicted usage probability (suitability multiplied by availability)
    cell_areas <- area(x)
    vals <- raster::values(x * cell_areas) ## km^2
    total_area <- sum(cell_areas[!is.na(x)]) ## non-land
    cell_areas <- raster::values(cell_areas)
    smallish <- quantile(vals, 0.1, na.rm = TRUE) ## cumulative area increases rapidly for small values, so put more points at small values
    tst <- c(seq(0, smallish, length.out = N), seq(smallish, max(vals, na.rm = TRUE), length.out = N)[-1])
    ## calculate the percentage of area corresponding to each of these tst values
    s2p <- function(z) sum(cell_areas[which(vals <= z)]) / total_area * 100
    arp <- vapply(tst, s2p, FUN.VALUE = 1)
    ## and use that relationship to transform the whole raster
    values(x) <- approx(tst, arp, vals)$y
    x
}

# Peform the conversion for each column of output (corresponding to each bootstrap iteration),
# in parallel

# Clean up
crds <- as.data.frame(pred.dat[ , c("x", "y")])
rm(pred.dat)
rm(dist)
rm(scamMod)

removeTmpFiles(h = 0)

# Set up parallel processing if specificied above,
# otherwise run sequentially
if (runInPar) {
  ncr <- detectCores() - 1
  ncr <- min(ncores, ncr)
  cluster <- makeCluster(ncr, type = "FORK")
  registerDoParallel(cluster)
} else {
  registerDoSEQ()
}


hold <- foreach (i=1:K,
                .combine='stack',
                .inorder=FALSE,
                .packages = c("raster")) %dopar% {
  
  # gbm * availability
  d <- crds
  d$z <- pred_gbm.av[ , i]
  d <- rasterFromXYZ(d, crs = CRS("+proj=longlat +ellps=WGS84"))
  
  # availability
  v <- crds
  v$z <- pred.scam[ , i]
  v <- rasterFromXYZ(v, crs = CRS("+proj=longlat +ellps=WGS84"))
  
  # importance
  ras.trans <- suitability_to_percentiles(x = d)
  
  # Clean up
  rm(d, v)
  removeTmpFiles(h = 0)
  gc()
  
  # ras.trans <- rasterToPoints(ras.trans, na.rm = FALSE)
  ras.trans
}

#--------------------------------------------------------
# Calculate summaries on the output

# Get the outputs as dataframes
hold <- rasterToPoints(hold)
hold.data <- as.data.frame(hold[ , 3:(K+2)])
hold <- as.data.frame(hold[ , c("x", "y")])

# Calculate 5th and 95th percentile, and the predcition interval
hold$lower <- parRapply(cl = cluster, x = hold.data, FUN = quantile, probs = c(0.05))
hold$upper <- parRapply(cl = cluster, x = hold.data, FUN = quantile, probs = c(0.95))
hold$interval <- hold$upper - hold$lower

# Variance
hold$var <- parRapply(cl = cluster, x = hold.data, FUN = var)

# SD
hold$sd <- parRapply(cl = cluster, x = hold.data, FUN = sd)

# Mean (for interest...)
hold$mean <- parRapply(cl = cluster, x = hold.data, FUN = mean)

# Stop the cluster
if (runInPar) {
  stopCluster(cluster)
}

#--------------------------------------------------------
## Create a raster of the given error statistic

# Define CRS
wgs <- CRS("+proj=longlat +ellps=WGS84")

# Create raster
rst <- rasterFromXYZ(xyz = hold[, c("x", "y", "sd")], crs=wgs)

plot(rst, col = ocean.thermal(125),
main = paste0(this.species, " ", this.stage, "\n",
              "Standard deviation of habitat importance (", K, " bootstraps)"))

# Save the error raster
writeRaster(rst, file = paste0("../predictionsError/", this.species, "_", this.stage, "_gbm_error_rast.grd"),
            format = "raster", overwrite = T)
