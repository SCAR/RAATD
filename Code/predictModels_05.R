## Predict the already-run models

library(sp)
library(raster)
library(caret)
library(caretEnsemble)
library(raster)
library(viridis)

#--------------------------------------------
# Species and stage specific stuff:

this.species <- "ADPE"
this.stage<- "chick-rearing"
multiMod <- FALSE

this.species.low <- tolower(this.species)

setwd(paste0("~/RAATD_01/RAATD/", this.species))

# Load the raster values for prediction
load(paste0("./envdata/", this.species.low, "_envdata_", this.stage, ".Rdata"))
stk <- stack(envdata)
rm(envdata)

dat <- rasterToPoints(stk)
# rm(stk) # keep for masking later, will not be neccessary to mask in final version
#--------------------------------------------

# If combining multiple small models, get them each, predict, and take the mean,
# otherwise, load the single model and predict

if (multiMod) {
  
  foo <- data.frame("bar" = 1:nrow(dat))
  
  for (i in 1:20) {
    # Load the saved model
    gbm <- readRDS(paste0("fittedMods/", this.species, "_", this.stage, "_gbm_", i, ".RDS"))
    
    #-------
    ## Predict each model
    
    # gbm
    pred.gbm <- predict.train(gbm, newdata = dat, type = "prob")
    pred.gbm <- pred.gbm$Y
    
    foo <- cbind(foo, pred.gbm)
    
    print(paste0("Done: ", i))
  }
  
  foo$bar <- NULL
  pred.gbm <- rowMeans(foo[ , -1], na.rm = T)
  rm(foo)
  
} else {
  # Load the saved model
  gbm <- readRDS(paste0("fittedMods/", this.species, "_", this.stage, "_gbm.RDS"))
  
  #-------
  ## Predict each model
  
  # gbm
  pred.gbm <- predict.train(gbm, newdata = dat, type = "prob")
  pred.gbm <- pred.gbm$Y
}

#-------
## Combine all the model predictions into a single dataframe
# (here only gbm)
grid <- cbind.data.frame(dat, gbm = pred.gbm)

#-------
# Get availability predictions and extract them
if (this.species != "HUWH") {
  av <- raster( paste0("~/../../../mnt/extra_storage/availabilityRasters/", this.species, "_", this.stage, "_availSCAM_rast.grd"))
  grid$av <- raster::extract(av, grid[ , c("x", "y")])
} else {
  grid$av <- 1
}

#-------
# Multiply
grid$gbm.av <- grid$gbm * grid$av

# Save the dataframe
saveRDS(grid, paste0("modPreds/predFrame_", this.species, "_", this.stage, ".RDS"))

#-------
## Create a raster of the given model predictions

# Get boundaries from the config file
source(paste0("/perm_storage/home/shared/github/raatd_data/data_split/config/", this.species, "config.r"))

# Define CRS
wgs <- CRS("+proj=longlat +ellps=WGS84")

# Define limits
lims <- c(-180, 180, -80, -40)

# Create empty raster
rast <- raster(ext = extent(lims), res = c(0.1,0.1), crs = wgs)

#-------
# Rasterize raw predictions on new raster
rast.ens <- rasterize(x = dat[ ,c("x", "y")],
                               y = rast,
                               field = grid[ ,"gbm"],
                               fun = mean)

# Rasterize availability
rast.av <- rasterize(x = dat[ ,c("x", "y")],
                      y = rast,
                      field = grid[ ,"av"],
                      fun = mean)
  

# Rasterize raw predictions given availability on new raster
rast.gbm.av <- rasterize(x = dat[ ,c("x", "y")],
                      y = rast,
                      field = grid[ ,"gbm.av"],
                      fun = mean)

# Mask
mask <- raster("../antarcticMask/mask.grd")
rast.ens <- raster::mask(x = rast.ens, mask = mask)
rast.av<- raster::mask(x = rast.av, mask = mask)
rast.gbm.av <- raster::mask(x = rast.gbm.av, mask = mask)

# Plot to check
plot(rast.ens, main = paste0(this.species, " - ", this.stage, " - PREFERENCE"), col = viridis(125))
plot(rast.av, main = paste0(this.species, " - ", this.stage, " - AVAILABILITY"), col = viridis(125))
plot(rast.gbm.av, main = paste0(this.species, " - ", this.stage, " - PREFERENCE|AVAILABILITY"), col = viridis(125))

# Save the raster
writeRaster(rast.ens, file = paste0("modPreds/", this.species, "_", this.stage, "_gbm_rast.grd"), format = "raster", overwrite = T)
writeRaster(rast.av, file = paste0("modPreds/", this.species, "_", this.stage, "_av.grd"), format = "raster", overwrite = T)
writeRaster(rast.gbm.av, file = paste0("modPreds/", this.species, "_", this.stage, "_gbm-av_rast.grd"), format = "raster", overwrite = T)

#-------
# Convert to percentile
#total_area <- sum(cell_areas[!is.na(rast.gbm.av)])

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

ras.trans <- suitability_to_percentiles(x = rast.gbm.av)
plot(ras.trans, main = paste0(this.species, " - ", this.stage, " - TRANSFORMED"), col = viridis(125))
writeRaster(ras.trans, file = paste0("modPreds/", this.species, "_", this.stage, "_gbm_rast_trans.grd"), format = "raster", overwrite = T)


#-------
# Colony weighted

# Get colony weighted availability
if (this.species %in% c("CRAS", "HUWH", "WESE")) {
  col <- rast.av } else {
    col <- raster(paste0("/perm_storage/home/shared/colonyWeighting/colonyWeights_", this.species, "_", this.stage, ".grd"))
  }

# Multiply them, percentile transform and save the output
wHab <- rast.ens*col
wHab <- suitability_to_percentiles(x = wHab)
writeRaster(wHab, paste0("modPreds/", this.species, "_", this.stage, "_gbmcw_rast_trans.grd"), overwrite = T)

# Add values to the prediction dataframe
# TODO, but actually not neccessary

# Compare with non colony-weighted

par(mfrow = c(2, 1))
plot(ras.trans, main = "Not weighted by colony size", col = viridis(125))
plot(wHab, main = "Weighted by colony size", col = viridis(125))
