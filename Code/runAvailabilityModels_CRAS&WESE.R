## Availability models for non-Central-Place-Foragers
## CRAS

# Ryan Reisinger
# April 2018

#--------------------------------------------------------
# Prepare

this.species <- "WESE"
this.stage <- "no-stage"

ice.thresh <- 15

too.big <- FALSE # Sub-sample large data?

library(raster)
library(raadtools)
library(dplyr)
library(scam)
library(ggplot2)

#-------------------
# Species and stage specific stuff:

setwd(paste0("~/RAATD_01/RAATD/", this.species))

sp <- this.species

# Load availability
load(paste0(this.species, "-avail-", this.stage, ".Rdata")) # Object should be called "usage", re-name if not (stage-specific example)

# Plot to check
ggplot(data = avail, aes(x = lon, y = lat, fill = avail)) + geom_raster() + coord_quickmap() +
  scale_fill_distiller(type = "seq", palette = "Blues")

# Keep only the required columns
avail <- avail[ , c("lon", "lat", "avail", "dist_deploy")]

#-------------------
# Expand the dataframe to represent one row per cell availability

require(splitstackshape)
isnot.avail <- avail[avail$avail == 0, ] # Not available cells
is.avail <- avail[avail$avail > 0, ] # Available cells
is.avail <- expandRows(dataset = is.avail, count = "avail", drop = FALSE) # Expand
is.avail$avail <- 1
rm(avail)

# If the model can't be fitted because the data are too large
if (too.big) {
  is.avail <- sample_frac(is.avail, size = 0.5)
  isnot.avail <- sample_frac(isnot.avail, size = 0.5)
}

avail <- rbind(is.avail, isnot.avail)
rm(is.avail, isnot.avail)

avail$is.avail <- "Y"
avail[avail$avail == 0, "is.avail"] <- "N"
avail$is.avail <- as.factor(avail$is.avail)

#-------------------
# Instead of using the the distance columns already in the data,
# here we extract a custom layer on which availability should be based

# First, define dates to use (same dates as for environmental data extractions)
# defined in 'extractEnvars_0X.R'

dates <- seq(as.POSIXct("2004-01-01"), as.POSIXct("2014-12-31"), 86400*7)

# Make an exmpty raster
lims <- c(-180, 180, -80, -40)
grd <- raster(extent(lims), resolution = c(0.1, 0.1))

# Next, extract and calculate the data of interest

#--------
# Ice distance

# land mask
bath <- readderivaadc("bathymetry", xylim = lims)
bath <- resample(bath, grd, method = "bilinear")
bath <- crop(bath, extent(grd))
landmask <- bath
landmask[landmask >= 0] <- NA

# Get ice
ice <- mean(readice(dates, setNA = T), na.rm = T)
rd3dummy <- grd
projection(grd) <- "+proj=longlat +datum=WGS84"
rd3dummy <- grd
rd3dummy[] <- 0
rd3points <-  rasterToPoints(rd3dummy, spatial = TRUE)
mice <- extract(ice, rd3points, method = "bilinear")
mice <- setValues(grd, mice)

# Apply the threshold
mice[mice < ice.thresh] <- NA

# Aggregate for faster calculation
mice <- aggregate(mice, 5)

# Distance
ice_distance <- distance(mice)

# Interpolate back to original resolution
ice_distance <- resample(ice_distance, grd, method = "bilinear")

# Mask land
ice_distance <- mask(ice_distance, landmask)

# Mask ice-sheets

# Clean up
rm(bath, ice, rd3dummy, rd3points)

# Save for use in av models later
writeRaster(ice_distance, paste0("/perm_storage/home/shared/distanceColony/distanceColony_", this.species, ".grd"),
            overwrite = TRUE)

#-------------------
# Extract
# wrap coordinatest to -180:+180
locs <- avail[ , c("lon", "lat")]
locs[locs$lon > 180, "lon"] <- locs[locs$lon > 180, "lon"] - 360

avail$dist_ice <- extract(ice_distance, locs)

#-------------------
# Fit model

scamMod <- scam(formula = is.avail ~ s(dist_ice, bs="mpd"),
                family = binomial,
                data = avail)

summary(scamMod)
plot(scamMod)
saveRDS(scamMod, paste0("fittedMods/", this.species, "_", this.stage, "_availSCAM.RDS"))

# Load model if neccessary
# scamMod <- readRDS(paste0("fittedMods/", this.species, "_", this.stage, "_availSCAM.RDS"))

#-------------------
# Predict the models

# Get data for prediction
dist <- as.data.frame(rasterToPoints(ice_distance))
names(dist) <- c("lon", "lat", "dist_ice") # names as in model input data

# Predict
pred <- predict.scam(scamMod, newdata = dist, type = "response")
dist$Pscam <- as.vector(pred)

# Plot to check
ggplot(data = dist, aes(x = lon, y = lat, fill = Pscam)) + geom_raster() + coord_quickmap() +
  scale_fill_distiller(type = "seq", palette = "Blues")

#-------------------
# Rasterize predictions
preds <- rasterize(x = dist[ c("lon", "lat")], y = ice_distance, field = dist$Pscam, fun = mean)

# Plot to check
plot(preds)

# Save the prediction raster for multiplication later
writeRaster(preds, file = paste0("/mnt/extra_storage/availabilityRasters/", this.species, "_", this.stage, "_availSCAM_rast.grd"),
            format = "raster", overwrite = T)
