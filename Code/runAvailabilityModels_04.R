# Availability models

# This code runs availability models on the ouput of the Marseille package (S. Wotherspoon)
# Script is species and stage-specific

# Ryan Reisinger & Ben Raymond

# last updated 2018-04-12

library(scam)
library(raster)
library(ggplot2)
library(dplyr)

#--------------------------------------------------------
# Prepare

this.species <- "WAAL"
this.stage <- "chick-rearing"

too.big <- FALSE # Sub-sample large data?

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
# Shape constrained additive model

#-------------------
# Fit model

scamMod <- scam(formula = is.avail ~ s(dist_deploy, bs="mpd"),
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
distCol <- raster(paste0("/perm_storage/home/shared/distanceColony/distanceColony_", this.species, ".grd"))
dist <- as.data.frame(rasterToPoints(distCol))
names(dist) <- c("lon", "lat", "dist_deploy") # names as in model input data
#dist$dist_deploy <- dist$dist_deploy * 1000 # convert to m, as used in the model input

# Predict
pred <- predict.scam(scamMod, newdata = dist, type = "response")
dist$Pscam <- as.vector(pred)

# Plot to check
ggplot(data = dist, aes(x = lon, y = lat, fill = Pscam)) + geom_raster() + coord_quickmap() +
  scale_fill_distiller(type = "seq", palette = "Blues")

#-------------------
# Rasterize predictions
preds <- rasterize(x = dist[ c("lon", "lat")], y = distCol, field = dist$Pscam, fun = mean)

# Plot to check
plot(preds)

# Save the prediction raster for multiplication later
writeRaster(preds, file = paste0("~/../../../mnt/extra_storage/availabilityRasters/", this.species, "_", this.stage, "_availSCAM_rast.grd"),
format = "raster", overwrite = T)
