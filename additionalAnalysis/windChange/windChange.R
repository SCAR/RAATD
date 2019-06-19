# Change in wind

setwd("/perm_storage/home/shared/data_extra/env_change/")

library(raster)
library(raadtools)
library(RColorBrewer)

# Dates

# 1987-1998
d1 <- seq(as.POSIXct("1987-01-01"), as.POSIXct("1998-12-31"), 86400*1)

# 2007-2017
d2 <- seq(as.POSIXct("2007-01-01"), as.POSIXct("2017-12-31"), 86400*1)


# Get Wind
wind1 <- mean(readwind(d1, magonly = TRUE, lon180 = TRUE, latest = F), na.rm = T) # don't use xylim, crops edge
wind2 <- mean(readwind(d2, magonly = TRUE, lon180 = TRUE, latest = F), na.rm = T) # don't use xylim, crops edge

# Difference
dWind <- wind2 - wind1

# Crop and mask
landmask <- raster("/perm_storage/home/ryan/RAATD_01/RAATD/antarcticMask/mask.grd")
dWind <- resample(dWind, landmask)
dWind <- mask(dWind, landmask)

# Save the output
fldr <- paste0("/perm_storage/home/shared/data_extra/env_change/")
writeRaster(dWind, paste0(fldr, "wind.grd"), format = "raster", overwrite = T)

pal <- c(rev(RColorBrewer::brewer.pal(3, "Reds")), "#F7F7F7", RColorBrewer::brewer.pal(6, "Blues"))
plot(clamp(dWind, lower=-2, upper=2), col=pal, axes=FALSE, box=FALSE)
