## Fast run-through to produce quick discussion results,
## based on species groups

library(rgdal)
library(raster)
library(sp)
library(ggplot2)
library(viridis)
library(sf)

setwd("~/RAATD_01/RAATD/additionalAnalysis/fastRun/")

# ---------------------------
# Specific stuff
group.1 <- c("ADPE", "ANPE", "CRAS", "EMPE", "HUWH", "SOES", "WESE")
group.2 <- c("ANFS", "BBAL", "DMSA", "GHAL", "HUWH", "KIPE", "LMSA", "MAPE.ROPE", "SOES", "WAAL", "WHCP")

species.groups <- list(group.1, group.2)

# ---------------------------
## Source functions

# Function to calculate AES
source("~/RAATD_01/RAATD/Code/function_aesPoly.R")

# Function to combine different species groups
source("~/RAATD_01/RAATD/Code/function_meanR.R")

# ---------------------------
## Mapping stuff

# Get Antarctic mask
landmask <- raster("../../antarcticMask/mask.grd")

# Projections
projW <- "+proj=longlat +datum=WGS84 +no_defs"
projP <- "+proj=laea +lat_0=-90  +lon_0=0 +datum=WGS84 +ellps=WGS84 +no_defs"

# Custom theme
theme_ryan <- function () { 
  theme_bw(base_size=9, base_family="") %+replace% 
    theme(
      axis.text = element_text(colour = "black"),
      axis.ticks = element_line(colour = "black"),
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA),
      axis.line = element_line(colour = "black")
    )
}

# ---------------------------
## MPAs
fldr <- paste0("/perm_storage/home/shared/data_extra/mpa")

existing_mpa <- readOGR(dsn = fldr, layer = "existing_mpa")
proposed_mpa <- readOGR(dsn = fldr, layer = "proposed_mpa")

# Rasterize

# Create reference grid
lims <- c(-180, 180, -80, -40)
grd <- raster(extent(lims), resolution = c(0.1, 0.1))

# Get area
a <- area(grd)

existing_rast <- rasterize(existing_mpa, grd)
existing_rast[!is.na(existing_rast)] <- 1
existing_rast[is.na(existing_rast)] <- 0
existing_rast <- mask(existing_rast, landmask)

proposed_rast <- rasterize(proposed_mpa, grd)
proposed_rast[!is.na(proposed_rast)] <- 2
proposed_rast[is.na(proposed_rast)] <- 0
proposed_rast <- mask(proposed_rast, landmask)

# ---------------------------
## Generate figures with AES

# ---------------------------
## Function to plot cumulative area

plotGraph <- function(this.data) {
  
  d1 <- this.data
  
  this.q75 <- quantile(d1$MEAN, probs = c(0.75), na.rm = T)
  this.q90 <- quantile(d1$MEAN, probs = c(0.90), na.rm = T)
  
  d1.Current <- subset(d1, MPA == 1)
  d1.Proposed <- subset(d1, MPA == 1 | MPA == 2)
  d1.Outside <- subset(d1, MPA == 0 | MPA == 1 | MPA == 2)
  
  d1.Current <- d1.Current[order(d1.Current$MEAN, decreasing = T), ]
  d1.Proposed <- d1.Proposed[order(d1.Proposed$MEAN, decreasing = T), ]
  d1.Outside <- d1.Outside[order(d1.Outside$MEAN, decreasing = T), ]
  
  d1.Current$cumArea <- cumsum(d1.Current$area)
  d1.Proposed$cumArea <- cumsum(d1.Proposed$area)
  d1.Outside$cumArea <- cumsum(d1.Outside$area)
  
  d1.Current$What <- "C.Current"
  d1.Proposed$What <- "B.Proposed"
  d1.Outside$What <- "A.Outside"
  
  # Calculate proportions
  
  d1.All <- rbind(d1.Outside, d1.Proposed, d1.Current)
  
  hold <- d1.All[d1.All$MEAN >= this.q75, ]
  q75.outside <- max(hold[hold$What == "A.Outside", "cumArea"], na.rm = T)
  q75.proposed <- max(hold[hold$What == "B.Proposed", "cumArea"], na.rm = T)
  q75.current <- max(hold[hold$What == "C.Current", "cumArea"], na.rm = T)
  
  hold <- d1.All[d1.All$MEAN >= this.q90, ]
  q90.outside <- max(hold[hold$What == "A.Outside", "cumArea"], na.rm = T)
  q90.proposed <- max(hold[hold$What == "B.Proposed", "cumArea"], na.rm = T)
  q90.current <- max(hold[hold$What == "C.Current", "cumArea"], na.rm = T)
  
  hold <- d1.All
  q100.outside <- max(hold[hold$What == "A.Outside", "cumArea"], na.rm = T)
  q100.proposed <- max(hold[hold$What == "B.Proposed", "cumArea"], na.rm = T)
  q100.current <- max(hold[hold$What == "C.Current", "cumArea"], na.rm = T)
  
  cat(paste0("For 90th percentile:", "\n",
             "Designated = ", (q90.current/q90.outside)*100, " % of area", "\n",
             "Designated + Proposed = ", (q90.proposed/q90.outside)*100, " % of area", "\n"))
  
  cat(paste0("For 75th percentile:", "\n",
             "Designated = ", (q75.current/q75.outside)*100, " % of area", "\n",
             "Designated + Proposed = ", (q75.proposed/q75.outside)*100, " % of area", "\n"))
  
  cat(paste0("For 0th percentile:", "\n",
             "Designated = ", (q75.current/q100.outside)*100, " % of area", "\n",
             "Designated + Proposed = ", (q100.proposed/q100.outside)*100, " % of area", "\n"))
  
  # d1.All$What <- as.factor(d1.All$What)
  
  #d1.All$What=factor(d1.All$What , levels = c("Outside", "Proposed", "Current"))
  
  p <- ggplot(data = d1.Outside, aes(x = cumArea, y = MEAN)) +
    geom_area(fill = "#BBBBBB") +
    geom_area(data = d1.Proposed, aes(x = cumArea, y = MEAN), fill = "#EE3377") +
    geom_area(data = d1.Current, aes(x = cumArea, y = MEAN), fill = "#EE7733") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "Cumulative area (km2)", y = "Mean habitat importance") +
    geom_hline(yintercept = this.q75) +
    geom_hline(yintercept = this.q90) +
    theme_ryan()
  
  return(p)
  
}

# ---------------------------
# a. Weighted
datW <- readRDS("../../meanPredictions/griddedScoresEnvarsTransformedColonyWeighted.RDS")
datW <- meanR(datW, species.groups)

datW$area <- extract(a, datW[ , c("x", "y")])

datW$existing <- extract(existing_rast, datW[ , c("x", "y")])
datW$proposed <- extract(proposed_rast, datW[ , c("x", "y")])
datW$MPA <- datW$existing + datW$proposed

p <- plotGraph(this.data = datW)

# png("mpaPropWeighted.png", width = 600, height = 600, units = "px", res = 150)
# print(p)
# dev.off()

# ---------------------------
# a. Non-weighted
datNW <- readRDS("../../meanPredictions/griddedScoresEnvarsTransformed.RDS")
datNW <- meanR(datNW, species.groups)

datNW$area <- extract(a, datNW[ , c("x", "y")])

datNW$existing <- extract(existing_rast, datNW[ , c("x", "y")])
datNW$proposed <- extract(proposed_rast, datNW[ , c("x", "y")])
datNW$MPA <- datNW$existing + datNW$proposed

p <- plotGraph(this.data = datNW)

# png("mpaPropUnweighted.png", width = 600, height = 600, units = "px", res = 150)
# print(p)
# dev.off()

# ---------------------------
## Generate maps with AES

prj <- "+proj=laea +lat_0=-90  +lon_0=0 +datum=WGS84 +ellps=WGS84 +no_defs"

# Project MPAs
existing_mpa_proj <- spTransform(existing_mpa, CRSobj = prj)
proposed_mpa_proj <- spTransform(proposed_mpa, CRSobj = prj)

# Prepare world map
# Create polar circle
library(dplyr)
library(rgeos)
library(geosphere)

data(countriesLow, package = "rworldmap")
polar.crop <- SpatialPoints(cbind(0, 0), proj4string = CRS(prj)) %>%
  gBuffer(., quadsegs = 1000, width = distGeo(cbind(0, -90), cbind(0, -41.5)))

## use a zero width buffer to ensure there's no topology problems (long story)
wrld <- spTransform(countriesLow, prj) %>%
  gBuffer(., width = 0, quadsegs = 1000) %>%
  crop(., polar.crop)

## Crop the MPAs
existing_mpa_proj <-gBuffer(existing_mpa_proj, width = 0, quadsegs = 1000) %>%
  crop(., polar.crop)
proposed_mpa_proj <-gBuffer(proposed_mpa_proj, width = 0, quadsegs = 1000) %>%
  crop(., polar.crop)

# Colour
trans.none <- rgb(0, 0, 0, alpha = 255/255, maxColorValue = 255)

trans.red <- rgb(231, 51, 119, alpha = 255/2, maxColorValue = 255)
trans.orange <- rgb(238, 119, 51, alpha = 255/2, maxColorValue = 255)

# ---------------------------
# a. Weighted

# Get the mean habitat importance
d <- rasterFromXYZ(datW[ , c("x", "y", "MEAN")])

# Get and project the raster
projection(d) <- "+proj=longlat +datum=WGS84"
d <- projectRaster(d, crs = prj)

# Polygons
poly75 <- readRDS("polyW75.RDS")
poly90 <- readRDS("polyW90.RDS")
poly75 <- st_transform(poly75, projP)
poly90 <- st_transform(poly90, projP)

# Plot
tiff("mpaMapWeighted.tiff", height = 8, width = 8, units = "in", res = 150) 
plot(d, axes = FALSE, bty = "n", box = FALSE, col = viridis(125),
     horizontal = TRUE,
     legend.shrink = 0.7,
     legend.args = list(text = 'Overall habitat importance'))
plot(wrld,
     col = "darkgrey",
     border = FALSE,
     add = TRUE)
plot(existing_mpa_proj, add = T, col = trans.orange, border = "#EE7733")
plot(proposed_mpa_proj, add = T, col = trans.red, border = "#EE3377")
#plot(poly75, add = T, col = trans.none, border = "grey75")
plot(poly90, add = T, col = trans.none, border = "grey100")
lines(polar.crop, lwd = 2, add = T)
dev.off()

# ---------------------------
# b. Unweighted

# Get the mean habitat importance
d <- rasterFromXYZ(datNW[ , c("x", "y", "MEAN")])

# Get and project the raster
projection(d) <- "+proj=longlat +datum=WGS84"
d <- projectRaster(d, crs = prj)

# Polygons
poly75 <- readRDS("polyNW75.RDS")
poly90 <- readRDS("polyNW90.RDS")
poly75 <- st_transform(poly75, projP)
poly90 <- st_transform(poly90, projP)

# Plot
tiff("mpaMapUnweighted.tiff", height = 8, width = 8, units = "in", res = 150) 
plot(d, axes = FALSE, bty = "n", box = FALSE, col = viridis(125),
     horizontal = TRUE,
     legend.shrink = 0.7,
     legend.args = list(text = 'Overall habitat importance'))
plot(wrld,
     col = "darkgrey",
     border = FALSE,
     add = TRUE)
#plot(poly75, add = T, col = trans.none, border = "grey75")
plot(poly90, add = T, col = trans.none, border = "grey100")
plot(existing_mpa_proj, add = T, col = trans.orange, border = "#EE7733")
plot(proposed_mpa_proj, add = T, col = trans.red, border = "#EE3377")
lines(polar.crop, lwd = 2, add = T)
dev.off()