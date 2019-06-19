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

## EEZ
eez <- readOGR(dsn = fldr, layer = "eez")

## CCAMLR
ccamlr <- readOGR(dsn = fldr, layer = "ccamlr")

# Rasterize

# Create reference grid
lims <- c(-180, 180, -80, -40)
grd <- raster(extent(lims), resolution = c(0.1, 0.1))
crs(grd) <- "+proj=longlat +datum=WGS84 +no_defs"

# Get area
a <- area(grd)
a <- mask(a, landmask)

existing_rast <- rasterize(existing_mpa, grd)
existing_rast[!is.na(existing_rast)] <- 1
existing_rast[is.na(existing_rast)] <- 0
existing_rast <- mask(existing_rast, landmask)

proposed_rast <- rasterize(proposed_mpa, grd)
proposed_rast[!is.na(proposed_rast)] <- 2
proposed_rast[is.na(proposed_rast)] <- 0
proposed_rast <- mask(proposed_rast, landmask)

eez_rast <- rasterize(eez, grd)
eez_rast[!is.na(eez_rast)] <- 1
eez_rast[is.na(eez_rast)] <- 0
eez_rast <- mask(eez_rast, landmask)


# CCAMLR polygons are not closed after projecting
# to WGS184,
# so use an alternative method,
# rasterizing a projected layer first

library(CCAMLRGIS)
ra <- load_ASDs(format  = "RDATA")
ra_grd <- raster(extent(ra), resolution = c(1000, 1000))
ra_rast <- rasterize(ra, ra_grd)
ccamlr_rast <- projectRaster(ra_rast, grd)

ccamlr_rast[!is.na(ccamlr_rast)] <- 1
ccamlr_rast[is.na(ccamlr_rast)] <- 0
ccamlr_rast <- mask(ccamlr_rast, landmask)

# only ccamlr outside eezs
ccamlr_rast <- ccamlr_rast-eez_rast
ccamlr_rast[ccamlr_rast < 0] <- 0

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
             "Designated = ", (q100.current/q100.outside)*100, " % of area", "\n",
             "Designated + Proposed = ", (q100.proposed/q100.outside)*100, " % of area", "\n"))
  
  hold <- data.frame("Where" = rep(c("B.AES", "A.Total"), 3),
                     "What" = c("C.Current.MPA", "C.Current.MPA", "B.Proposed.MPA", "B.Proposed.MPA", "A.Outside.MPA", "A.Outside.MPA"),
                     "Area" = c(q90.current, q100.current, q90.proposed-(q90.current), q100.proposed-(q100.current), q90.outside-(q90.current+q90.proposed), q100.outside-(q100.current+q100.proposed)))
  
  hold$Area <- hold$Area/1000000
  
  p <- ggplot(data = hold, aes(y = Area, x = Where, fill = What)) +
    geom_col(position = "stack") +
    coord_flip() +
    scale_y_continuous(expand = c(0, 0)) +
    labs(y = "Cumulative area (million km2)", x = "") +
    scale_fill_manual(values = c("#0077BB", "#EE3377", "#EE7733")) +
    theme_ryan()
    
  return(p)
  
}


# ---------------------------
## Generate figures with Jurisdiction

# ---------------------------
## Function to plot jurisdiction areas

plotJurisdiction <- function(this.data) {
  
  d1 <- this.data
  
  this.q90 <- quantile(d1$MEAN, probs = c(0.90), na.rm = T)
  
  d1.ccamlr <- subset(d1, ccamlr == 1)
  d1.eez <- subset(d1, eez == 1)
  d1.highseas <- subset(d1, ccamlr == 0 & eez == 0)
  
  d1.ccamlr <- d1.ccamlr[order(d1.ccamlr$MEAN, decreasing = T), ]
  d1.eez <- d1.eez[order(d1.eez$MEAN, decreasing = T), ]
  d1.highseas <- d1.highseas[order(d1.highseas$MEAN, decreasing = T), ]
  
  d1.ccamlr$cumArea <- cumsum(d1.ccamlr$area)
  d1.eez$cumArea <- cumsum(d1.eez$area)
  d1.highseas$cumArea <- cumsum(d1.highseas$area)
  
  d1.ccamlr$What <- "C.ccamlr"
  d1.eez$What <- "B.eez"
  d1.highseas$What <- "A.highseas"
  
  # Calculate proportions
  
  d1.All <- rbind(d1.highseas, d1.eez, d1.ccamlr)
  
  hold <- d1.All[d1.All$MEAN >= this.q90, ]
  q90.highseas <- max(hold[hold$What == "A.highseas", "cumArea"], na.rm = T)
  q90.eez <- max(hold[hold$What == "B.eez", "cumArea"], na.rm = T)
  q90.ccamlr <- max(hold[hold$What == "C.ccamlr", "cumArea"], na.rm = T)
  
  hold <- d1.All
  q100.highseas <- max(hold[hold$What == "A.highseas", "cumArea"], na.rm = T)
  q100.eez <- max(hold[hold$What == "B.eez", "cumArea"], na.rm = T)
  q100.ccamlr <- max(hold[hold$What == "C.ccamlr", "cumArea"], na.rm = T)
  
  cat(paste0("Proportions:", "\n",
             "CCAMLR (Excluding EEZ in CCAMLR area), AES = ", (q90.ccamlr/q100.ccamlr)*100, " % of area", "\n",
             "EEZ (All), AES = ", (q90.eez/q100.eez)*100, " % of area", "\n",
             "High Seas (Outside CCAMLR area), AES = ", (q90.highseas/q100.highseas)*100, " % of area", "\n"))
  
  hold <- data.frame("Where" = rep(c("B.AES", "A.Total"), 3),
                     "What" = c("B.ccamlr", "B.ccamlr", "C.eez", "C.eez", "A.highseas", "A.highseas"),
                     "Area" = c(q90.ccamlr, q100.ccamlr-q90.ccamlr,
                                q90.eez, q100.eez-q90.eez,
                                q90.highseas, q100.highseas-q90.highseas))
  
  hold$Area <- hold$Area/1000000
  
  p <- ggplot(data = hold, aes(y = Area, x = What, fill = Where)) +
    geom_col(position = "stack") +
    coord_flip() +
    scale_y_continuous(expand = c(0, 0)) +
    labs(y = "Cumulative area (million km2)", x = "") +
    scale_fill_manual(values = viridis(2)) +
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

pdf("mpaPropWeighted.pdf", width = 4.5, height = 2)
print(p)
dev.off()

# ---------------------------
# b. Non-weighted
datNW <- readRDS("../../meanPredictions/griddedScoresEnvarsTransformed.RDS")
datNW <- meanR(datNW, species.groups)

datNW$area <- extract(a, datNW[ , c("x", "y")])

datNW$existing <- extract(existing_rast, datNW[ , c("x", "y")])
datNW$proposed <- extract(proposed_rast, datNW[ , c("x", "y")])
datNW$MPA <- datNW$existing + datNW$proposed

p <- plotGraph(this.data = datNW)

pdf("mpaPropUnweighted.pdf", width = 4.5, height = 2)
print(p)
dev.off()

# ---------------------------
# c. Jurisdiction
datW <- readRDS("../../meanPredictions/griddedScoresEnvarsTransformedColonyWeightedGroup.RDS")
datW <- meanR(datW, species.groups)

datW$area <- extract(a, datW[ , c("x", "y")])

datW$ccamlr <- extract(ccamlr_rast, datW[ , c("x", "y")])
datW$eez <- extract(eez_rast, datW[ , c("x", "y")])

p <- plotJurisdiction(this.data = datW)

pdf("jurisdictionWeighted.pdf", width = 4.5*0.96, height = 2*0.96)
print(p)
dev.off()

# ---------------------------
## Generate maps with AES

prj <- "+proj=laea +lat_0=-90  +lon_0=0 +datum=WGS84 +ellps=WGS84 +no_defs"

# Project MPAs
existing_mpa_proj <- spTransform(existing_mpa, CRSobj = prj)
proposed_mpa_proj <- spTransform(proposed_mpa, CRSobj = prj)

eez_proj <- spTransform(eez, CRSobj = prj)
ccamlr_proj <- spTransform(ccamlr, CRSobj = prj)

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

eez_proj <-gBuffer(eez_proj, width = 0, quadsegs = 1000) %>%
  crop(., polar.crop)

ccamlr_proj <-gBuffer(ccamlr_proj, width = 0, quadsegs = 1000) %>%
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
tiff("mpaMapWeighted.tiff", height = 7.5, width = 7.5, units = "in", res = 300) 
plot(d, axes = FALSE, bty = "n", box = FALSE, col = viridis(125), legend = F)
plot(ccamlr_proj, add = T, col = trans.none, border = "#33BBEE", lwd = 1.5)
plot(eez_proj, add = T, col = trans.none, border = "black", lwd = 1.5)
plot(existing_mpa_proj, add = T, col = trans.orange, border = "#EE7733", lwd = 0.5)
plot(proposed_mpa_proj, add = T, col = trans.red, border = "#EE3377", lwd = 0.5)
#plot(poly75, add = T, col = trans.none, border = "grey75")
plot(poly90, add = T, col = trans.none, border = "grey100", lwd = 1.5)
plot(wrld,
     col = "darkgrey",
     border = FALSE,
     add = TRUE)
lines(polar.crop, lwd = 2, add = T)
plot(d,
     col = viridis(125),
     legend.only=TRUE,
     horizontal = TRUE,
     legend.shrink=0.7,
     legend.args = list(text = 'Overall habitat importance')
)
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
tiff("mpaMapUnweighted.tiff", height = 7.5, width = 7.5, units = "in", res = 300) 
plot(d, axes = FALSE, bty = "n", box = FALSE, col = viridis(125), legend = F)
plot(ccamlr_proj, add = T, col = trans.none, border = "#33BBEE", lwd = 1.5)
plot(eez_proj, add = T, col = trans.none, border = "black", lwd = 1.5)
plot(existing_mpa_proj, add = T, col = trans.orange, border = "#EE7733", lwd = 0.5)
plot(proposed_mpa_proj, add = T, col = trans.red, border = "#EE3377", lwd = 0.5)
#plot(poly75, add = T, col = trans.none, border = "grey75")
plot(poly90, add = T, col = trans.none, border = "grey100")
plot(wrld,
     col = "darkgrey",
     border = FALSE,
     add = TRUE)
lines(polar.crop, lwd = 2, add = T)
plot(d,
     col = viridis(125),
     legend.only=TRUE,
     horizontal = TRUE,
     legend.shrink=0.7,
     legend.args = list(text = 'Overall habitat importance')
)
dev.off()