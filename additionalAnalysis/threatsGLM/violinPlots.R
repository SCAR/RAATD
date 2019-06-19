## GLMs of threats in core and non-core areas

library(raster)
library(vioplot)

setwd("~/RAATD_01/RAATD/additionalAnalysis/threatsGLM/")

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

# Calculate mean habitat importance
# Get data
dat <- readRDS("../../meanPredictions/griddedScoresEnvarsTransformedColonyWeighted.RDS")
dat <- meanR(dat, species.groups)

# Mask
landmask <- raster("../../antarcticMask/mask.grd")

# Ice mask
msk <- raster("/perm_storage/home/shared/data_extra/env_change/ice_mask.grd")
msk[msk <= 5] <- NA
msk <- projectRaster(msk, landmask)
msk[msk <= 5] <- NA

#-------------------------------
# Get threat data
#-------------------------------

fldr <- paste0("/perm_storage/home/shared/data_extra/env_change/")

#-------------------------------
# Fishing
load(paste0(fldr, "fishing_effort_.5deg.Rdata"))
r[is.na(r)] <- 0
fish <- r
# fish <- mask(fish, landmask)
rm(r)

dat$fish <- raster::extract(fish, dat[ , c("x", "y")])
dat$fishlog <- log(dat$fish + 1) # log transform
dat$fishcube <- dat$fish^(1/3) # Cube root transform
dat$fishsinh <- asinh(dat$fish) # inverse hyperbolic sine transformation

#-------------------------------
# SST change
sst <- raster(paste0(fldr, "sst_change.grd"))
dat$sst <- raster::extract(sst, dat[ , c("x", "y")])

#-------------------------------
# Ice change
ice <- raster(paste0(fldr, "ice_duration_change.grd"))
ice <- projectRaster(ice, to = landmask)
ice <- mask(ice, msk)
dat$ice <- raster::extract(ice, dat[ , c("x", "y")])

#-------------------------------
# Wind change
wind <- raster(paste0(fldr, "wind.grd"))
dat$wind <- raster::extract(wind, dat[ , c("x", "y")])

#-------------------------------
# Halpern cumulative impact
impact <- raster(paste0(fldr, "HalpernCumulativeThreats.grd"))
dat$impact <- raster::extract(impact, dat[ , c("x", "y")])

#-------------------------------
# Calculate quantiles for each region and
# flag the upper and lower quantiles as core and non-core

q75 <- stats::quantile(dat$MEAN, probs = c(0.75), na.rm = T)
q90 <- stats::quantile(dat$MEAN, probs = c(0.90), na.rm = T)

dat$core <- 0
dat[dat$MEAN >= q90 & !is.na(dat$MEAN), "core"] <- 1



#-------------------------------
# Violin plots
#-------------------------------

zerofish <- TRUE # Include nonzero fishing effort or not?

pdf("violins.pdf", width = 3.5*(10/6), height = 3.5*(10/6))

par(mfrow=c(2,2), mar=c(4,5,2,2), font.main = 1, cex.lab = 1, cex.main = 1)

#-----------
# Fishing

# the inverse hyperbolic sine transformation

if(zerofish) {
vioplot(x = dat[dat$core == 0, "fishsinh"],
        yaxt = "n",
        col = "darkgrey",
        plotCentre = "line",
        side = "left",
        add = F,
        # General:
        border = NA,
        names = "Outside AES | Inside AES   "
)

vioplot(x = dat[dat$core == 1, "fishsinh"],
        yaxt = "n",
        col = "#0077BB",
        plotCentre = "line",
        side = "right",
        add = T,
        # General:
        border = NA,
        names = "Outside AES | Inside AES   "
)

title(ylab = "Fishing effort (hours)")

axis(at = c(0, asinh(10), asinh(100), asinh(1000), asinh(10000), asinh(100000)),
  labels = c(0, 10, 100, 1000, 10000, "100000"),
  side = 2)

} else {
  vioplot(x = dat[dat$fishsinh > 0 & dat$core == 0, "fishsinh"],
          yaxt = "n",
          col = "darkgrey",
          plotCentre = "line",
          side = "left",
          add = F,
          # General:
          border = NA,
          names = "Outside AES | Inside AES   "
  )
  
  vioplot(x = dat[dat$fishsinh > 0 & dat$core == 1, "fishsinh"],
          yaxt = "n",
          col = "#0077BB",
          plotCentre = "line",
          side = "right",
          add = T,
          # General:
          border = NA,
          names = "Outside AES | Inside AES   "
  )
  
  title(ylab = "Nonzero fishing effort (hours)")
  
  axis(at = c(0, asinh(10), asinh(100), asinh(1000), asinh(10000), asinh(100000)),
       labels = c(0, 10, 100, 1000, 10000, "100000"),
       side = 2)
}

#-----------
# Ice

vioplot(x = dat[dat$core == 0, "ice"],
        # ylim = c(-50, 50),
        col = "darkgrey",
        plotCentre = "line",
        side = "left",
        add = F,
        # General:
        border = NA,
        names = "Outside AES | Inside AES   "
)

vioplot(x = dat[dat$core == 1, "ice"],
        # ylim = c(-50, 50),
        col = "#0077BB",
        plotCentre = "line",
        side = "right",
        add = T,
        # General:
        border = NA,
        names = "Outside AES | Inside AES   "
)

title(ylab = "Change in mean ice duration (days)")

#-----------
# SST

vioplot(x = dat[dat$core == 0, "sst"],
        # ylim = c(-1.5, 1.5),
        col = "darkgrey",
        plotCentre = "line",
        side = "left",
        add = F,
        # General:
        border = NA,
        names = "Outside AES | Inside AES   "
)

vioplot(x = dat[dat$core == 1, "sst"],
        # ylim = c(-1.5, 1.5),
        col = "#0077BB",
        plotCentre = "line",
        side = "right",
        add = T,
        # General:
        border = NA,
        names = "Outside AES | Inside AES   "
)

title(ylab = "Change in mean SST (C)")

#-----------
# Wind

vioplot(x = dat[dat$core == 0, "wind"],
        # ylim = c(-2, 2),
        col = "darkgrey",
        plotCentre = "line",
        side = "left",
        add = F,
        # General:
        border = NA,
        names = "Outside AES | Inside AES   "
)

vioplot(x = dat[dat$core == 1, "wind"],
        # ylim = c(-2, 2),
        col = "#0077BB",
        plotCentre = "line",
        side = "right",
        add = T,
        # General:
        border = NA,
        names = "Outside AES | Inside AES   "
)

title(ylab = "Change in mean wind speed (m/s)")

dev.off()
