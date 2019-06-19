## Threats in core and non-core areas

library(raster)
library(lme4)
library(beanplot)

library(tidyverse)

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
# Test for nonzero change
#-------------------------------

foo <- function(which.var, k = 1000) {
  
  # which.var is a vector of the change values
  # k is the number of repeats
  
  out <- vector(length = k)
  
  for(i in 1:k) {
    this.dat <- sample(which.var, 10000)
    pos <- length(this.dat[this.dat > 0])/10000
    out[i] <- pos
  }
  
  return(out)
}

# Apply the function

hold <- foo(dat$ice)
# Plot the result
hist(hold)

# Test value (extreme values, greater than 0.95 or less than 0.05 indicate the hypothesis of no change is rejected)
length(hold[hold<0.5])/k

#BR had in mind:
hold <- foo(dat$ice, k = 1)
min(hold, 1-hold)

#-------------------------------
# Simple histogram plots
#------------------------------

library(ggplot2)
library(gridExtra)

# Custom theme
theme_rr <- function () { 
  theme_linedraw() %+replace% 
    theme(
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        strip.background = element_blank(),
        strip.text = element_blank()
    )
}

# Calculate group medians
compare_dat <- dat %>%
  group_by(core) %>%
  summarise(fish.m = asinh(median(fish, na.rm = T)),
            ice.m = median(ice, na.rm = T),
            sst.m = median(sst, na.rm = T),
            wind.m = median(wind, na.rm = T))


p1 <- ggplot(dat[!is.na(dat$core), ], aes(x = fishsinh, colour = as.factor(core), fill = as.factor(core))) +
  geom_histogram(bins = 100) +
  facet_wrap(~core, ncol = 1, scales = "free_y") +
  labs(x = "Fishing effort (hours)",
       y = "Number of cells") +
  scale_x_continuous(breaks = c(0, asinh(10), asinh(100), asinh(1000), asinh(10000), asinh(100000)),
                     labels = c(0, 10, 100, 1000, 10000, "100000")) +
  scale_y_continuous(breaks = c(0, 10, 100, 1000, 10000),
    trans = "log1p") +
  geom_vline(data = compare_dat, aes(xintercept = fish.m), size = 1, colour = "#CC3311") +
  scale_colour_manual(values = c("darkgrey", "#0077BB")) +
  scale_fill_manual(values = c("darkgrey", "#0077BB")) +
  guides(colour = "none", fill = "none") +
  theme_rr()

p2 <- ggplot(dat[!is.na(dat$core), ], aes(x = ice, colour = as.factor(core), fill = as.factor(core))) +
  geom_histogram(bins = 100) +
  facet_wrap(~core, ncol = 1, scales = "free_y") +
  geom_vline(xintercept = 0, size = 0.5) +
  geom_vline(data = compare_dat, aes(xintercept = ice.m), size = 1, colour = "#CC3311") +
  labs(x = "Change in mean ice duration (days)",
       y = "Number of cells") +
  scale_colour_manual(values = c("darkgrey", "#0077BB")) +
  scale_fill_manual(values = c("darkgrey", "#0077BB")) +
  guides(colour = "none", fill = "none") +
  theme_rr()

p3 <- ggplot(dat[!is.na(dat$core), ], aes(x = sst, colour = as.factor(core), fill = as.factor(core))) +
  geom_histogram(bins = 100) +
  facet_wrap(~core, ncol = 1, scales = "free_y") +
  geom_vline(xintercept = 0, size = 0.5) +
  geom_vline(data = compare_dat, aes(xintercept = sst.m), size = 1, colour = "#CC3311") +
  labs(x = "Change in mean SST (C)",
       y = "Number of cells") +
  scale_colour_manual(values = c("darkgrey", "#0077BB")) +
  scale_fill_manual(values = c("darkgrey", "#0077BB")) +
  guides(colour = "none", fill = "none") +
  theme_rr()


p4 <- ggplot(dat[!is.na(dat$core), ], aes(x = wind, colour = as.factor(core), fill = as.factor(core))) +
  geom_histogram(bins = 100) +
  facet_wrap(~core, ncol = 1, scales = "free_y") +
  geom_vline(xintercept = 0, size = 0.5) +
  geom_vline(data = compare_dat, aes(xintercept = wind.m), size = 1, colour = "#CC3311") +
  labs(x = "Change in mean wind speed (m/s)",
       y = "Number of cells") +
  scale_colour_manual(values = c("darkgrey", "#0077BB")) +
  scale_fill_manual(values = c("darkgrey", "#0077BB")) +
  guides(colour = "none", fill = "none") +
  theme_rr()



pdf("histograms.pdf", width = 3.5/0.77, height = (2/0.77)*2)
grid.arrange(p1, p2, p3, p4, ncol=2)
dev.off()




#-------------------------------
# Permutation tests
#-------------------------------

library(coin)

dat$core <- as.factor(dat$core)

# Ice
ice.test <- independence_test(formula = ice ~ core, data = dat[complete.cases(dat[ , c("ice", "core")]), ],
                  distribution = approximate(nresample = 10000))

wilcox_test(ice ~ core, data = dat[complete.cases(dat[ , c("ice", "core")]), ],
             distribution = approximate(nresample = 10000))

# SST
sst.test <- independence_test(formula = sst ~ core, data = dat[complete.cases(dat[ , c("sst", "core")]), ],
                  distribution = approximate(nresample = 10000))

wilcox_test(sst ~ core, data = dat[complete.cases(dat[ , c("sst", "core")]), ],
            distribution = approximate(nresample = 10000))

# Fishing
fish.test <- independence_test(formula = fish ~ core, data = dat[complete.cases(dat[ , c("fish", "core")]), ],
                  distribution = approximate(nresample = 10000))

wilcox_test(fish ~ core, data = dat[complete.cases(dat[ , c("fish", "core")]), ],
            distribution = approximate(nresample = 10000))

# Wind
wind.test <- independence_test(formula = wind ~ core, data = dat[complete.cases(dat[ , c("wind", "core")]), ],
                  distribution = approximate(nresample = 10000))

wilcox_test(wind ~ core, data = dat[complete.cases(dat[ , c("wind", "core")]), ],
            distribution = approximate(nresample = 10000))


#-------------------------------
# Bean plots
#-------------------------------

pdf("beans6.pdf", width = 3.5*(10/6), height = 3.5*(10/6))

par(mfrow=c(2,2), mar=c(4,5,2,2), font.main = 1, cex.lab = 1, cex.main = 1)

#-----------
# Fishing

fish.label <-
  paste0("Outside AES | Inside AES", "\n",
         "Z = ", round(statistic(fish.test), 2), ", p = ", round(pvalue(fish.test)[1], 3))

# Log transform
#dat$fishlog[!is.finite(dat$fishlog)] <- NA

# beanplot(fishlog ~ core, data = dat[dat$fishlog > 0, ],
#          main = "", ylab = "Fishing effort (log (hours))", xlab="Outside AES | Inside AES", side = "both",
#          names = "",
#          bw="nrd0",
#          log = "",
#          #ylim = c(0, 10),
#          overallline = "mean",
#          what=c(0,1,0,0), frame.plot = F,
#          maxwidth = 0.8,
#          border = NA, col = list("darkgrey", "#0077BB"))

# Cube root transform
# beanplot(fishcube ~ core, data = dat[complete.cases(dat[ , c("fishcube", "core")]), ],
#          main = "", ylab = "Fishing effort (cube root (hours))", xlab="Outside AES | Inside AES", side = "both",
#          names = "",
#          bw="nrd0",
#          log = "",
#          ylim = c(0, 10),
#          overallline = "mean",
#          what=c(0,1,0,0), frame.plot = F,
#          maxwidth = 0.8,
#          border = NA, col = list("darkgrey", "#0077BB"))

# the inverse hyperbolic sine transformation
beanplot(fishsinh~ core, data = dat[dat$fishsinh > 0, ],
         main = "",
         ylab = "Nonzero fishing effort (hours)",
         # xlab = fish.label,
         xlab = "Outside AES | Inside AES   ",
         side = "both",
         names = "",
         bw="nrd0",
         log = "",
         yaxt = "n",
         #ylim = c(0, 10),
         # General params from here
         overallline = "median",
         what=c(0,1,1,0), # the total average line, the beans, the bean average, and the beanlines
         frame.plot = F,
         maxwidth = 0.8,
         ll = 0.1, # beanline width
         maxstripline = 0.1, # max width of the beanlines
         beanlines = "quantiles", # what overall line for each bean
         beanlinewd = 1, # width of the quantile/average lines
         border = NA,
         # cols = area of the beans, the lines inside the bean, the lines outside the bean, and the average line per bean
         col = list(c("darkgrey", adjustcolor( "black", alpha.f = 0.008), adjustcolor( "black", alpha.f = 0.008), "grey30"), #cols for first group
                    c("#0077BB", adjustcolor( "black", alpha.f = 0.008), adjustcolor( "black", alpha.f = 0.008), "grey30"))) #cols for second group

axis(at = c(0, asinh(10), asinh(100), asinh(1000), asinh(10000), asinh(100000)),
  labels = c(0, 10, 100, 1000, 10000, "100000"),
  side = 2)

#-----------
# Ice
ice.label <-
  paste0("Outside AES | Inside AES", "\n",
         "Z = ", round(statistic(ice.test), 2), ", p = ", round(pvalue(ice.test)[1], 3))

beanplot(ice ~ core, data = dat[complete.cases(dat[ , c("ice", "core")]), ],
         main = "",
         ylab = "Change in mean ice duration (days)",
         # xlab = ice.label,
         xlab = "Outside AES | Inside AES   ",
         side = "both",
         names = "",
         # ylim = c(-50, 50),
         bw="nrd",
         # General params from here
         overallline = "median",
         what=c(0,1,1,0), # the total average line, the beans, the bean average, and the beanlines
         frame.plot = F,
         maxwidth = 0.8,
         ll = 0.1, # beanline width
         maxstripline = 0.1, # max width of the beanlines
         beanlines = "quantiles", # what overall line for each bean
         beanlinewd = 1, # width of the quantile/average lines
         border = NA,
         # cols = area of the beans, the lines inside the bean, the lines outside the bean, and the average line per bean
         col = list(c("darkgrey", adjustcolor( "black", alpha.f = 0.008), adjustcolor( "black", alpha.f = 0.008), "grey30"), #cols for first group
                                 c("#0077BB", adjustcolor( "black", alpha.f = 0.008), adjustcolor( "black", alpha.f = 0.008), "grey30"))) #cols for second group


#-----------
# SST

sst.label <-
  paste0("Outside AES | Inside AES", "\n",
         "Z = ", round(statistic(sst.test), 2), ", p = ", round(pvalue(sst.test)[1], 3))

beanplot(sst ~ core, data = dat[complete.cases(dat[ , c("sst", "core")]), ],
         main = "", ylab = "Change in mean SST (C)",
         # xlab = sst.label,
         xlab = "Outside AES | Inside AES   ",
         side = "both",
         names = "",
         bw="nrd",
         # ylim = c(-1.5, 1.5),
         # General params from here
         overallline = "median",
         what=c(0,1,1,0), # the total average line, the beans, the bean average, and the beanlines
         frame.plot = F,
         maxwidth = 0.8,
         ll = 0.1, # beanline width
         maxstripline = 0.1, # max width of the beanlines
         beanlines = "quantiles", # what overall line for each bean
         beanlinewd = 1, # width of the quantile/average lines
         border = NA,
         # cols = area of the beans, the lines inside the bean, the lines outside the bean, and the average line per bean
         col = list(c("darkgrey", adjustcolor( "black", alpha.f = 0.008), adjustcolor( "black", alpha.f = 0.008), "grey30"), #cols for first group
                    c("#0077BB", adjustcolor( "black", alpha.f = 0.008), adjustcolor( "black", alpha.f = 0.008), "grey30"))) #cols for second group


#-----------
# Wind

wind.label <-
  paste0("Outside AES | Inside AES", "\n",
         "Z = ", round(statistic(wind.test), 2), ", p = ", round(pvalue(wind.test)[1], 3))

beanplot(wind ~ core, data = dat[complete.cases(dat[ , c("wind", "core")]), ],
         main = "",
         ylab = "Change in mean wind speed (m/s)",
         # xlab = wind.label,
         xlab = "Outside AES | Inside AES   ",
         side = "both",
         names = "",
         bw="nrd",
         # ylim = c(-2, 2),
         # General params from here
         overallline = "median",
         what=c(0,1,1,0), # the total average line, the beans, the bean average, and the beanlines
         frame.plot = F,
         maxwidth = 0.8,
         ll = 0.1, # beanline width
         maxstripline = 0.1, # max width of the beanlines
         beanlines = "quantiles", # what overall line for each bean
         beanlinewd = 1, # width of the quantile/average lines
         border = NA,
         # cols = area of the beans, the lines inside the bean, the lines outside the bean, and the average line per bean
         col = list(c("darkgrey", adjustcolor( "black", alpha.f = 0.008), adjustcolor( "black", alpha.f = 0.008), "grey30"), #cols for first group
                    c("#0077BB", adjustcolor( "black", alpha.f = 0.008), adjustcolor( "black", alpha.f = 0.008), "grey30"))) #cols for second group

dev.off()

#-------------------------------
# Fishing proportions
#-------------------------------
library(viridis)

fishdat <- dat[!is.na(dat$fish), ]

# Fishing categories
fishdat$anyfish <- NA
fishdat[fishdat$fish == 0, "anyfish"] <- "nofishing"
fishdat[fishdat$fish > 0, "anyfish"] <- "fishing"

cairo_pdf("fishing_proportions.pdf", width = 2/0.68, height = 1.4/0.68)
ggplot(dat = fishdat, aes(x = as.factor(core), fill = anyfish)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = rev(viridis(2))) +
  theme_rr()
dev.off()

# Proportions
# In AES
nrow(fishdat[fishdat$core == 1 & fishdat$anyfish == "fishing", ]) / nrow(fishdat[fishdat$core == 1, ])

# Outside AES
nrow(fishdat[fishdat$core == 0 & fishdat$anyfish == "fishing", ]) / nrow(fishdat[fishdat$core == 0, ])


# Prepare .csv for publishing
fishpub <- select(fishdat,
                  x,
                  y,
                  core,
                  fish,
                  fishsinh)
names(fishpub) <- c("lon", "lat", "AES", "fishing.effort", "fishing.effort.asinh")

write.csv(fishpub, "Data_Fig2.csv", row.names = F)
