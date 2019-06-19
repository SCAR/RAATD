## Compare untransformed with percentile transformed predictions

## Ryan Reisinger

library(ggplot2)

setwd("~/RAATD_01/RAATD/additionalAnalysis/compareTransformation/")

# Get the data
datU <- readRDS("../../meanPredictions/griddedScoresEnvarsUntransformed.RDS")
datT <- readRDS("../../meanPredictions/griddedScoresEnvars.RDS")

datT$t <- "Transformed"
datU$t <- "Untransformed"

datT <- datT[complete.cases(datT), ]

ggplot(dat = datT, aes(y = BBAL, x = rank(datT$BBAL))) +
  geom_point(color = "blue") +
  geom_point(dat = datU, aes(y = BBAL*100, x = rank(datU$BBAL)), color = "red", inherit.aes = F)
