# Plots showing covariate mean, variance and sd against number of simulations
# for all species together

library(ggplot2)

## Set working directory
setwd("/perm_storage/home/ryan/RAATD_01/RAATD/")

## Get a list of all the files containing data for plotting
filesToGet <- list.files("~/../../../mnt/extra_storage/simDiscovery/plotData/", full.names = T)

## Read in and bind
dat <- do.call(rbind, lapply(filesToGet, readRDS))
dat$spec.stage <- paste(dat$species, dat$stage)

ggplot(dat, aes(x = nsim, y = mean, group = spec.stage, colour = species)) +
  geom_line() +
  facet_wrap(~Covariate, ncol = 2, scales = "free_y") +
  labs(x = "Number of simulations", y = "Covariate mean", title = "Mean - All species & stages")

ggplot(dat, aes(x = nsim, y = variance, group = spec.stage, colour = species)) +
  geom_line() +
  facet_wrap(~Covariate, ncol = 2, scales = "free_y") +
  labs(x = "Number of simulations", y = "Covariate variance", title = "Variance - All species & stages")

ggplot(dat, aes(x = nsim, y = sd, group = spec.stage, colour = species)) +
  geom_line() +
  facet_wrap(~Covariate, ncol = 2, scales = "free_y") +
  labs(x = "Number of simulations", y = "Covariate standard deviation", title = "Standard deviation - All species & stages")
