library(ggplot2)
library(raster)
library(dplyr)


this.species <- "ADPE"
this.stage <- "incubation"

## Set working directory
setwd(paste0("/perm_storage/home/ryan/RAATD_01/RAATD/", this.species, "/"))
source(paste0("/perm_storage/home/shared/github/raatd_modelling/simConfig/", this.species, "simConfig.R"))

## Load tracks and environmental data
real_tracks <- readRDS(paste0("~/../../../mnt/extra_storage/simDiscovery/sims/",this.species,"/simBatchTracks_",this.species,"_real.RDS"))
env <- envdata(this.stage) # load env data

## real tracks and env data
real_tracks <- real_tracks %>% dplyr::filter(stage==this.stage)
real_env <- extract(x=env, y=real_tracks[,c("lon","lat")])
all_env <- colnames(real_env)
##real_env <- bind_cols(real_tracks, as_tibble(real_env))

## simulated tracks and env data
sim_tracks <- lapply(seq_len(50), function(si) readRDS(paste0("~/../../../mnt/extra_storage/simDiscovery/sims/",this.species,"/simBatchTracks_",this.species,"_sim_",si,".RDS")) %>% dplyr::filter(stage==this.stage))
sim_env <- lapply(sim_tracks, function(z) extract(x=env, y=z[,c("lon","lat")]))##bind_cols(z, as_tibble(extract(x=env, y=z[,c("lon","lat")]))))

## ---------------------------------------------
## Simple plots of covariate variance (or sd) against nsims

hold <- data.frame()

for (i in 1:length(sim_env)) {
  print(paste(i, "of", length(sim_env)))
  this.sim <- sim_env[1:i]
  this.sim <- do.call(rbind, this.sim)
  vari <- data.frame(mean = apply(this.sim, 2, mean, na.rm = T),
                     variance = apply(this.sim, 2, var, na.rm = T),
                     sd = apply(this.sim, 2, sd, na.rm = T),
                     nsim = i*20)
  vari$Covariate <- row.names(vari)
  hold <- rbind(hold, vari)
  rm(this.sim)
}

hold$species <- this.species
hold$stage <- this.stage

## Save the plots
pdf(file = paste0("~/../../../mnt/extra_storage/simDiscovery/plots/envarSim_", this.species, "_", this.stage, ".pdf"), paper = "a4")
ggplot(hold, aes(x = nsim, y = mean, group = Covariate, colour = Covariate)) +
  geom_line() +
  facet_wrap(~Covariate, ncol = 2, scales = "free_y") +
  labs(x = "Number of simulations", y = "Covariate mean", title = paste(this.species, this.stage))

ggplot(hold, aes(x = nsim, y = variance, group = Covariate, colour = Covariate)) +
  geom_line() + facet_wrap(~Covariate, ncol = 2, scales = "free_y") +
  labs(x = "Number of simulations", y = "Covariate variance", title = paste(this.species, this.stage))

ggplot(hold, aes(x = nsim, y = sd, group = Covariate, colour = Covariate)) +
  geom_line() + facet_wrap(~Covariate, ncol = 2, scales = "free_y") +
  labs(x = "Number of simulations", y = "Covariate standard deviation", title = paste(this.species, this.stage))

dev.off()

## Save the output
saveRDS(hold, paste0("~/../../../mnt/extra_storage/simDiscovery/plotData/envarSim_", this.species, "_", this.stage, ".RDS"))