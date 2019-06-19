## Simulate tracks and write output to file, in batches
## B. Raymond, R. Reisinger, S. Wotherspoon and I. Jonsen
## Updated 2018-03-14
##----------------------------------------

## Which species
this.species <- "WESE"

##----------------------------------------

## Set working directory
setwd(paste0("~/RAATD_01/RAATD/", this.species, "/"))

## Must load dplyr after raster
options(stringsAsFactors=FALSE)
library(raster)
library(dplyr)
library(multidplyr)
library(TMB)
library(Marseille)

## Species specific configuration data
source(paste0("/perm_storage/home/shared/github/raatd_data/data_split/config/", this.species, "config.r"))
source(paste0("/perm_storage/home/shared/github/raatd_modelling/simConfig/", this.species, "simConfig.R"))

## Read species metadata
meta.csv <- "/perm_storage/home/shared/github/raatd_data/metadata/SCAR_Metadata_2017_forWEBDAV.csv"
df.meta <- read.csv(meta.csv,header=TRUE,stringsAsFactors=FALSE)

## Load the stage-wise fitted tracks discarding errors
fits <- local({
  ssm_by_stage_final <- readRDS(stage.rdata)
  ssm_by_stage_final %>%
    prepare() %>%
    filter(length(ssm[[1]]) != 1 & !inherits(ssm[[1]],"character") & !inherits(ssm[[1]],"try-error")) %>%
    ungroup()
})

## Merge with metadata
fits <- left_join(
  fits,
  df.meta %>%
    select(id=individual_id,
           type=device_type,
           deploy_site=deployment_site,
           deploy_lon=deployment_decimal_longitude,
           deploy_lat=deployment_decimal_latitude),
  by=c("id"))

## Merge with buffer information
fits <- left_join(fits,buff,by=c("type"))

## Some species still  have lingering duplicate stage names
fits$stage[fits$stage == "incubation.3"] <- "incubation" # WHCP
fits$stage[fits$stage == "post-breeding.3.4"] <- "post-breeding" # WAAL

## Cleanup
rm(df.meta)

## The prediction region
bound

## Create raster that defines region of interest
grd <- raster(extent(bound),resolution=c(0.1,0.1))
crs(grd) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

## Mask - prevent tracks on land and restrict latitudes to [-90,0]
isSea <- gshhsMask(res="0.1",land=FALSE,latmin=-90,latmax=0)

## Define fixed points
df.fixed <- fits %>%
  dplyr::rename(site=deploy_site) %>%
  group_by(site) %>%
  summarise(lon=mean(wrapLon(deploy_lon)),lat=mean(deploy_lat))

## The first point on the track, any points along the track that are
## within the prescribed range of any point in the table, and any points
## on land within 200km of of a point in the table are marked as fixed
library(geosphere)
fixedPoints <- function(trk,buff) {
  land <- !isSea(trk$date,trk$lon,trk$lat)
  ## Wrap the longitudes into -180,180
  trk <- cbind(wrapLon(trk$lon),trk$lat)
  ## Distance to nearest point in df.fixed
  R <- Inf
  for(k in seq_len(NROW(df.fixed)))
    R <- pmin(R,distCosine(trk,df.fixed[k,c("lon","lat")]))
  ## Force first point to be fixed
  fixed <- logical(NROW(trk))
  fixed[1] <- TRUE
  ## Test the track against each site.
  fixed | (R < buff) | (land & (R < 200000))
}


## Define the track simulation function
simulateTracks <- function(fit, nsim=50) {
  ## Determine fixed points
  fixed <- fixedPoints(fit$ssm[[1]]$predicted, fit$tbuff)
  ## Attempt to simulate a single track
  simulate <- function(retry=20) {
    for(k in seq_len(retry)) {
      trk <- dcrwSimulate(fit$ssm[[1]], fixed=fixed, point.check=isSea)
      if(!is.null(trk)) return(cbind(fixed=fixed, trk))
    }
    ## Return zero length track on failure
    data.frame()
  }
  ## Simulate nsims tracks
  trks <- data.frame(sim=seq_len(nsim)) %>% group_by(sim) %>% do(simulate())
  ## Adjoin descriptors if successful
  if(nrow(trks)>0)
    data.frame(id=fit$id,
               track=fit$track,
               stage=fit$stage,
               deploy_site=fit$deploy_site,
               deploy_lon=fit$deploy_lon,
               deploy_lat=fit$deploy_lat,
               trks)
  else
    trks
}

## Attempt to simulate `nsim` new tracks from each fitted track
cl <- parallel::detectCores() - 1

grp <- rep(1:cl, length.out = nrow(fits))
fits.tmp <- bind_cols(tibble(grp = grp), fits)

cl <- create_cluster(cl)
cl %>% 
  cluster_copy(simulateTracks) %>%
  cluster_copy(isSea) %>%
  cluster_copy(df.fixed) %>%
  cluster_copy(fixedPoints) %>%
  cluster_library("geosphere") %>%
  cluster_library("Marseille") %>%
  cluster_library("dplyr")



# ---------------------------
# Run simulations in batches, writing output for later use
# Run 50 x batches of 20 sims
for (k in 1:50) {
  print(paste("Simulating batch", k, "of 50"))
  siml <- fits.tmp %>%
    partition(., grp, track, cluster = cl) %>%
    do(simulateTracks(., nsim = 20)) %>%
    collect()
  
  # Remove tracks which don't meet the land mask
  siml <- siml[siml$fixed | isSea(siml$date,siml$lon,siml$lat),]
  
  # Add batch ID
  siml$batch <- k
  
  # Save output
  ifelse(!dir.exists(paste0("~/../../../mnt/extra_storage/simDiscovery/sims/", this.species, "/")),
         dir.create(paste0("~/../../../mnt/extra_storage/simDiscovery/sims/", this.species, "/")),
         FALSE) # check that dir exists
  saveRDS(siml, paste0("~/../../../mnt/extra_storage/simDiscovery/sims/", this.species, "/simBatchTracks_", this.species, "_sim_", k, ".RDS"))
  
  ## Extract the actual tracks
  actl <- fits %>%
    rowwise() %>%
    do(data.frame(id=.$id,
                  track=.$track,
                  stage=.$stage,
                  deploy_site=.$deploy_site,
                  deploy_lon=.$deploy_lon,
                  deploy_lat=.$deploy_lat,
                  sim=0, ## real tracks have sim=0
                  fixed=fixedPoints(.$ssm$predicted[,2:3],.$tbuff),
                  .$ssm$predicted[,1:3])) %>%
    ungroup()
  
  if (k == 1) {
    saveRDS(actl, paste0("~/../../../mnt/extra_storage/simDiscovery/sims/", this.species, "/simBatchTracks_", this.species, "_real.RDS"))
  }
  
  # -----------------------------------------
  # Write per-stage dataframes with environmental data, for each batch
  
  for (j in 1:nrow(stages)) {
    this.stage <- stages$stage[j]
    print(this.stage)
    this.siml <- siml[siml$stage == this.stage, ]
    this.actl <- actl[actl$stage == this.stage, ]
    
    if (nrow(this.siml) > 0) {
      env <- envdata(this.stage) # load env data
      env.avail <- extract(x = env,  y = this.siml[ , c("lon", "lat")]) # extract envdata
      saveRDS(env.avail, paste0("~/../../../mnt/extra_storage/simDiscovery/sims/", this.species, "/simBatchEnv_", this.species, "_sim_", k, "_", this.stage, ".RDS"))
      
      # ----
      # Extract envars for actual tracks
      # Do this only once
      if (k == 1) {
        env.used <- extract(x = env,  y = this.actl[ , c("lon", "lat")]) # extract envdata
        saveRDS(env.used, paste0("~/../../../mnt/extra_storage/simDiscovery/sims/", this.species, "/simBatchEnv_", this.species, "_real_", this.stage, ".RDS"))
      }
    }
  }
}