## Simulate tracks and create gridded output for habitat modelling
## B. Raymond, R. Reisinger, S. Wotherspoon and I. Jonsen
## Updated 2017-12-12
##----------------------------------------

## Which species
this.species <- "KIPE"
##----------------------------------------

## Set working directory
setwd(paste0("~/RAATD_01/RAATD/", this.species, "/"))

## Must load dplyr after raster
memory.limit(128000)
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

## Environmental Data

## It is assumed that environmental data has been extracted for each stage
## as a separate step, and that the species configuration file defines a
## function `envdata` that loads the environmental data corresponding to a particular stage
envdata

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

siml <- fits.tmp %>%
  partition(., grp, track, cluster = cl) %>%
  do(simulateTracks(., nsim = 50)) %>%
  collect()

#parallel::stopCluster(cl)

## Save the simulated tracks, if neccessary
# saveRDS(siml, paste0(this.species, "_simTracks.RDS"))

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

## Show map of approx region

## Map tracks by departure sites
library(ggmap)
ggmap(get_map(location=rbind(pmax(-180,pmin(179.99,bound[1:2])),bound[3:4]),
              maptype="toner-lite",source="stamen",zoom=3))+
  geom_path(mapping=aes(x=lon,y=lat,group=track),colour="dodgerblue",alpha=0.1,data=siml)+
  geom_path(mapping=aes(x=lon,y=lat,group=track),colour="darkorange",alpha=0.1,data=actl)+
  geom_point(mapping=aes(x=lon,y=lat),col="black",cex=0.5,inherit.aes=FALSE,
             data=actl %>% select(deploy_site,lon=deploy_lon,lat=deploy_lat) %>% distinct())+
  facet_wrap(~deploy_site)+theme(legend.position="none")

## Cleanup
rm(fits)

## Function to bin tracks
if(this.species != "BBAL"){
  
  binTracks <- function(trks) {
    
    trkcount <- function(x,...) length(unique(na.omit(x)))
    
    ## Remove non-fixed points that do not meet the mask
    trks <- trks[trks$fixed | isSea(trks$date,trks$lon,trks$lat),]
    ## Separate into actual and simulated tracks
    actl <- trks[trks$sim==0,]
    siml <- trks[trks$sim>0,]
    ## Compute mean deployment location
    deploy <- c(mean(wrapLon(trks$deploy_lon,xmin(grd))),mean(trks$deploy_lat))
    
    list(
      avail=rasterize(trks[,c("lon","lat")],grd,
                      fun="count",
                      background=0),
      usage=rasterize(actl[,c("lon","lat")],grd,
                      fun="count",
                      background=0),
      nactl=rasterize(actl[,c("lon","lat")],grd,
                      field=actl$track,
                      fun=trkcount,
                      background=0),
      nsiml=rasterize(siml[,c("lon","lat")],grd,
                      field=interaction(siml$track,siml$sim,drop=TRUE),
                      fun=trkcount,
                      background=0),
      dist=distanceFromPoints(grd,deploy))
  }
} else {
  ## BBAL needs a different function to prevent hanging
  source("../functions_binSplitBBAL_02.R")
}

## Create rasters for each site and stage
rstr <- bind_rows(siml,actl) %>%
  filter(!fixed) %>%
  mutate(lon=wrapLon(lon,xmin(grd))) %>%
  group_by(deploy_site,stage) %>%
  do(rasters=binTracks(.))

## Cleanup
rm(actl,siml)

## Data

## The binned tracks are converted back to dataframes for model fitting.
## This and the binning step could be merged.

## Availability

## The availability data contains availability, distance from departure
## site, site, stage, and environmental data

mkAvail <- function(bin) {
  r <- bin$rasters
  
  ## INSERT range limitation here
  ## Missing values in any raster used to construct df
  ## cause the cell to excluded
  R <- 1.2*quantile(r$dist[r$avail>0],0.95)
  R <- max(R,5000)
  r$dist[r$dist > R] <- NA
  
  ## Convert avail and dist rasters to a dataframe.
  df <- as.data.frame(setNames(stack(r$avail,r$dist),
                               c("avail","dist_deploy")),
                      xy=TRUE,na.rm=TRUE)
  if(nrow(df)) {
    colnames(df)[1:2] <- c("lon","lat")
    ## Add site and environmental predictors.
    env <- envdata(bin$stage)
    df <- cbind(df,
                site=bin$deploy_site,
                stage=bin$stage,
                extract(env,cbind(wrapLon(df$lon),df$lat)))
  }
  df
}

## Data can be generated for all stages
local({
  ## Make data for stage
  avail <- rstr %>%
    rowwise() %>%
    do(mkAvail(.)) %>%
    ungroup()
  ## Write as RData
  save(avail,file=paste0(paste(spcode,"avail",sep="-"),".Rdata"))
})

## Or stage by stage
local({
  for(s in unique(rstr$stage)) {
    ## Make data for stage
    avail <- rstr %>%
      filter(stage==s) %>%
      rowwise() %>%
      do(mkAvail(.)) %>%
      ungroup()
    ## Write as RData
    save(avail,file=paste0(paste(spcode,"avail",s,sep="-"),".Rdata"))
  }
})

## Usage

## The usage data contains availability, usage, distance from departure
## site, site, stage, and environmental data for cells that were available
mkUsage <- function(bin) {
  r <- bin$rasters
  
  ## INSERT range limitation here
  ## Missing values in any raster used to construct df
  ## cause the cell to excluded
  R <- 1.2*quantile(r$dist[r$avail>0],0.95)
  R <- max(R,5000)
  r$dist[r$dist > R] <- NA
  
  ## Limit to available cells
  r$avail[r$avail==0] <- NA
  
  ## Convert avail, usage and dist rasters to a dataframe.
  df <- as.data.frame(setNames(stack(r$avail,r$usage,r$dist),
                               c("avail","usage","dist_deploy")),
                      xy=TRUE,na.rm=TRUE)
  
  if(nrow(df)) {
    colnames(df)[1:2] <- c("lon","lat")
    ## Add site and environmental predictors.
    env <- envdata(bin$stage)
    df <- cbind(df,
                site=bin$deploy_site,
                stage=bin$stage,
                extract(env,cbind(wrapLon(df$lon),df$lat)))
  }
  df
}

## For all stages
local({
  ## Make data for stage
  usage <- rstr %>%
    rowwise() %>%
    do(mkUsage(.)) %>%
    ungroup()
  ## Write as RData
  save(usage,file=paste0(paste(spcode,"usage",sep="-"),".Rdata"))
})

## Or stage by stage
local({
  for(s in unique(rstr$stage)) {
    ## Make data for stage
    usage <- rstr %>%
      filter(stage==s) %>%
      rowwise() %>%
      do(mkUsage(.)) %>%
      ungroup()
    ## Write as RData
    save(usage,file=paste0(paste(spcode,"usage",s,sep="-"),".Rdata"))
  }
})

## Spent

## The spent data uses the number of track points in each cell as a proxy
## for time spent, and contains time spent, distance from departure
## site, site, stage, and environmental data for cells visited by a track
mkSpent <- function(bin) {
  r <- bin$rasters
  
  ## Convert to time spent per track retaining only visited cells
  r$nactl[r$nactl==0] <- NA
  spent <- r$usage/r$nactl
  
  ## Convert spent and dist rasters to a dataframe.
  df <- as.data.frame(setNames(stack(spent,r$dist),
                               c("spent","dist_deploy")),
                      xy=TRUE,na.rm=TRUE)
  
  if(nrow(df)) {
    colnames(df)[1:2] <- c("lon","lat")
    
    ## Add site and environmental predictors.
    env <- envdata(bin$stage)
    df <- cbind(df,
                site=bin$deploy_site,
                stage=bin$stage,
                extract(env,cbind(wrapLon(df$lon),df$lat)))
  }
  df
}

## For all stages
local({
  ## Make data for stage
  spent <- rstr %>%
    rowwise() %>%
    do(mkSpent(.)) %>%
    ungroup()
  ## Write as RData
  save(spent,file=paste0(paste(spcode,"spent",sep="-"),".Rdata"))
})

## Or stage by stage
local({
  for(s in unique(rstr$stage)) {
    ## Make data for stage
    spent <- rstr %>%
      filter(stage==s) %>%
      rowwise() %>%
      do(mkSpent(.)) %>%
      ungroup()
    ## Write as RData
    save(spent,file=paste0(paste(spcode,"spent",s,sep="-"),".Rdata"))
  }
})

#--------------------------------------------------------------------
