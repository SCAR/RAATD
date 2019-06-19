## Mean habitat importance and AES based on
## pre-defined species groups

library(raster)
library(ggplot2)
library(viridis)

# library(Rcmdr)

setwd("~/RAATD_01/RAATD/additionalAnalysis/speciesGroups/")

#----------------------------
# Function to calculate AES
source("~/RAATD_01/RAATD/Code/function_aesPoly.R")

#----------------------------
# 1. Unweighted
dat <- readRDS("../../meanPredictions/griddedScoresEnvarsTransformed.RDS")

# Antarctic
dat$ant <- rowMeans(dat[, c("ADPE", "ANPE", "CRAS", "EMPE", "HUWH", "SOES", "WESE")])
ggplot(dat, aes(x, y, fill = ant)) + geom_raster() + coord_quickmap() +
  scale_fill_distiller(palette = "Spectral") +
  ggtitle("Antarctic - Not Weighted")

# Subantarctic
dat$sub <- rowMeans(dat[, c("ANFS", "BBAL", "DMSA", "GHAL", "HUWH", "KIPE", "LMSA", "MAPE.ROPE", "SOES", "WAAL", "WHCP")])
ggplot(dat, aes(x, y, fill = sub)) + geom_raster() + coord_quickmap() +
  scale_fill_distiller(palette = "Spectral") +
  ggtitle("Subantarctic - Not Weighted")

# Create rasters
antR <- rasterFromXYZ(dat[ , c("x", "y", "ant")])
subR <- rasterFromXYZ(dat[ , c("x", "y", "sub")])

# Create AES polygons
antRpoly <- aesPoly(antR, 0.75)
subRpoly <- aesPoly(subR, 0.75)

# Plot
trans.blue <- rgb(0, 119, 187, alpha = 255/3, maxColorValue = 255)

plot(antR)
plot(antRpoly, col = trans.blue, border = "black", add = T)

plot(subR)
plot(subRpoly, col = trans.blue, border = "black", add = T)

# What about max of the two values?
dat$max <- do.call(pmax, dat[c("ant", "sub")])
maxR <- rasterFromXYZ(dat[ , c("x", "y", "max")])

maxPoly75 <- aesPoly(rast = maxR, percentile.value = 0.75)
maxPoly90 <- aesPoly(rast = maxR, percentile.value = 0.90)

trans.none <- rgb(0, 119, 187, alpha = 255/3, maxColorValue = 255)

png("MaxMeansUnweighted.png", width = 2000, height = 800, units = "px", res = 150)
plot(maxR, col = viridis(125))
plot(maxPoly75, add = T, col = trans.none, border = "black")
plot(maxPoly90, add = T, col = "gold", border = FALSE)
dev.off()

#----------------------------
# 2. Weighted
dat <- readRDS("../../meanPredictions/griddedScoresEnvarsTransformedColonyWeighted.RDS")

# Antarctic
dat$ant <- rowMeans(dat[, c("ADPE", "ANPE", "CRAS", "EMPE", "HUWH", "SOES", "WESE")])
ggplot(dat, aes(x, y, fill = ant)) + geom_raster() + coord_quickmap() +
  scale_fill_distiller(palette = "Spectral") +
  ggtitle("Antarctic - Not Weighted")

# Subantarctic
dat$sub <- rowMeans(dat[, c("ANFS", "BBAL", "DMSA", "GHAL", "HUWH", "KIPE", "LMSA", "MAPE.ROPE", "SOES", "WAAL", "WHCP")])
ggplot(dat, aes(x, y, fill = sub)) + geom_raster() + coord_quickmap() +
  scale_fill_distiller(palette = "Spectral") +
  ggtitle("Subantarctic - Not Weighted")

# Create rasters
antR <- rasterFromXYZ(dat[ , c("x", "y", "ant")])
subR <- rasterFromXYZ(dat[ , c("x", "y", "sub")])

# Create AES polygons
antRpoly <- aesPoly(antR, 0.75)
subRpoly <- aesPoly(subR, 0.75)

# Plot
trans.blue <- rgb(0, 119, 187, alpha = 255/3, maxColorValue = 255)

plot(antR)
plot(antRpoly, col = trans.blue, border = "black", add = T)

plot(subR)
plot(subRpoly, col = trans.blue, border = "black", add = T)

# What about max of the two values?
dat$max <- do.call(pmax, dat[c("ant", "sub")])
maxR <- rasterFromXYZ(dat[ , c("x", "y", "max")])

maxPoly75 <- aesPoly(rast = maxR, percentile.value = 0.75)
maxPoly90 <- aesPoly(rast = maxR, percentile.value = 0.90)

trans.none <- rgb(0, 119, 187, alpha = 255/3, maxColorValue = 255)

png("MaxMeansWeighted.png", width = 2000, height = 800, units = "px", res = 150)
plot(maxR, col = viridis(125))
plot(maxPoly75, add = T, col = trans.none, border = "black")
plot(maxPoly90, add = T, col = "gold", border = FALSE)
dev.off()

#----------------------------
# What about cluster analysis
datT <- dat[complete.cases(dat[ , 3:18]), ]
datT <- as.data.frame(t(datT[, 3:18]))

# Results of k = 2 are rubbish; k = 3 looks okay
clusters <- Kmeans(datT, centers = 3, method = "manhattan", iter.max = 100)

# Dendrogram of hierarchical clustering, average linkage gives the highest CCC
hcl <- hclust(dist(datT, method = "manhattan"), method = "average")
cor(cophenetic(hcl), dist(datT, method = "manhattan"))
pdf("cluster.pdf", height = 6.1, width = 6.1)
plot(hcl, hang = -1)
dev.off()