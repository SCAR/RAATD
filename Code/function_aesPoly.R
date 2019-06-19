## Function to calculate polygon of a given percentile value
## from a raster

# Ryan Reisinger
# September 2018

aesPoly <- function(rast, percentile.value = 0.75) {
  
  # rast = raster to calculate the polygon from
  # percentile.value = numeric giving the percentile value to use as threshold
  
  require(raster)
  require(sf)
  require(lwgeom) ## was liblwgeom, but you must mean lwgeom - BR Dec 2018
  require(smoothr)
  
  # Copy
  foo <- rast
  
  # Calculate quantiles
  q <- quantile(foo, p = percentile.value)
  
  # Threshold raster based on quantiles
  foo[foo < q] <- NA
  foo[!is.na(foo)] <- 1
  
  # Convert to polygons
  c <- rasterToPolygons(foo, dissolve = T)
  
  # Set CRS
  proj4string(c) <- "+proj=longlat +datum=WGS84 +no_defs"
  
  # Convert to SF object 
  r_poly <- st_as_sf(c)
  
  # Drop small, isolated cells
  area_thresh <- units::set_units(100*100, km^2)
  r_poly_dropped <- drop_crumbs(r_poly, area_thresh)
  
  # Fill small holes
  r_poly_filled <- fill_holes(r_poly_dropped, area_thresh)
  
  # And smooth the edges
  r_poly_smooth <- smoothr::smooth(r_poly_filled, method = "ksmooth", smoothness = 20)
  
  # Plot to check
  plot(c, col = NA, border = NA) # set up plot extent
  plot(r_poly_smooth, col = "#4DAF4A", border = "grey20", lwd = 1.5, add = TRUE)
  
  return(r_poly_smooth)
  
}
