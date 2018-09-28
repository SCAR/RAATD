##' wrapLon
##'
##' Internal function not normally called by user
##'
##' @title wrapLon
##' @export

wrapLon <- function(lon, lmin = -180)
  (lon - lmin) %% 360 + lmin
