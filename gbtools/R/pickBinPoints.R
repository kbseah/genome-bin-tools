#' Define polygon for genome bin on plot
#' 
#' Wrapper for locator() and polygon() to perform interactive binning on the
#' current plot. Returns the polygon vertices which can be used in
#' get.bin.stats()
#'
#' @param num.points Number of points in polygon (integer)
#' @param draw.polygon Draw polygon as overlay on plot (logical)
#'
#' @return numerical vector containing polygon vertices coordinates
#' @keywords internal
#'
pickBinPoints <- function(num.points=6,  # How many points in polygon?
                          draw.polygon=TRUE  # Overlay polygon on plot?
                          ) {
## Wrapper for locator() and polygon() to perform interactive binning on the current
## plot. Returns the polygon vertices which can be used in get.bin.stats()
    thepoints <- locator(num.points,pch=20,type="p")
    if (draw.polygon) { polygon(thepoints) }
    return(thepoints)
}