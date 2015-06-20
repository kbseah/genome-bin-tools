#' Add custom user annotations to gbt object
#'
#' Custom user annotations for each scaffold can be added to existing gbt
#' objects. The annotations should be in a data.frame, with at least column
#' "scaffold" that matches scaffold IDs in the gbt object. Pass the name of the
#' data.frame to the userTab parameter. Give a unique name for this annotation
#' to the userSource parameter.
#'
#' @param x Object of class gbt
#' @param userTab data.frame with user annotations, see Details
#' @param userSource Name for this annotation table
#' @return Object of class gbt
#' @seealso \code{\link{gbt}} \code{\link{plot.gbt}}
#' @export
userAdd <- function(x,userTab,userSource) UseMethod("userAdd")
