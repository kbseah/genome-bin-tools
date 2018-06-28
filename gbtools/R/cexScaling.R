#' Function for different types of plot character scaling

#' @seealso \code{\link{plot.gbt}}

cexScaling <- function (x, # Numeric vector (contig lengths)
                        type, # Scaling type (options: "length","area")
                        const # Linear scaling constant (numeric)
                        ) {
    if (type == "area") {
        result <- sqrt(x)/const
    }
    else if (type == "length") {
        result <- x / (10*const)
    }
    return(result)
}