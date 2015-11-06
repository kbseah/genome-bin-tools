#' Compare two sets of bins
#'
#' Given two lists of gbtbin objects, compare them in terms of their shared
#' contigs and output a table.
#'
#' @param x List of gbtbin objects (either list(bin1,bin2) or mget(vector_of_bin_names)
#' @param y Second list of gbtbin objects
#' @param binNames.x Vector of names for each bin in list x
#' @param binNames.y Vector of names for each bin in list y
#' @param weight Compare by number of bases shared (default TRUE), else by number of contigs
#' @param by Which list of bins to compute fraction of overlap: "x" (default), "y", or "raw" (raw counts)
#'
#' @return Matrix reporting quantity of overlap between two sets of bins
#' @seealso \code{\link{mergeOverlapBins}}
#' @seealso \code{\link{importBins}}
#' @seealso \code{\link{multiBinPlot}}
#' @export
#'


tabOverlapBins <- function(x, # First list of gbtbin objects 
                           y, # Second list of gbtbin objects
                           binNames.x="", # list of names for x. Default -- number sequentially
                           binNames.y="", # list of names for y
                           weight=TRUE, # weight overlaps by number of bases if true (default), else by number of contigs
                           by="x" # Output type: "x" (default), "y", or "raw" -- by "x" means 'as fraction of rows', and by "y" means 'as fraction of columns'
                           ) {
    # Check that inputs are lists of gbtbin objects
    if (!is.list(x) | !is.list(y)) {
        cat ("gbtools ERROR: The input variables x and y must be character vectors giving names of gbtbin objects \n")
    }
    else {
        if (binNames.x=="" || binNames.y=="") {
            binNames.x <- 1:length(x)
            binNames.y <- 1:length(y)
            binNames.x <- sapply(binNames.x,
                                 function(X) paste("x",X,sep=""))
            binNames.y <- sapply(binNames.y,
                                 function(X) paste("y",X,sep=""))
        }
        if (!weight) {
            counts <- NULL # Initialize the vector of counts x in y
            countsA <- NULL # Initialize vector of fractional counts x in y divided by total of x
            countsB <- NULL # Initialize vector of fractional counts x in y divided by total of y
            for (i in 1:length(x)) {
                for (j in 1:length(y)) {
                    curr.count <- length(which(x[[i]]$scaff$ID %in% y[[j]]$scaff$ID))
                    counts <- c(counts, curr.count) # Append count of overlaps
                    countsA <- c(countsA, curr.count / length(x[[i]]$scaff$ID)) # As fraction of bin x
                    countsB <- c(countsB, curr.count / length(y[[j]]$scaff$ID)) # As fraction of bin y
                    
                }
            }
            counts.matrix <- matrix(counts, byrow=TRUE, nrow=length(x), dimnames=list(binNames.x,binNames.y))
            countsA.matrix <- matrix(countsA, byrow=TRUE, nrow=length(x), dimnames=list(binNames.x,binNames.y))
            countsB.matrix <- matrix(countsB, byrow=TRUE, nrow=length(x), dimnames=list(binNames.x,binNames.y))
            if (by=="x") {
                return (countsA.matrix)
            } else if (by=="y") {
                return (countsB.matrix)
            } else if (by=="raw") {
                return (counts.matrix)
            } else {
                cat ("gbtools ERROR: Invalid value for \"by\" -- should be \"x\", \"y\", or \"raw\" \n" )
            }
        } else {
            lens <- NULL # Initialize the vector of bases x in y
            lensA <- NULL # Initialize the vector of fractional bases x in y divided by total of x
            lensB <- NULL
            for (i in 1:length(x)) {
                for (j in 1:length(y)) {
                    curr.len <- sum(x[[i]]$scaff$Length[which(x[[i]]$scaff$ID %in% y[[j]]$scaff$ID)])
                    lens <- c(lens, curr.len)
                    lensA <- c(lensA, curr.len / sum(x[[i]]$scaff$Length))
                    lensB <- c(lensB, curr.len / sum(y[[j]]$scaff$Length))
                    
                }
            }
            lens.matrix <- matrix(lens,byrow=TRUE,nrow=length(x),dimnames=list(binNames.x,binNames.y))
            lensA.matrix <- matrix(lensA, byrow=TRUE,nrow=length(x),dimnames=list(binNames.x,binNames.y))
            lensB.matrix <- matrix(lensB, byrow=TRUE,nrow=length(x),dimnames=list(binNames.x,binNames.y))
            if (by=="x") {
                return (lensA.matrix)
            } else if (by=="y") {
                return (lensB.matrix) 
            } else if (by=="raw") {
                return(lens.matrix)
            } else {
                cat ("gbtools ERROR: Invalid value for \"by\" -- should be \"x\", \"y\", or \"raw\" \n" )
            }
        }
    }
}

#tabOverlapBins(x=list(bin77,bin96),y=list(overlapbin96,overlapbin77))