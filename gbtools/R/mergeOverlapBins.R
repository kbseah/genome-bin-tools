#' Merge two sets of bins
#'
#' Given two lists of gbtbin objects, compare them in terms of their shared
#' contigs and merge bins which exceed a certain overlap threshold.
#'
#' @inheritParams tabOverlapBins
#' @param by Condition by fraction in y (default) or x
#' @param mergeto Merge bins in list x to list y ("x" - default) or vice-versa
#'        ("y")
#' @param threshold Threshold for fraction that must be shared before bins are
#'        merged (default 0.8)
#' @param out Name prefix of output bins (default "mergedBin")
#'
#' @return List of gbtbin objects containing the merged bins (if any)
#' @seealso \code{\link{tabOverlapBins}}
#' @seealso \code{\link{importBins}}
#' @seealso \code{\link{multiBinPlot}}
#'
#' @export

mergeOverlapBins <- function(x, # List of gbtbin objects
                             y, # Second list of gbtbin objects
                             binNames.x="", # Names for gbtbin objects
                             binNames.y="", # Names of gbtbin objects
                             weight=TRUE, # Count overlaps by number of bases (default), else number of contigs
                             by="y", # Condition by fraction in y (default) or x
                             mergeto="x", # Merge y into x (default) or x into y; usually the opposite to "by"
                             threshold="0.8", # Threshold of shared fraction for merging bins
                             out="mergedBin" # Name prefix for output bins
                             ) {
    # Check that inputs are lists of gbtbin objects
    if (!is.list(x) || !is.list(y)) {
        stop("Input should be lists of gbtbin objects")
    } else {
        if (binNames.x=="" || binNames.y=="") {
            binNames.x <- 1:length(x)
            binNames.y <- 1:length(y)
            binNames.x <- sapply(binNames.x,
                                 function(X) paste("x",X,sep=""))
            binNames.y <- sapply(binNames.y,
                                 function(X) paste("y",X,sep=""))
        }
        # Calculate the overlap table
        overlapTab <- tabOverlapBins(x=x,y=y,
                                     binNames.x=binNames.x,
                                     binNames.y=binNames.y,
                                     weight=weight,
                                     by=by)
        # Initialize the list of merged bins
        outList <- NULL
        counter <- 1 # Counter to name the merged bins
        # If merging by "x"
        if (mergeto=="x") {
            for (i in 1:length(x)) {
                hits <- which(overlapTab[i,] > threshold)
                if ( length(hits) > 0 ) {
                    # merge the first y object into the x object
                    newbin <- add(x[[i]],y[[hits[1]]])
                    # If more than one y object exceeds threshold, merge those in sequentially
                    if ( length(hits) > 1) {
                        for (j in 2:length(hits)) {
                            newbin <- add(newbin,y[[hits[j]]])
                        }
                    }
                    # Report which bins were merged
                    cat ("Merged:",
                         as.character(binNames.x[i]),
                         as.character(binNames.y[hits]),
                         "\n",
                         sep=" ")
                    # Add the new bin to list of merged bins
                    outList[[paste(as.character(out),as.character(counter),sep="")]] <- newbin
                    counter <- counter+1
                }
            }
        } else if (mergeto=="y") {   # If merging by "y"
            for (i in 1:length(y)) {
                hits <- which(overlapTab[,i] > threshold)
                if ( length(hits) > 0 ) {
                    # merge the first y object into the x object
                    newbin <- add(y[[i]],x[[hits[1]]])
                    # If more than one y object exceeds threshold, merge those in sequentially
                    if ( length(hits) > 1) {
                        for (j in 2:length(hits)) {
                            newbin <- add(newbin,x[[hits[j]]])
                        }
                    }
                    # Report which bins were merged
                    cat ("Merged:",
                         as.character(binNames.y[i]),
                         as.character(binNames.x[hits]),
                         "\n",
                         sep=" ")
                    # Add the new bin to list of merged bins
                    outList[[paste(as.character(out),as.character(counter),sep="")]] <- newbin
                    counter <- counter+1
                }
            }
        }
        # Return list of gbtbin objects
        return(outList)
    }
}