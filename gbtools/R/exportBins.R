#' Export a list of gbtbin objects to external file
#'
#' Takes a list of gbtbin objects, and exports files, each corresponding
#' to a single bin, with a list of contig names in that bin
#'
#' @param x List of gbtbin objects (e.g. list(bin1,bin2,bin3))
#' @param useBinNames Use the bin names associated with the list (default TRUE)
#' @param prefix If useBinNames is FALSE, then use this prefix followed by a number for the exported file names
#'
#' @return External files corresponding to the individual bins
#' @seealso \code{\link{write.gbtbin}}
#' @seealso \code{\link{importBins}}
#' @seealso \code{\link{summaryLOB}}
#' @export

exportBins <- function (x, # List of gbtbin objects
                        useBinNames=TRUE,
                        prefix="gbtExport" # Filename prefix for exported contig lists -- used if useBinNames is FALSE 
                        ) {
    if (!is.list(x)) {
        stop("Object x must be a list of gbtbin objects")
    } else {
        if (is.null(names(x))) {
            outnames <- paste(rep(prefix,length(x)),1:length(x),sep="_")
        } else {
            if (useBinNames) {
                outnames <- names(x)
            } else {
                outnames <- paste(rep(prefix,length(x)),1:length(x),sep="_")
            }
        }
        for (i in 1:length(x)) {
            if (class(x[[i]]) != "gbtbin") { stop("Object x mustb e a list of gbtbin objects")}
            else {
                write(as.character(x[[i]]$scaff$ID),file=paste(as.character(outnames[i]),"list",sep="."))
            }
        }
    }
}