#' Create bin objects from taxonomic affiliation of markers
#'
#' Given a gbt object with taxonomic markers, create a set of bin objects
#' corresponding to the affiliations, at a given taxonomic level (e.g. Order or
#' Family). 
#'
#' @param x gbt object
#' @param marksource Name of marker source to use for the taxon assignment
#' @param taxon Taxonomic level to perform the binning (e.g. "Class", "Family")
#' @param out Prefix for output bin names
#' @param to.list Output bins to list of gbtbin objects (default TRUE)
#'
#' @return List of gbtbin objects 
#' @seealso \code{\link{summaryLOB}}
#' @seealso \code{\link{importBins}}
#' @seealso \code{\link{exportBins}}
#' @seealso \code{\link{winnow}}
#' @export

binsFromMarkers <- function(x, # Object of class gbt
                            marksource="amphora2", # Taxonomic marker set to use
                            taxon="Order", # Taxonomic level at which to perform the binning
                            out="taxonBins", # Prefix for output bin names
                            to.list=TRUE # Output bins to list of gbtbin objects? (default TRUE)
                           ) {
    if (class(x) != "gbt") {
        stop("Object x must be of class gbtbin")
    } else {
        theout <- NULL # Initialize list for output
        thesubset <- subset(x$markTab,source==as.character(marksource))
        bintaxa <- levels(factor(thesubset[[as.character(taxon)]]))
        for (i in 1:length(bintaxa)) {
            theshortlist <- as.vector(thesubset[which(thesubset[[taxon]]==bintaxa[i]),"scaffold"])
            theout[[paste(out,i,sep="_")]] <- gbtbin(shortlist=theshortlist, x=x, slice=1)
        }
        return(theout)
    }
    
}