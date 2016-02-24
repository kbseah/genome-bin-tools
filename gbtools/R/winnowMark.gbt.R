#' Subset a gbt or gbtbin object by marker gene taxonomy
#'
#' @param marksource Which marker set to use (default -- use the first one)
#' @param param Taxonomic level to make subset (default "Class")
#' @param value Value of the taxon to make subset (default
#'               "Gammaproteobacteria")
#' @inheritParams winnow
#' @return Object of class gbtbin
#' @seealso \code{\link{winnow}}, \code{\link{gbt}}
#' @export
#'
winnowMark.gbt <- function(x,  # Object of class gbt
                           marksource=NA, # Which marker set to use
                           param="Class",  # Which taxonomic level to choose?
                           value="Gammaproteobacteria",  # Which taxon to choose?
                           save=FALSE,  # Save list of contigs to external file?
                           file="bin_scaffolds.list"  # File to save list of contigs
                           ) {
## Winnow a gbt object by its marker table values
    if (is.na(marksource)) {
        cat("gbtools WARNING: No marksource supplied, using the first marker set by default...\n")
        marksource=levels(x$markTab$source)[1]
    }
    markTab.subset <- subset(x$markTab,source==marksource)
    scafflist <- as.character(markTab.subset$scaffold[which(markTab.subset[,which(names(markTab.subset)==param)]==value)])
    bin <- gbtbin (shortlist=scafflist,
                   x=x,
                   points=NA,
                   slice=NA,
                   save=save,
                   file=file)
    bin$call[[length(bin$call)+1]] <- match.call()  # Record function call
    return(bin)
}
