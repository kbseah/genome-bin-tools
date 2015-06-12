#' Subset a gbt or gbtbin object by marker gene taxonomy
#'
#' @param param Taxonomic level to make subset (default "Class")
#' @param value Value of the taxon to make subset (default
#'               "Gammaproteobacteria")
#' @inheritParams winnow
#' @return Object of class gbtbin
#' @seealso \code{\link{winnow}}, \code{\link{gbt}}
#' @export
#'
winnowMark <- function(x,param,value,save,file) UseMethod("winnowMark")
winnowMark.gbt <- function(x,  # Object of class gbt
                           param="Class",  # Which taxonomic level to choose?
                           value="Gammaproteobacteria",  # Which taxon to choose?
                           save=FALSE,  # Save list of contigs to external file?
                           file="bin_scaffolds.list"  # File to save list of contigs
                           ) {
## Winnow a gbt object by its marker table values
    scafflist <- as.character(x$markTab$scaffold[which(x$markTab[,which(names(x$markTab)
                                                                        ==param)]
                                                       ==value)])
    bin <- gbtbin (shortlist=scafflist,
                   x=x,
                   points=NA,
                   slice=NA,
                   save=save,
                   file=file)
    bin$call[[length(bin$call)+1]] <- match.call()  # Record function call
    return(bin)
}
winnowMark.gbtbin <- winnowMark.gbt