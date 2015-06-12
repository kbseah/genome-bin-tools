#' Subset a gbt or gbtbin object by GC, length, or coverage cutoffs
#'
#' @param x Object of class gbt or gbtbin
#' @param gc Vector giving min and max GC values (default c(0,Inf))
#' @param len Vector giving min and max contig lengths (default c(0,Inf))
#' @param covmin Vector giving minimum coverage values per sample (default NA)
#' @param covmax Vector giving maximum coverage values per sample (default NA)
#' @param slice Vector of sample numbers for the coverage cutoffs (default NA)
#' @param save Save list of contigs to external file? (logical, defautl FALSE)
#' @param file File name for export of contig list.
#' @return Object of class gbtbin
#' @seealso \code{\link{winnowMark}}, \code{\link{gbt}}
#' @export
#'
winnow <- function (x, gc, len, covmin, covmax, slice, save, file) UseMethod ("winnow")
winnow.gbt <- function (x,
                        gc=c(0,1),
                        len=c(0,Inf),
                        covmin=NA,
                        covmax=NA,
                        slice=NA,
                        save=FALSE,
                        file="bin_scaffolds.list") {
## "Winnow" a gbt object by GC%, Length, and/or coverage cutoffs
    ## Shortlist scaffolds that match GC and Length criteria #######################
    scafflist <- as.character(x$scaff$ID[which(x$scaff$Ref_GC > gc[1]
                                               & x$scaff$Ref_GC < gc[2]
                                               & x$scaff$Length > len[1]
                                               & x$scaff$Length < len[2])])
    ## Check if coverage cutoffs given, with relevant slices #######################
    if ( is.numeric(slice) ) {
        if (length(covmin)==length(covmax)
            && length(covmin) == length(slice)) {  # Check that values match
            covslist <- as.character(x$covs$ID[which(x$covs[slice[1]+1] > covmin[1]
                                                     & x$covs[slice[1]+1] < covmax[1])])
            if (length(slice)> 1) {
                for (i in 2:length(slice)) {
                    covslist2 <- as.character(x$covs$ID[which(x$covs[slice[i]+1]> covmin[i]
                                                              & x$covs[slice[i]+1]< covmax[i])])
                    covslist <- intersect(covslist,covslist2)
                }
            }
            scafflist <- intersect(scafflist, covslist)  # Update scaffolds shortlist
        } else {
            cat("gbtools ERROR: Lengths of covmin, covmax, and slice parameters
                do not match\n")
        }
    }
    ## Package and return result ###################################################
    bin <- gbtbin(shortlist=scafflist,
                  x=x,
                  slice=NA,
                  taxon=taxon,
                  points=NA,
                  save=save,
                  file=file)
    bin$call <- match.call()  # Record function call that produced this winnow bin
    return(bin)
}
winnow.gbtbin <- winnow.gbt  # Inherit behavior of gbt method