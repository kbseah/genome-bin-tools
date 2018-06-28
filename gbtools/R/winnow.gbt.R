#' @export
#'
winnow.gbt <- function (x,
                        gc=c(0,1),
                        len=c(0,Inf),
                        covmin=0,
                        covmax=Inf,
                        slice=1,
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
            stop("Lengths of covmin, covmax, and slice parameters do not match")
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
    bin$call[[length(bin$call)+1]] <- match.call()  # Record function call 
    return(bin)
}
