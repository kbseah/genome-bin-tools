#' @export

winnowMark.gbt <- function(x,  # Object of class gbt
                           marksource, # Which marker set to use
                           param,  # Which taxonomic level to choose?
                           value,  # Which taxon to choose?
                           save,  # Save list of contigs to external file?
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
