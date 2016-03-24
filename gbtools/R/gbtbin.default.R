#' Create gbtbin from gbt object
#'
#' Creates object of class gbtbin from existing gbt object, representing a
#' genome bin defined from a metagenome assembly.
#'
#' This function is called by \code{\link{choosebin}}, but can be called
#' directly by the user, given a shortlist of contigs in the bin.
#'
#' @param shortlist List of contigs in the new bin (required)
#' @param x Object of class gbt
#' @param slice Numeric vector describing which coverage values were used to
#'               draw the plot from which the bin was selected by choosebin.
#'               Use NA if not calling this function from choosebin.
#' @param taxon For internal use
#' @param points Numeric vector describing points used to select the bin from
#'                a plot by choosebin. Use NA if not calling this function
#'                from choosebin.
#' @param save Save list of contigs to external file? (logical, default FALSE)
#' @param file File name to save contigs, if save=TRUE.
#'
#' @return Object of class gbtbin
#' @seealso \code{\link{gbt}}, \code{\link{choosebin}}
#' @export
gbtbin.default <- function(shortlist,  # Character vector, contigs to extract from gbt object
                           x,  # Object of class gbt
                           slice,  # Which slice used by choosebin()?
                           taxon,  # Deprecated - user don't change
                           points=NA,  # Number of points in polygon for choosebin()
                           save=FALSE,  # Save contig list to external file?
                           file="interactive_bin.list"  # File name to save contig list
                           ) {
    scaff.subset <- subset(x$scaff, ID%in%shortlist)
    covs.subset <- subset(x$covs, ID%in%shortlist)
    markTab.subset <- NA
    marksource <- x$markSource
    ssuTab.subset <- NA
    trnaTab.subset <- NA
    thecall <- x$call
    userTab <- list()
    userSource <- x$userSource
    ## Summary statistics initialize ##########################################
    bin.nummarkers <- NA
    bin.uniqmarkers <- NA
    bin.numtRNAs <- NA
    bin.uniqtRNAs <- NA
    bin.numSSUs <- NA
    bin.singlemarkers <- NA
    marker.tab <- list()
    tRNAs.tab <- NA
    ## Take subset of markTab #################################################
    if (is.data.frame(x$markTab)) {
        markTab.subset <- subset(x$markTab,scaffold%in% shortlist)
        for (j in 1:length(marksource)) {
            markTab.subsubset <- droplevels(subset(markTab.subset,source==marksource[j]))
            bin.nummarkers[j] <- dim(markTab.subsubset)[1]
            marker.tab[[j]] <- table(markTab.subsubset$gene)
            bin.uniqmarkers[j] <- length(which(marker.tab[[j]] > 0))
            bin.singlemarkers[j] <- length(which(marker.tab[[j]] ==1))
        }
    }
    ## Take subset of ssuTab ##################################################
    if (is.data.frame(x$ssuTab)) {
        ssuTab.subset <- subset(x$ssuTab,scaffold%in%shortlist)
        bin.numSSUs <- dim(ssuTab.subset)[1]
    }
    ## Take subset of trnaTab #################################################
    if (is.data.frame(x$trnaTab)) {
        trnaTab.subset <- subset(x$trnaTab,scaffold%in% shortlist)
        bin.numtRNAs <- dim(trnaTab.subset)[1]
        tRNAs.tab <- table (trnaTab.subset$tRNA_type)
        bin.uniqtRNAs <- length(which(tRNAs.tab > 0))
    }
    ## Take subset of userTab #################################################
    if (length(x$userTab) > 0) {
        for (i in 1:length(x$userTab)) {
            userTabSubset <- subset(x$userTab[[i]],scaffold %in% shortlist)
            userTab[[i]] <- userTabSubset
        }
    }
    ## Summary statistics #####################################################
    bin.length <- sum(scaff.subset$Length)
    bin.numscaffolds <- dim(scaff.subset)[1]
    bin.summary <- list(Total_length=bin.length,
                        Num_scaffolds=bin.numscaffolds,
                        Scaff_length_max=max(scaff.subset$Length),
                        Scaff_length_min=min(scaff.subset$Length),
                        Scaff_length_median=median(scaff.subset$Length),
                        Scaff_length_N50=getN50(scaff.subset$Length),
                        Marker_sources=marksource,
                        Num_markers=bin.nummarkers,
                        Num_unique_markers=bin.uniqmarkers,
                        Num_singlecopy_markers=bin.singlemarkers,
                        Num_SSUs=bin.numSSUs,
                        Num_tRNAs=bin.numtRNAs,
                        Num_tRNAs_types=bin.uniqtRNAs)
    
    ## Write to file, if save option is used ##################################
    if (save) {
        write(as.vector(scaff.subset$ID),file=file)
    }
    
    ## Package and return result ##############################################
    result <- list(scaff=scaff.subset,
                   covs=covs.subset,
                   markTab=markTab.subset,
                   markSource=marksource,
                   ssuTab=ssuTab.subset,
                   trnaTab=trnaTab.subset,
                   userTab=userTab,
                   userSource=userSource,
                   summary=bin.summary,
                   marker.table=marker.tab,
                   tRNA.table=tRNAs.tab,
                   points=points,
                   slice=slice,
                   call=thecall)
    class(result) <- "gbtbin"
    return(result)
}
