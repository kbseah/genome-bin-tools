#' Generic operation for merging and subsetting two gbtbin objects
#' 
#' @param x1 Object of class gbtbin
#' @param x2 Object of class gbtbin
#' @param shortlist Vector of contig IDs to make new bin
#' @return Object of class gbtbin
#' @keywords internal

setOperation.gbtbin <- function(x1,x2,shortlist) {
## Generic operation for merging and subsetting two gbtbin objects
    ## Merge the two scaff tables ##################################################
    scaff.add <- unique (rbind(x1$scaff,x2$scaff))
    ## Catch special value for shortlist ###########################################
    if (shortlist[1] == "all") {
        shortlist <- as.character(scaff.add$ID)
    }
    ## Subset scaff and coverage tables ############################################
    scaff.add <- subset (scaff.add, ID %in% shortlist)
    covs.add <- unique (rbind(x1$covs,x2$covs))
    covs.add <- subset(covs.add, ID %in% shortlist)
    ## Initialize marker, ssu, and tRNA tables #####################################
    markTab.add <- NA
    ssuTab.add <- NA
    trnaTab.add <- NA
    marksource <- x1$markSource  # markSource should be identical in x1 and x2 -- TODO: Implement a check?
    ## Initialize variables for summary stats ######################################
    bin.nummarkers <- NA    # In case the marker.list is not supplied
    bin.uniqmarkers <- NA
    bin.numtRNAs <- NA      # Likewise for number of tRNAs
    bin.uniqtRNAs <- NA
    bin.numSSUs <- NA
    bin.singlemarkers <- NA
    marker.tab <- list()  # Make empty list object (NA will throw error later)
    tRNAs.tab <- NA
    ## Combine marker tables #######################################################
    #if (!is.na(x1$markTab) || !is.na(x2$markTab)) {
    if (is.data.frame(x1$markTab) && is.data.frame(x2$markTab)) {
        if (dim(x1$markTab)[1] > 0 || dim(x2$markTab)[1] > 0) { 
            markTab.add <- unique (rbind(x1$markTab,x2$markTab))
            markTab.add <- subset(markTab.add, scaffold %in% shortlist)
            for (j in 1:length(marksource)) {
                markTab.subsubset <- subset(markTab.add,source==marksource[j])
                bin.nummarkers[j] <- dim(markTab.subsubset)[1]
                marker.tab[[j]] <- table(markTab.subsubset$gene)
                bin.uniqmarkers[j] <- length(which(marker.tab[[j]] > 0))
                bin.singlemarkers[j] <- length(which(marker.tab[[j]] ==1))
            }
        }
    }
    ## Combine SSU Tables ##########################################################
    #if (!is.na(x1$ssuTab) || !is.na(x2$ssuTab)) {
    if (is.data.frame(x1$ssuTab) && is.data.frame(x2$ssuTab)) {
        if ( dim(x1$ssuTab)[1]>0 || dim(x2$ssuTab)[1]>0 ) {
            ssuTab.add <- unique (rbind(x1$ssuTab,x2$ssuTab))
            ssuTab.add <-subset(ssuTab.add, scaffold %in% shortlist)
            bin.numSSUs <- dim(ssuTab.add)[1]
        }
    }
    ## Combine trna Tables #########################################################
    #if (!is.na(x1$trnaTab) || !is.na(x2$trnaTab)) {
    if (is.data.frame(x1$trnaTab) && is.data.frame(x2$trnaTab)) {
        if ( dim(x1$trnaTab)[1] > 0 || dim(x2$trnaTab)[1] > 0 ) {
            trnaTab.add <- unique (rbind(x1$trnaTab,x2$trnaTab))
            trnaTab.add <- subset(trnaTab.add, scaffold %in% shortlist)
            bin.numtRNAs <- dim(trnaTab.add)[1]
            tRNAs.tab <- table(trnaTab.add$tRNA_type)
            bin.uniqtRNAs <- length(which(tRNAs.tab > 0))
        }
    }
    ## Summary statistics ##########################################################
    bin.length <- sum(scaff.add$Length)
    bin.numscaffolds <- dim(scaff.add)[1]
    bin.summary <- list (Total_length=bin.length,
                         Num_scaffolds=bin.numscaffolds,
                         Scaff_length_max=max(scaff.add$Length),
                         Scaff_length_min=min(scaff.add$Length),
                         Scaff_length_median=median(scaff.add$Length),
                         Scaff_length_N50=getN50(scaff.add$Length),
                         Marker_sources=marksource,
                         Num_markers=bin.nummarkers,
                         Num_unique_markers=bin.uniqmarkers,
                         Num_singlecopy_markers=bin.singlemarkers,
                         Num_SSUs=bin.numSSUs,
                         Num_tRNAs=bin.numtRNAs,
                         Num_tRNAs_types=bin.uniqtRNAs)
    ## Package and return results ##################################################
    result <- list (scaff=scaff.add,
                    covs=covs.add,
                    markTab=markTab.add,
                    markSource=marksource,
                    ssuTab=ssuTab.add,
                    trnaTab=trnaTab.add,
                    summary=bin.summary,
                    marker.table=marker.tab,
                    tRNA.table=tRNAs.tab,
                    points=NA,
                    slice=NA)
    class(result) <- "gbtbin"
    return(result)
}
