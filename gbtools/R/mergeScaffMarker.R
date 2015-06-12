#' Merge scaffolds table with marker table
#'
#' Merge table of scaffold statistics (output from pileup.sh in BBMap package)
#' and table of marker statistics parsed by parse_phylotype_result.pl
#'
#' @inheritParams generateLegendColors
#'
#' @return data.frame with marker data and scaffold statistics
#' @keywords internal
#' @importFrom plyr ddply
#'
mergeScaffMarker <- function(scaffold.stats,marker.list,taxon,consensus=TRUE) {
## Merge table of scaffold statistics (output from pileup.sh in BBMap package)
## and table of marker statistics parsed by parse_phylotype_result.pl
    require(plyr)
    marker.list[,"taxon"] <- marker.list[,which(names(marker.list)==taxon)]
    marker.stats <- merge(scaffold.stats,
                          marker.list,
                          by.x="ID",
                          by.y="scaffold")
    ## For scaffolds with multiple marker genes, take majority consensus of marker #
    ## taxon assignment ############################################################
    if (consensus) {    
        #scaffs.with.multi <- as.vector(names(table(marker.stats$ID)[which(table(marker.stats$ID)>1)]))
        #consensus.list <- ddply(marker.list,
        #                       .(scaffold),
        #                       function(x) levels(x$taxon)[which.max(tabulate(x$taxon))])
        consensus.list <- plyr::ddply(marker.list,
                                      .(scaffold),
                                      summarize,
                                      taxon=levels(taxon)[which.max(tabulate(taxon))])
        marker.stats <- merge(scaffold.stats,
                              consensus.list,
                              by.x="ID",
                              by.y="scaffold")
    }
    return(marker.stats)
}