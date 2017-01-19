#' Generates colors for marker gene phylotypes in plot by cumulative weight
#'
#' For each taxon, calculates the total contigLength*coverage, and assigns
#' colors for the top taxa which in total account for more than a specified
#' minimum weight.
#'
#' @param scaffold.stats Scaff table from gbt object
#' @param marker.list markTab table from gbt object
#' @param taxon Taxonomic level to do coloring
#' @param consensus Logical - if taxon assignments conflict, take consensus?
#' @param weightCutoff Cutoff quantile for contig length*coverage weight (between 0 and 1)
#' @param markCustomPalette (Optional) data.frame with taxon names (column 1) and specific colors (column 2) to assign to them
#' @return data.frame with color assignments for each scaffold
#' @keywords internal
#'
generatePlotColors2 <- function(scaffold.stats,  # scaff table from gbt object
                               marker.list,  # markTab table from gbt object
                               taxon,  # Taxonomic level to do the coloring
                               consensus,  # Logical-if taxon assgs conflict, take consens?
                               weightCutoff, # Cutoff quantile for assigning colors (between 0 and 1)
                               markCustomPalette # Custom palette for specific taxa
                               ) {           # This took a very long time to get it right
## Generates colors for marker gene phylotypes in plot
    ## Merge tables to have points to plot for the markers #########################
    marker.stats <- mergeScaffMarker(scaffold.stats,
                                     marker.list,
                                     taxon,
                                     consensus)
    ## Calculate total weight for each taxon #######################################
    taxon.agg <- aggregate(marker.stats$Length*marker.stats$Avg_fold,
                           by=list(marker.stats$taxon),
                           FUN=sum
                           )
    total.weight <- sum(taxon.agg$x)
    taxon.agg.order <- taxon.agg[order(taxon.agg$x,decreasing=TRUE),] # Sort descending
    ## Assign colors to each taxon #################################################
    if (is.data.frame(markCustomPalette)) {
        # If a custom palette is specified, use those colors
        colnames(markCustomPalette) <- c("taxon","colors")
        # Vector of all taxon names
        taxon.vector <- as.character(taxon.agg.order[,1])
        # Check that all taxon names in markCustomPalette are in taxon.vector
        taxonIgnore <- markCustomPalette$taxon[which(!markCustomPalette$taxon %in% taxon.vector)]
        if (length(taxonIgnore) > 0) {
            cat("gbtools WARNING: Some taxa in custom palette not found in gbt object\n")
        }
        # Which taxa are not in the custom palette should be in grey
        taxon.grey <- taxon.vector[which(!taxon.vector %in% markCustomPalette[,1])]
        repgrey <- rep("grey50",length(taxon.grey))
        greyframe <- data.frame(taxon=taxon.grey,
                                colors=repgrey)
        # Combine the custom palette with the taxa that are greyed out
        colorframe <- rbind (markCustomPalette,greyframe)
        # Merge marker.stats with color values
        marker.stats <- merge(marker.stats,
                              colorframe,
                              by="taxon")
        # Return marker.stats
        return(marker.stats)
    } else {
        # Otherwise, by default assign rainbow colors to taxa in order of weights
        if (weightCutoff > 1 || weightCutoff < 0) { # Catch errors for weightCutoff
            cat ("gbtools ERROR: weightCutoff parameter must be between 0 and 1\n")
        }
        else {
            # Count taxa which have the highest weight, until weightCutoff
            numAboveCutoff <- length(
                                     which(
                                           cumsum(taxon.agg.order$x) <= weightCutoff*total.weight
                                           )
                                     )
            numBelowCutoff <- length(
                                     which(
                                           cumsum(taxon.agg.order$x) > weightCutoff*total.weight
                                           )
                                     )
            # Generate colors from red to violet for taxa above cutoff
            thecolors <- rainbow (numAboveCutoff,
                                  start=0,
                                  end=3/4)
            # Everything else is colored grey
            repgrey <- rep ("grey50", numBelowCutoff)
            # Combine the two vectors
            thecolors <- c(thecolors, repgrey)
            colorframe <- data.frame(taxon=taxon.agg.order[,1],
                                     colors=thecolors)
            # Merge into marker.stats df
            marker.stats <- merge(marker.stats,
                                  colorframe,
                                  by="taxon")
            # Return table
            return(marker.stats)
        }
    }
}
