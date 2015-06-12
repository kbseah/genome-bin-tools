#' Generates colors for plot legends when coloring by markers
#'
#' @param scaffold.stats scaff data.frame from a gbt object
#' @param marker.list List of markers from markTab data.frame of gbt object
#' @param taxon Taxon names for which to generate colors
#' @param consensus If scaffold has more than one marker, take majority-rule
#'                    consensus taxonomic assignment? (logical)
#' @return data.frame with color assignments for each scaffold
#' @keywords internal 
generateLegendColors <- function(scaffold.stats,  # Same params as generatePlotColors()
                                 marker.list,
                                 taxon,
                                 consensus) {
## Generates colors for plot legends when coloring by markers
    ## Some table merging to have points to plot for the markers ###################
    marker.stats <- mergeScaffMarker(scaffold.stats,
                                     marker.list,
                                     taxon,
                                     consensus) 
    marker.list[,"taxon"] <- marker.list[,which(names(marker.list)==taxon)]
    ## Count how many taxa are only supported by one marker gene ###################
    singleton.taxa <- names(table(marker.list$taxon)[which(table(marker.list$taxon)==1)])
    ## Which taxon has the most marker genes? ######################################
    top.taxon <- names(table(marker.list$taxon)[which.max(table(marker.list$taxon))]) 
    taxnames <- names(table(marker.stats$taxon))  # Names of taxa
    ## Vector to hold color names ##################################################
    taxcolors <- rep("",length(names(table(marker.stats$taxon))))
    taxcolors[which(names(table(marker.stats$taxon))  # Give singleton taxa color "grey50" 
                    %in% singleton.taxa)] <- "grey50"
    numsingletons <- length(taxcolors[which(names(table(marker.stats$taxon))
                                            %in% singleton.taxa)]) # Count singleton taxa
    taxcolors[which(names(table(marker.stats$taxon))
                    ==top.taxon)] <- "red"  # Give taxon with most marker genes color "red"
    ## How many other colors do we need, given that all singletons have same color? #
    numcolors <- length(table(marker.stats$taxon)) - 1 - numsingletons
    ## Generate needed colors, from yellow to magenta, giving red a wide berth ######
    thecolors <- rainbow(numcolors,
                         start=1/6,
                         end=5/6)
    taxcolors[which(!(names(table(marker.stats$taxon)) %in% singleton.taxa)
                    & names(table(marker.stats$taxon))!=top.taxon)] <- thecolors
    ## Data frame containing which colors correspond to which taxa ##################
    colorframe <- data.frame(taxon=taxnames,
                             colors=taxcolors) 
    ## Merge this by taxon into the marker.stats table for plotting (this works #####
    ## even when consensus option is called) ########################################
    marker.stats <- merge(marker.stats,
                          colorframe,
                          by="taxon")
    return(colorframe)
}