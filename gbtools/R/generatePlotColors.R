#' Generates colors for marker gene phylotypes in plot
#'
#' @inheritParams generateLegendColors
#' @return data.frame with color assignments for each scaffold
#' @keywords internal
#'
generatePlotColors <- function(scaffold.stats,  # scaff table from gbt object
                               marker.list,  # markTab table from gbt object
                               taxon,  # Taxonomic level to do the coloring
                               consensus  # Logical-if taxon assgs conflict, take consens?
                               ) {           # This took a very long time to get it right
## Generates colors for marker gene phylotypes in plot
    ## Merge tables to have points to plot for the markers #########################
    marker.stats <- mergeScaffMarker(scaffold.stats,
                                     marker.list,
                                     taxon,
                                     consensus) 
    marker.list[,"taxon"] <- marker.list[,which(names(marker.list)==taxon)]
    ## Count taxa supported only by one marker gene ################################
    singleton.taxa <- names(table(marker.list$taxon)[which(table(marker.list$taxon)==1)])
    ## Which taxon has the most marker genes? ######################################
    top.taxon <- names(table(marker.list$taxon)[which.max(table(marker.list$taxon))])
        # Use which.max() because it breaks ties. Otherwise all genomes with same
        # number of marker genes will have same color!
        # Important: Identification of singleton taxa uses the original marker.list
        #            because after "consensus", each scaffold has only one taxon
        #            assignment and scaffolds with >1 marker will be undercounted
    ## For plot colors - identify singleton taxa and the taxon with the highest ####
    ## marker counts, and assign them special colors ###############################
    taxnames <- names(table(marker.stats$taxon))  # Names of taxa
    taxcolors <- rep("",length(names(table(marker.stats$taxon))))  # Vector to hold color names
    taxcolors[which(names(table(marker.stats$taxon))
                    %in% singleton.taxa)] <- "grey50"  # Give singletons color "grey50"
    numsingletons <- length(taxcolors[which(names(table(marker.stats$taxon))
                                            %in% singleton.taxa)]) # Count singleton taxa
    taxcolors[which(names(table(marker.stats$taxon))
                    ==top.taxon)] <- "red"  # Give taxon with most marker genes color "red"
    ## How many other colors do we need, given that all singletons have same color?#
    numcolors <- length(table(marker.stats$taxon)) - 1 - numsingletons
    ## Generate needed colors, from yellow to magenta, giving red a wide berth #####
    thecolors <- rainbow(numcolors,
                         start=1/6,
                         end=5/6) 
    taxcolors[which(!(names(table(marker.stats$taxon))%in% singleton.taxa)
                    & names(table(marker.stats$taxon))!=top.taxon)] <- thecolors
    ## Data frame containing which colors correspond to which taxa #################
    colorframe <- data.frame(taxon=taxnames,
                             colors=taxcolors)
    ## Merge this by taxon into the marker.stats table for plotting (this works ####
    ## even when consensus option is called) #######################################
    marker.stats <- merge(marker.stats,
                          colorframe,
                          by="taxon")
    return(marker.stats)
}