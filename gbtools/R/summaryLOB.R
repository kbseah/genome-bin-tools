#' Summarize a list of gbtbin objects
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

summaryLOB <- function (x, # List of gbtbin objects
                         marksource="amphora2" # Marker source to summarize
                         ) {
    if (length(x[[1]]$summary[[7]]) > 1) { # If there is more than one marksource
        theout <- sapply(x, function (y) {s <- which(y$summary[[7]] == marksource)
                                          c(y$summary[1:6],
                                            y$summary[[7]][s],
                                            y$summary[[8]][s],
                                            y$summary[[9]][s],
                                            y$summary[[10]][s],
                                            y$summary[11:13])
                                          }
                        )
    } else if (length(x[[1]]$summary[[7]]) == 1 ) { # If there is only one marksource
        theout <- sapply(x, function(y) y$summary) 
    }
    rownames(theout) <- c("Total_length",
                        "Num_scaffolds",
                        "Scaff_length_max",
                        "Scaff_length_min",
                        "Scaff_length_median",
                        "Scaff_length_N50",
                        "Marker_sources",
                        "Num_markers",
                        "Num_unique_markers",
                        "Num_singlecopy_markers",
                        "Num_SSUs",
                        "Num_tRNAs",
                        "Num_tRNAs_types")
    return(theout)
}

