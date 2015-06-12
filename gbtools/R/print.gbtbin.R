#' Print summary of gbt objectbin
#'
#' @param x Object of class gbtbin
#' @return data.frame summarizing contents of x
#' @seealso \code{\link{gbtbin}}, \code{\link{summary.gbtbin}}
#' @export
print.gbtbin <- function(x) {
    cat("Object of class gbtbin\n")
    cat ("\nScaffolds:\n")
    lengthdf <- data.frame(x$summary$Total_length,
                           x$summary$Num_scaffolds,
                           x$summary$Scaff_length_max,
                           x$summary$Scaff_length_min,
                           x$summary$Scaff_length_median,
                           x$summary$Scaff_length_N50)
    names(lengthdf) <- c("Total", "Scaffolds", "Max", "Min", "Median", "N50")
    print(lengthdf)
    cat("\nMarkers:\n")
    markerdf <- data.frame(x$summary$Marker_sources,
                           x$summary$Num_markers,
                           x$summary$Num_unique_markers,
                           x$summary$Num_singlecopy_markers)
    names(markerdf) <- c("Source","Total","Unique","Singlecopy")
    print(markerdf)
    cat("\nSSU markers:\n")
    print(x$summary$Num_SSUs)
    cat("\ntRNA_markers:\n")
    print(x$summary$Num_tRNAs)
    cat("\nCall:\t")
    print(x$call)
}
