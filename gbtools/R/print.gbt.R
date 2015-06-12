#' Print summary of gbt object
#'
#' @param x Object of class gbt
#' @return data.frame summarizing contents of x
#' @seealso \code{\link{gbt}}, \code{\link{summary.gbt}}
#' @export
print.gbt <- function(x) {
    cat("*** Object of class gbt ***\n")
    cat("\n*** Scaffolds ***\n")
    lengthdf <- data.frame(x$summary$Total_length,
                           x$summary$Num_scaffolds,
                           x$summary$Scaff_length_max,
                           x$summary$Scaff_length_min,
                           x$summary$Scaff_length_median,
                           x$summary$Scaff_length_N50)
    names(lengthdf) <- c("Total", "Scaffolds", "Max", "Min", "Median", "N50")
    print(lengthdf)
    cat("\n*** Marker counts by source ***")
    print(x$summary$Num_markers)
    cat("\n*** SSU markers ***\n")
    print(x$summary$Num_SSU)
    cat("\n*** tRNA markers ***\n")
    print(x$summary$Num_tRNAs)
    cat("\n*** User-supplied variables ***\n")
    print(x$userSource)
    cat("\n*** Function call history ***\n")
    print(x$call)
}