#' Print summary of gbtbin object
#'
#' For a gbtbin object, summary() differs from print() in displaying more
#' information: if available, table of counts per marker gene (for each marker
#' source), and table of counts per tRNA type.
#'
#' @param x Object of class gbtbin
#' @return data.frame summarizing contents of x
#' @seealso \code{\link{gbtbin}}, \code{\link{print.gbtbin}}
#' @export
summary.gbtbin <- function (x) {
    print(x)  # Print the standard summary
    ## Show the marker tables
    cat ("\nPolygon for choosebin (if applicable):\n")
    print(x$points)
    cat ("\nTable(s) of marker genes\n")
    print(x$marker.table)
    cat ("\nTable of tRNA genes\n")
    print(x$tRNA.table)
}