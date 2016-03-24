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
    ## Show summary of coverage
    cat ("\n***Summary of coverage values ***\n")
    weighted <- as.numeric(
                           unlist(
                                  apply(x[["scaff"]],
                                        1,
                                        function(y) rep(y[["Avg_fold"]],
                                                        y[["Length"]]
                                                        )
                                        )
                                  )
                           )
    covtable <- t(data.frame(quantile(x$scaff$Avg_fold,probs=c(0,0.25,0.5,0.75,1)),
                             quantile(weighted,probs=c(0,0.25,0.5,0.75,1))
                             )
                  )
    rownames(covtable) <- c("unweighted","weighted")
    colnames(covtable) <- c("min","q25","median","q75","max")
    print(covtable)
    ## Print polygon for choosebin
    cat ("\n*** Polygon for choosebin (if applicable) ***\n")
    print(x$points)
    ## Show the marker tables
    cat ("\n*** Table(s) of marker genes ***\n")
    print(x$marker.table)
    cat ("\n*** Table of tRNA genes ***\n")
    print(x$tRNA.table)
}