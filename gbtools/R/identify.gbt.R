#' Identify points in gbtools plot
#'
#' Click on GC-coverage plot or differential coverage plot to identify contigs.
#' In cluttered plots it may not be very accurate! Will write the contig ID
#' as a label overlay on the plot.
#'
#' @return Object of class gbtbin containing identified contigs
#' 
#' @seealso \code{\link{plot.gbt}}
#' @seealso \code{\link{choosebin}}
#'
#' @export

identify.gbt <- function(d,
#' @param d Object of class gbt or gbtbin, in plot
                         slice=1,
#' @param slice Which sample data was plotted? (Default: 1)
                         ...
#' @param ... Further arguments passed to identify.default()
                         ) {
    # Catch invalid "slice" parameters
    if (is.na(slice) || !is.numeric(slice)) {
        stop("Please specify valid value for slice parameter")
    } else {
        # Data frame for GC-coverage plots
        if (length (slice) == 1) {
            X <- merge(data.frame(ID=d$scaff$ID,
                                  Ref_GC=d$scaff$Ref_GC,
                                  Length=d$scaff$Length,
                                  Avg_fold=d$scaff$Avg_fold,
                                  xVals=d$scaff$Ref_GC),
                       data.frame(ID=d$covs$ID,
                                  yVals=d$covs[[slice[1]+1]]
                                  ),
                       by="ID")
            names(X) <- c("ID","Ref_GC","Length","Avg_fold","xVals","yVals")
        }
        # Data frame for differential coverage plots
        else if (length(slice) == 2) {
            X <- merge (data.frame (ID=d$scaff$ID,
                                    Ref_GC=d$scaff$Ref_GC,
                                    Length=d$scaff$Length,
                                    Avg_fold=d$scaff$Avg_fold),
                        data.frame (ID=d$covs$ID,
                                    xVals=d$covs[[slice[1]+1]],
                                    yVals=d$covs[[slice[2]+1]]),
                        by="ID")
            names(X) <- c("ID","Ref_GC","Length","Avg_fold","xVals","yVals")
        }
        # Implement the identify parameter...
        shortlist <- identify(x=X$xVals,
                              y=X$yVals,
                              labels=X$ID,
                              ...
                              )
        shortlist.contigs <- X$ID[shortlist]
        return(gbtbin(as.character(shortlist.contigs),d,slice=slice))
        #return(shortlist)
    }
}