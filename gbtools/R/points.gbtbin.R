#' Add points to plot of gbt or gbtbin object
#'
#' Add points to GC-coverage or differential coverage plots from a gbt object
#'
#' See \code{\link{plot.gbt}} for more details on gbt plots.
#'
#' @inheritParams plot.gbt
#' @param slice
#'
#' @seealso \code{\link{plot.gbt}}, \code{\link{plot.gbtbin}}
#' @export
#'
points.gbtbin <- function(x,  # Object of class gbtbin
                          col="black",  # Overlay plot points color
                          slice="default",  # Which slice to use for plotting?
                          cutoff=0,  # Min contig length to plot
                          pch=20,
                          symbolScale="area",
                          symbolScaleParam=100,
                          ...) {
    ## Defaults to the same slice used to choose the bin ###########################
    if (slice == "default") { 
        slice <- x$slice
    }
    ## Catch invalid slice values ##################################################
    if (is.na(slice) || !is.numeric(slice) || length(slice) > 2) { 
        stop("Please specify valid value for slice option for this bin")
    }
    ## Add points to GC-coverage plot ##############################################
    else if (is.numeric(slice) && length(slice)==1) {
        X <- merge(data.frame(ID=x$scaff$ID,
                              Ref_GC=x$scaff$Ref_GC,
                              Length=x$scaff$Length),
                   data.frame(ID=x$covs$ID,
                              Avg_fold=x$covs[slice[1]+1]),
                   by="ID")
        names(X) <- c("ID","Ref_GC","Length","Avg_fold")
        if (cutoff > 0) {
            X <- subset(X, Length >= cutoff)
        }
        points(x=X$Ref_GC,
               y=X$Avg_fold,
               pch=pch,
               cex=cexScaling(as.numeric(X$Length),type=symbolScale,const=symbolScaleParam),
               col=col, ...)
    }
    ## Add points to differential coverage plot #####################################
    else if (is.numeric(slice) && length(slice)==2) {
        X <- merge(data.frame(ID=x$scaff$ID,
                              Ref_GC=x$scaff$Ref_GC,
                              Length=x$scaff$Length),
                   data.frame(ID=x$covs$ID,
                              Avg_fold_1=x$covs[slice[1]+1],
                              Avg_fold_2=x$covs[slice[2]+1]),
                   by="ID")
        names(X) <- c("ID","Ref_GC","Length","Avg_fold_1","Avg_fold_2")
        if (cutoff > 0 ) {
            X <- subset(X,Length >= cutoff)
        }
        points(x=X$Avg_fold_1,
               y=X$Avg_fold_2,
               pch=pch,
               cex=cexScaling(as.numeric(X$Length),type=symbolScale,const=symbolScaleParam),
               col=col, ...)
    }
    ## Throw error message for invalid slice options ################################
    else { stop("Please specify valid value for slice option for this bin")}
}