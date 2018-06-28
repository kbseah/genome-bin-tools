#' @export

choosebinPolygon.gbt <- function(x,  # Object of class gbt
                          slice,  # Which slices used for the plot from which points to be chosen?
                          taxon="Class",  # Deprecated - user don't change this
                          binpolygon=NA, # The polygon
                          save=FALSE,  # Save list of contigs in bin to external file?
                          file="interactive_bin.list"  # Name of file to save list of contigs in bin
                          ) {
    require(sp)
## Wrapper for picking bin interactively from GC-cov or diff-cov plot
    if (!is.numeric(slice) || length(slice) > 2) {
        cat ("gbtools ERROR: Please specify the library(-ies) used to make the plot in focus\n")
    } else {
        if (length(slice)==1) {  # Pick bin from GC-coverage plot
            X <- merge(data.frame(ID=x$scaff$ID,
                                  Ref_GC=x$scaff$Ref_GC),
                       data.frame(ID=x$covs$ID,
                                  Avg_fold=x$covs[slice[1]+1]),
                       by="ID")
            names(X) <- c("ID","Ref_GC","Avg_fold")
            inpolygon <- sp::point.in.polygon(X$Ref_GC,
                                              X$Avg_fold,
                                              binpolygon$x,
                                              binpolygon$y)
        }
        else if (length(slice)==2) {  # Pick bin from differential coverage plot
            X <- merge(data.frame(ID=x$scaff$ID,
                                  Ref_GC=x$scaff$Ref_GC),
                       data.frame(ID=x$covs$ID,
                                  Avg_fold_1=x$covs[slice[1]+1],
                                  Avg_fold_2=x$covs[slice[2]+1]),
                       by="ID")
            names(X) <- c("ID","Ref_GC","Avg_fold_1","Avg_fold_2")
            inpolygon <- sp::point.in.polygon(X$Avg_fold_1,
                                              X$Avg_fold_2,
                                              binpolygon$x,
                                              binpolygon$y)
        }
        X.subset <- X[which(inpolygon==1),]
        X.shortlist <- as.character(X.subset$ID)
        result <- gbtbin(shortlist=X.shortlist,
                         x=x,
                         slice=slice,
                         taxon=taxon,
                         points=binpolygon,
                         save=save,
                         file=file)
        result$call[[length(result$call)+1]] <- match.call()  # Record choosebin() function call
        return(result)
    }
}
