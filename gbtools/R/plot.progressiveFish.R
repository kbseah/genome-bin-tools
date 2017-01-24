plot.progressiveFish <- function(x,
                                 parent,
                                 plotType="gccov",
                                 slice=1,
                                 cutoff=1000,
                                 legend=TRUE,
                                 log="default") {
    if(class(parent) != "gbt") {
        cat("gbtools ERROR: parent must be of class gbt")
    } else {
        if (plotType == "gccov") {
            palette <- heat.colors(length(x$binList))
            plot.gbt(parent,
                     marker=F,
                     slice=slice,
                     cutoff=cutoff,
                     log=log)
            for (i in length(x$binList):1) {
                points.gbtbin(x$binList[[i]], 
                              col=palette[i], 
                              slice=slice)
            }
            if (legend==TRUE) {
                legend(x="topright",
                      legend=0:(length(x$binList)-1),
                      fill=palette)
            }
        } else if (plotType == "summary") {
            par(mfrow=c(3,1))
            plot(x=x$summary$Fishing_depth,
                 y=x$summary$Total_length,
                 pch=20,
                 col="black",
                 ylab="Total length (bp)",
                 xlab="Fishing Depth")
            points(x=x$summary$Fishing_depth,
                   y=x$summary$Total_length,
                   type="l",
                   col="black")
            plot(x=x$summary$Fishing_depth,
                   y=x$summary$Num_scaffolds,
                   col="blue",
                   pch=20,
                   ylab="Number scaffolds",
                   xlab="Fishing Depth")
            points(x=x$summary$Fishing_depth,
                   y=x$summary$Num_scaffolds,
                   type="l",
                   col="blue")
            plot(x=x$summary$Fishing_depth,
                 y=x$summary$Num_SSUs,
                 pch=20,
                 col="green",
                 ylab="Number SSUs",
                 xlab="Fishing Depth")
            points(x=x$summary$Fishing_depth,
                   y=x$summary$Num_SSUs,
                   type="l",
                   col="green")
        }
    }
}