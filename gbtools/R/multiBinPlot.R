#' Plot gbt object and overlay multiple bins
#'
#' Given a metagenome represented by a gbt object, and several bins derived
#' from that metagenome, this function allows quick plotting of the bins
#' simultaneously as overlays.
#' This is useful for checking the results of automatic binning software,
#' for example.
#' 
#' Note that plotting markers (the marker, ssu, and trna parameters) are
#' not allowed for this type of plot, unlike the usual plot.gbt function.
#' 
#' The gbtbin objects to be plotted must be specified as a list, e.g.
#' list(bin1,bin2,bin3). If you use the c() function there will be an error.
#' 
#' You can specify a custom color palette for the bins, or the default of
#' rainbow palette will be used.
#' To include a legend with which color corresponds to which name, you can
#' specify the bin names with the binNames parameter. if this is not given,
#' then the bins are simply numbered sequentially in the legend, according
#' to the order they were specified in the bins= parameter.
#' 
#' @param x Object of class gbt
#' @param bins List of gbtbin objects
#' @param binNames Vector of names for the bins (only used to draw legend, default none)
#' @param slice For plotting coverage data, which sample to use? 
#' @param cols Custom palette for bins overlay (default Rainbow)
#' @param legend Draw legend? (logical, default FALSE)
#' @param cutoff Minimum length to plot contigs (numeric, default 1000)
#' @param assemblyName Name of the metagenome, for plot title
#' @param log Which axes should be logarithmic scale
#' @param main Custom label for plot title
#' @param xlab Custom label for x-axis label
#' @param ylab Custom label for y-axis label
#' @param xlim Custom x-axis range
#' @param ylim Custom y-axis range
#'
#' @return New graphics object and plot, or error message
#' @export
#' @seealso \code{\link{plot.gbt}}
#'

multiBinPlot <- function (x, # Object of class gbt
                          bins, # List of gbtbin objects (must use list() instead of c() !!)
                          binNames="", # Vector of names for each bin
                          legend=FALSE, # Logical - plot legend identifying each bin?
                          slice=1, # which slice to use for plotting (default = 1)
                          assemblyName="", # Name of assembly, for plot header
                          cutoff=0, # Length cutoff for plotting
                          cols=NULL, # User-specified palette
                          log="default",  # Log scale for axis?
                          main="default",  # Custom title for plot
                          xlab="default",  # Custom x-axis label for plot
                          ylab="default",  # Custom y-axis label for plot
                          xlim=NULL, # Default xlim and ylim are automatic, but user can override
                          ylim=NULL
                          ) {
    # Check if user-supplied palette matches no. of bins, else use default rainbow palette
    if (length(cols) != length(bins)) {
        cat ("\nUsing default rainbow palette...\n")
        cols <- rainbow (length(bins))
    }
    # Check if bin names match number of bins, else ignore legend plotting
    if (legend==TRUE && binNames == "") {
        cat ("\nNo names for bins supplied, using sequential numbers...\n")
        # Sequentially name the bins "bin 1, bin 2, ..." if no custom names supplied
        binNames <- 1:length(bins)
        binNames <- sapply(binNames, function(x) paste("bin",x))
    }
    else if (legend==TRUE && length(binNames) != length(bins)) {
        cat ("\nNumber of bin names doesn't match number of bins supplied, ignoring legend...\n")
        legend <- FALSE
    }
    # Check that gbt object is really a gbt object, else abandon plotting
    if (class(x) != "gbt") {
        cat ("\nObject ")
        print (as.character (x))
        cat (" is not a gbt object!\n")
    }
    else {
        plot.gbt(x,
                 slice=slice, # Inherit user-specified options
                 assemblyName=assemblyName,
                 cutoff=cutoff,
                 log=log,main=main,xlab=xlab,ylab=ylab,
                 xlim=xlim,ylim=ylim,
                 col="grey",
                 marker=F,gc=F,ssu=F,trna=F # Turn off plot annotations
                 )
        for (i in 1:length(bins)) {
            if (class(bins[[i]]) != "gbtbin") {
                cat ("\nObject number ")
                print(as.character(i))
                cat (" is not a gbtbin object!\n")
            }
            else {
                points.gbtbin(x=bins[[i]],
                              col=cols[i],
                              slice=slice,
                              cutoff=cutoff)
            }
        }
        if (legend==TRUE) {
            legend(x="topright",
                   legend=binNames,
                   fill=cols
                   )
        }
    }
}