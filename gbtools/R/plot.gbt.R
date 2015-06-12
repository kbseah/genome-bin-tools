#' Plot object of class gbt
#'
#' Plot GC-coverage or differential coverage plots from a gbt object
#'
#' The plot method for gbt objects can produce both GC-coverage and
#' differential coverage plots. A gbt object contains scaffold data and
#' annotations for a given metagenome, along with the coverage of each
#' scaffold for one or more samples. To plot GC-coverage plots, the sample
#' to use for the coverage data is specified by giving the number of the
#' sample to the slice parameter. To plot differential-coverage plots, the
#' numbers of the two samples for comparison are given to slice as a vector
#' e.g. c(1,2). The first sample is plotted as the x-axis and the second as
#' the y-axis. Supplying a vector with more than 2 elements, or a non-
#' numeric value will return an error.
#' Genome bins can be interactively chosen from a gbt plot with the
#' \code{\link{choosebin}} function; the same slice argument must be passed
#' to the choosebin function as the plot function, otherwise the results will
#' be meaningless!
#' 
#' @param x Object of class gbt
#' @param slice For plotting coverage data, which sample to use? (see
#'               Details section below)
#' @param cutoff Minimum length to plot contigs (numeric, default 1000)
#' @param taxon Taxonomic level for coloring the taxonomic markers, e.g.
#'               "Class" or "Phylum". (default "Class")
#' @param assemblyName Name of the metagenome, for plot title
#' @param marker Color plot by taxon markers? (logical, default TRUE)
#' @param gc Color plot by GC% instead of taxon markers? Only used for
#'            differential coverage plots, i.e. when two values are supplied
#'            to the slice parameter. (logical, default FALSE)
#' @param ssu Draw markers for SSU genes? (logical, default FALSE)
#' @param trna Draw markers for tRNA genes? (logical, default FALSE)
#' @param consensus For contigs with more than one marker gene with conflict
#'                   in taxonomy, take majority rule consensus? (logical,
#'                   default TRUE)
#' @param legend Draw legend? (logical, default FALSE)
#' @param textlabel Label SSU markers with taxon? (logical, default FALSE)
#' @param col Color for points (default "grey")
#' @param log Which axes should be logarithmic scale
#' @param main Custom label for plot title
#' @param xlab Custom label for x-axis label
#' @param ylab Custom label for y-axis label
#' @param ... See par for more plot options
#'
#' @return New graphics object and plot, or error message
#' @export
#' @seealso \code{\link{gbt}}
#'
plot.gbt <- function(x,  # Object of class gbt
                     slice=1,  # Which coverage values to plot?
                     cutoff=1000,  # Minimum contig length to plot
                     taxon="Class",  # Taxonomic level for coloring markers
                     assemblyName="",  # Assembly name, for plot title only
                     marker=TRUE,  # Display marker color overlay
                     marksource="",  # Which marker source to plot; if empty - default first
                     gc=FALSE,  # Diffcov plot only: Color by GC instead of marker
                     ssu=FALSE,  # Overlay SSU markers?
                     trna=FALSE,  # Overlay tRNA markers?
                     consensus=TRUE,  # Conflicting marker taxon assgs: take majority consensus?
                     legend=FALSE,  # Add legend to color? (GC or marker taxon colors)
                     textlabel=FALSE,  # Add text labels to SSU markers?
                     col="grey",  # Color for contig plots
                     log="default",  # Log scale for axis?
                                     # Default y for GC-cov plot
                                     # Default xy for differential coverage plot
                     main="default",  # Custom title for plot
                     xlab="default",  # Custom x-axis label for plot
                     ylab="default",  # Custom y-axis label for plot
                     ...) {
## Plot method for gbt objects
    if (is.na(slice[1]) || !is.numeric(slice)) {
        cat("gbtools ERROR: Please supply valid value for slice option\n")
    } else if (!is.na(slice[1]) && length(slice)==1) {
    ## GC-Coverage plot ##############################################
        ## Make a new data.frame for plotting ################################
        X <- merge( data.frame(ID=x$scaff$ID,Ref_GC=x$scaff$Ref_GC,Length=x$scaff$Length),
                   data.frame(ID=x$covs$ID,Avg_fold=x$covs[slice[1]+1]),
                   by="ID")
        names(X) <- c("ID","Ref_GC","Length","Avg_fold")
        
        ## Set plot parameters #################################################
        if (cutoff > 0) {  # Minimum length cutoff for contigs to be plotted
            X <- subset(X,Length >= cutoff)
        }
        if (main=="default") {  # Plot title
            main <- paste("Coverage-GC plot for metagenome ",as.character(assemblyName))
        }
        if (xlab=="default") {  # X-axis label default
            xlab <- "GC"
        }
        if (ylab=="default") {  # Y-axis label default
            ylab <- paste("Coverage ",as.character(slice[1]))
        }
        if (log=="default") {  # Default y-axis on logarithmic scale
            log <- "y"
        }
        
        ## Do the basic plot ###############################################
        plot(x=X$Ref_GC,
             y=X$Avg_fold,
             pch=20,cex=sqrt(X$Length)/100,
             col=col,log=log,
             main=main, xlab=xlab,ylab=ylab, ...)
        
        ## Add marker taxonomy overlay ########################################
        if (marker && !is.na(x$markTab)) {
            if (marksource == "") {
                marksource <- x$markSource[1]     # In default, display markers from the first source
            } else if (is.character(marksource) && length(marksource)==1) {
                # catch cases where marksource is not matching entries in x$markSource
                if (!any(x$markSource==marksource) ) {
                    cat ("gbtools WARNING: marksource doesn't match any entries.
                         Defaulting to first...\n")
                    marksource <- x$markSource[1]
                }
            } else if (is.character(marksource) && length(marksource) > 1) {
                cat ("gbtools WARNING: Only one marker source can be plotted
                     (marksource parameter). Defaulting to first supplied... \n")
                marksource <- marksource[1]
            } else {
                cat ("gbtools ERROR: Please check marksource argument\n")
            }
            markTabTemp <- subset(x$markTab,source==marksource)
            mark.stats <- generatePlotColors(X,markTabTemp,taxon,consensus)
            points(x=mark.stats$Ref_GC,
                   y=mark.stats$Avg_fold,
                   pch=20,cex=sqrt(mark.stats$Length)/100,
                   col=as.character(mark.stats$colors))
            if (legend) {
                colorframe <- generateLegendColors(X,markTabTemp,taxon,consensus)
                new.colorframe <- subset(colorframe,colors!="grey50")
                newrow <- c("singletons","grey50")
                new.colorframe <-rbind (new.colorframe,newrow)
                legend("topright",
                       legend=new.colorframe$taxon,
                       cex=0.6,
                       fill=as.character(new.colorframe$colors))
            }
        }
        
        ## Add SSU marker overlay ################################################################
        if (ssu && !is.na(x$ssuTab)) {
            ssu.stats <- mergeScaffMarker(X,x$ssuTab, taxon, consensus=FALSE)
            points(ssu.stats$Ref_GC, ssu.stats$Avg_fold,pch=10,cex=2,col="black")
            if (textlabel==TRUE) {
                text(ssu.stats$Ref_GC,ssu.stats$Avg_fold,as.character(ssu.stats$taxon),
                     pos=3,offset=0.2,font=2)
            }
        }
        
        ## Add tRNA marker overlay ###################################################################
        if (trna && !is.na(x$trnaTab)) {
            trna.stats <- merge(X, x$trnaTab, by.x="ID",by.y="scaffold")
            points(trna.stats$Ref_GC,trna.stats$Avg_fold,
                   pch=4,cex=1,col="black")
        }
    } else if (!is.na(slice[1]) && length(slice)==2) {
    ## Differential coverage plot ############################################################
        ## Make a new data.frame for plotting ##########################################################
        X <- merge (data.frame (ID=x$scaff$ID,
                                Ref_GC=x$scaff$Ref_GC,
                                Length=x$scaff$Length),
                    data.frame (ID=x$covs$ID,
                                Avg_fold_1=x$covs[slice[1]+1],
                                Avg_fold_2=x$covs[slice[2]+1]),
                    by="ID")
        names(X) <- c("ID","Ref_GC","Length","Avg_fold_1","Avg_fold_2")
        ## Plot parameters #########################################################
        if (cutoff > 0 ) {
            X <- subset(X,Length>=cutoff)
        }
        if (main == "default") {  # Default plot title
            main <- paste("Differential coverage plot for metagenome ",as.character(assemblyName))
        }
        if (xlab=="default") {  # Default X-axis label
            xlab <- paste("Coverage ",slice[1])
        }
        if (ylab=="default") {  # Default Y axis label
            ylab <- paste("Coverage ",slice[2])
        }
        if (log=="default") {  # Default both x- and y-axis on logarithmic scale
            log <- "xy"
        }
        ## Define GC palette colors (from Albertsen scripts) #######################
        gbr <- colorRampPalette(c("green","blue","orange","red"))
        ## Catch invalid coloring option combination ###############################
        if (gc && marker & !is.na(x$markTab)) {  
            cat("plase choose to plot only with GC or marker coloring, but not both!\n")
        }
        ## Color plot by GC% values ################################################
        else if (gc && !marker) {  
            palette (adjustcolor(gbr(70)))
            plot(x=X$Avg_fold_1,
                 y=X$Avg_fold_2,
                 pch=20,cex=sqrt(X$Length)/100,col=X$Ref_GC*100,
                 main=main,log=log,xlab=xlab,ylab=ylab,...)
            if (legend) {
                legendcolors <- c("20","30","40","50","60","70","80")
                legend ("topright",
                        legend=as.character(legendcolors),
                        fill=as.numeric(legendcolors))
            }
        }
        ## Color plot by Marker taxonomy ############################################
        else if (!gc && !is.na(x$markTab) && marker) {  
            ## Identify which marker source to plot
            if (marksource == "") {
                marksource <- x$markSource[1]     # In default, display markers from the first source
            } else if (is.character(marksource) && length(marksource)==1) {
                # catch cases where marksource is not matching entries in x$markSource
                if (!any(x$markSource==marksource) ) {
                    cat ("gbtools WARNING: marksource doesn't match any entries.
                         Defaulting to first...\n")
                    marksource <- x$markSource[1]
                }
            } else {
                cat ("gbtools ERROR: Only one marker source can be plotted
                     (marksource parameter)\n")
            }
            markTabTemp <- subset(x$markTab,source==marksource)
            ## Produce the base plot
            plot (x=X$Avg_fold_1,
                  y=X$Avg_fold_2,
                  pch=20,cex=sqrt(X$Length)/100,
                  main=main,col=col,log=log,ylab=ylab,xlab=xlab, ...)
            ## Generate plot colors for markers and overlay with markers
            mark.stats <- generatePlotColors(X,markTabTemp,taxon,consensus)
            points(x=mark.stats$Avg_fold_1,
                   y=mark.stats$Avg_fold_2,
                   pch=20,cex=sqrt(mark.stats$Length)/100,
                   col=as.character(mark.stats$colors))
            if (legend) {
                colorframe <- generateLegendColors(X,markTabTemp,taxon,consensus)
                new.colorframe <- subset(colorframe,colors!="grey50")
                newrow <- c("singletons","grey50")
                new.colorframe <- rbind (new.colorframe,newrow)
                legend("topright",
                       legend=new.colorframe$taxon,
                       cex=0.6,
                       fill=as.character(new.colorframe$colors))
            }
        }
        ## Do not color plot ###########################################################################
        else if (!gc && !marker) {  
            plot (X$Avg_fold_1,X$Avg_fold_2,pch=20,cex=sqrt(X$Length)/100,
                  main=main,col=col,log=log,xlab=xlab,ylab=ylab, ...)
        }
        if (!(gc&&marker)) {
            
            ## Add SSU markers to plot ############################################################################
            if (ssu && !is.na(x$ssuTab)) { 
                ssu.stats <-mergeScaffMarker (X,x$ssuTab,taxon,consensus=FALSE)
                points(x=ssu.stats$Avg_fold_1,
                       y=ssu.stats$Avg_fold_2,
                       pch=10,cex=2,col="black")
                if (textlabel==TRUE) {  # Add text labels to SSU markers
                    text(x=ssu.stats$Avg_fold_1,
                         y=ssu.stats$Avg_fold_2,
                         as.character(ssu.stats$taxon),pos=3,offset=0.2,font=2)
                }
            }
            ## Add tRNA markers to plot ##########################################################################
            if (trna && !is.na(x$trnaTab)) {  
                trna.stats <- merge(X,x$trnaTab,by.x="ID",by.y="scaffold")
                points(x=trna.stats$Avg_fold_1,
                       y=trna.stats$Avg_fold_2,
                       pch=4,cex=1,col="black")
            }
        }
    }
}