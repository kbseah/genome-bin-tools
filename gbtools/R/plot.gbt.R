#' Plot object of class gbt or gbtbin
#'
#' Plot GC-coverage or differential coverage plots from a gbt or gbtbin object
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
#' If you wish to plot by user-supplied custom values (added to gbt object by
#' the userAdd() function), indicate which ones to use for x- or y- axis in
#' the userAxis= parameter. For example, to plot user-custom data set 1 in X-
#' and user-custom data set 2 in Y-axis, specify userAxis=c(1,2). To plot vs.
#' GC or coverage values, supply "gc" or "cov" to userAxis, e.g. GC as X-axis
#' and user-value set 2 as Y-axis: userAxis=c("gc",2). If "cov" specified, it
#' will take the coverage set that is given by the slice= parameter (default
#' is the first set of coverage data).
#'
#' @param x Object of class gbt or gbtbin
#' @param slice For plotting coverage data, which sample to use? (see Details
#'        section below)
#' @param cutoff Minimum length to plot contigs (numeric, default 1000)
#' @param taxonLevel Taxonomic level for coloring the taxonomic markers, e.g.
#'        "Class" or "Phylum". (default "Class")
#' @param assemblyName Name of the metagenome, for plot title
#' @param marker Color plot by taxon markers? (logical, default TRUE)
#' @param marksource Specify which marker set to plot (default: first supplied)
#' @param markCutoff Length x coverage weight cutoff for colored markers
#'        (default: 0.99)
#' @param markCustomPalette (Optional) Assign colors to use for markers
#'        belonging to specific taxa. This must be a data frame with taxon names
#'        in column 1, and color names (English or hex codes) in column 2. Taxa
#'        that are not represented in the gbt object will be ignored. Taxa not
#'        in this table will be colored grey. If custom colors are specified,
#'        then the markCutoff parameter will be ignored. (default: 0)
#' @param highlightTaxon Color markers affiliated to specified taxon only
#'        (default: "")
#' @param gc Color plot by GC% instead of taxon markers? Only used for
#'        differential coverage plots, i.e. when two values are supplied to the
#'        slice parameter. (logical, default FALSE)
#' @param userAxis Use user-custom values for axis. See Details
#' @param ssu Draw markers for SSU genes? (logical, default FALSE)
#' @param trna Draw markers for tRNA genes? (logical, default FALSE)
#' @param consensus For contigs with more than one marker gene with conflict
#'        in taxonomy, take majority rule consensus? (logical, default TRUE)
#' @param legend Draw legend? (logical, default FALSE)
#' @param textlabel Label SSU markers with taxon? (logical, default FALSE)
#' @param col Color for points (default "grey")
#' @param log Which axes should be logarithmic scale
#' @param main Custom label for plot title
#' @param xlab Custom label for x-axis label
#' @param ylab Custom label for y-axis label
#' @param symbolScale Scaling function to use for plot characters (options:
#'        "area", "length")
#' @param symbolScaleParam Scaling constant for plot chars (default: 100)
#' @param ... See par for more plot options
#'
#' @return New graphics object and plot, or error message
#' @export
#' @seealso \code{\link{gbt}}
#'
plot.gbt <- function(x,  # Object of class gbt
                     slice=1,  # Which coverage values to plot?
                     cutoff=1000,  # Minimum contig length to plot
                     taxon="Class",  # Taxonomic level for coloring markers (keep for legacy compatibility)
                     taxonLevel="", # Alias for parameter "taxon"
                     highlightTaxon="", # Overlay colored markers for this specific taxon only
                     assemblyName="",  # Assembly name, for plot title only
                     marker=TRUE,  # Display marker color overlay
                     marksource="",  # Which marker source to plot; if empty - default first
                     markCutoff=0.99, # Weight cutoff for marker overlay
                     markCustomPalette=0, # Custom palette for marker taxons (Experimental)
                     gc=FALSE,  # Diffcov plot only: Color by GC instead of marker
                     ssu=FALSE,  # Overlay SSU markers?
                     trna=FALSE,  # Overlay tRNA markers?
                     consensus=TRUE,  # Conflicting marker taxon assgs: take majority consensus?
                     legend=FALSE,  # Add legend to color? (GC or marker taxon colors)
                     textlabel=FALSE,  # Add text labels to SSU markers?
                     userAxis=NA,  # Use userTab variables for axis
                     col="grey",  # Color for contig plots
                     log="default",  # Log scale for axis?
                                     # Default y for GC-cov plot
                                     # Default xy for differential coverage plot
                     main="default",  # Custom title for plot
                     xlab="default",  # Custom x-axis label for plot
                     ylab="default",  # Custom y-axis label for plot
                     symbolScale="area", # Type of scaling to use for plot symbols
                     symbolScaleParam=100, # Scaling parameter for the plot symbols
                     ...) {
## Plot method for gbt objects
    ## Process "taxonLevel" alias for "taxon" parameter
    if (taxonLevel != taxon) {
        if (taxonLevel == "" && taxon != "") {
            taxonLevel <- taxon
        } else {
            taxon <- taxonLevel
        }
    }
    ## Catch missing "slice" parameter
    if (is.na(slice[1]) || !is.numeric(slice)) {
        cat("gbtools ERROR: Please supply valid value for slice option\n")
    }
    ## Generate data.frame subset for user-custom axis values #################
    else if (!is.na(userAxis[1])) {
        if (length(userAxis)>0 ){
            ## Catch exceptions ###################################################
            if (length (userAxis)>2 ) {
                cat ("gbtools ERROR: Cannot supply more than 2 userAxis values\n")
            }
            else {
                ## Define X-axis ##################################################
                if (userAxis[1]=="gc") {
                    userX.x <- data.frame(ID=x$scaff$ID,
                                          Ref_GC=x$scaff$Ref_GC,
                                          Length=x$scaff$Length,
                                          Avg_fold=x$scaff$Avg_fold,
                                          xVals=x$scaff$Ref_GC)
                } else if (userAxis[1]=="cov") {
                    userX.x <- merge(data.frame(ID=x$scaff$ID,
                                                Ref_GC=x$scaff$Ref_GC,
                                                Length=x$scaff$Length,
                                                Avg_fold=x$scaff$Avg_fold),
                                     data.frame(ID=x$covs$ID,
                                                xVals=x$covs[[slice[1]+1]]),
                                     by="ID")
                } else if (any(userAxis[1]==x$userSource)) {
                    sourcenum <- which(x$userSource==userAxis[1])
                    userX.x <- merge(data.frame(ID=x$scaff$ID,
                                                Ref_GC=x$scaff$Ref_GC,
                                                Length=x$scaff$Length,
                                                Avg_fold=x$scaff$Avg_fold),
                                     data.frame(ID=x$userTab[[sourcenum]]$scaffold,
                                                xVals=x$userTab[[sourcenum]][2]),
                                     by="ID")
                } else {
                    cat ("gbtools ERROR: Invalid value for userAxis\n")
                }
                ## Define Y-axis ##################################################
                if (userAxis[2]=="gc") {
                    userX.y <- data.frame(ID=x$scaff$ID,
                                          yVals=x$scaff$Ref_GC)
                } else if (userAxis[2] =="cov") {
                    userX.y <- data.frame(ID=x$covs$ID,
                                          yVals=x$covs[[slice[1]+1]])
                } else if (any(userAxis[2]==x$userSource)) {
                    sourcenum <- which(x$userSource==userAxis[2])
                    userX.y <- data.frame(ID=x$userTab[[sourcenum]]$scaffold,
                                          yVals=x$userTab[[sourcenum]][2])
                } else {
                    cat ("gbtools ERROR: Invalid value for userAxis\n")
                }
                ## Make plot data.frame ###########################################
                X <- merge(userX.x,
                           userX.y,
                           by="ID")
                names(X) <- c("ID","Ref_GC","Length","Avg_fold","xVals","yVals")
                if (cutoff > 0) {
                    X <- subset(X,Length>=cutoff)
                }
                ## Set plot parameters ############################################
                if (main=="default") {  # Plot title
                    main <- paste("Plot for metagenome ",
                                  as.character(assemblyName))
                }
                if (xlab=="default") {  # X-axis label default
                    if (userAxis[1]=="gc") {
                        xlab <- "GC"
                    } else if (userAxis[1]=="cov") {
                        xlab <- paste("Coverage ",
                                      as.character(slice[1]))
                    } else {
                        xlab <- paste("User supplied variable ",
                                      as.character(userAxis[1]))
                    }
                }
                if (ylab=="default") {  # Y-axis label default
                    if (userAxis[2]=="gc") {
                        ylab <- "GC"
                    } else if (userAxis[2]=="cov") {
                        ylab <- paste("Coverage ",
                                      as.character(slice[1]))
                    } else {
                        ylab <- paste("User supplied variable ",
                                      as.character(userAxis[2]))
                    }
                }
                if (log=="default") {  # Default log scale only for cov parameter
                    if (userAxis[1]=="cov" &&
                        userAxis[2]!="cov") {
                        log <- "x"
                    } else if (userAxis[2]=="cov" &&
                               userAxis[1]!="cov") {
                        log <- "y"
                    } else {
                        log <- ""
                    }
                }
            }
        }
    }
    ## Generate data.frame subset for GC-Coverage plot ########################
    else if (!is.na(slice[1]) && length(slice)==1) {
        ## Make a new data.frame for plotting #################################
        X <- merge(data.frame(ID=x$scaff$ID,
                              Ref_GC=x$scaff$Ref_GC,
                              Length=x$scaff$Length,
                              Avg_fold=x$scaff$Avg_fold,
                              xVals=x$scaff$Ref_GC),
                   data.frame(ID=x$covs$ID,
                              yVals=x$covs[[slice[1]+1]]),
                   by="ID")
        names(X) <- c("ID","Ref_GC","Length","Avg_fold","xVals","yVals")
        ## Do the basic plot ##################################################
        if (cutoff > 0) {  # Minimum length cutoff for contigs to be plotted
            X <- subset(X,Length >= cutoff)
        }
        ## Set plot parameters ################################################
        if (main=="default") {  # Plot title
            main <- paste("Coverage-GC plot for metagenome ",
                          as.character(assemblyName))
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
    }
    ## Differential coverage plot #############################################
    else if (!is.na(slice[1]) && length(slice)==2) {
        ## Make a new data.frame for plotting #################################
        X <- merge (data.frame (ID=x$scaff$ID,
                                Ref_GC=x$scaff$Ref_GC,
                                Length=x$scaff$Length,
                                Avg_fold=x$scaff$Avg_fold),
                    data.frame (ID=x$covs$ID,
                                xVals=x$covs[[slice[1]+1]],
                                yVals=x$covs[[slice[2]+1]]),
                    by="ID")
        names(X) <- c("ID","Ref_GC","Length","Avg_fold","xVals","yVals")
        ## Plot parameters ####################################################
        if (cutoff > 0 ) {
            X <- subset(X,Length>=cutoff)
        }
        if (main == "default") {  # Default plot title
            main <- paste("Differential coverage plot for metagenome ",
                          as.character(assemblyName))
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
    }
    ## Basic plot without GC coloring #########################################
    if (!gc) {
        ## Draw the basic plot ################################################
        plot(x=X$xVals,
             y=X$yVals,
             pch=20,cex=cexScaling(X$Length, type=symbolScale, const=symbolScaleParam),
             col=col,
             log=log,
             main=main,
             xlab=xlab,
             ylab=ylab,
             ...)
        ## Add marker taxonomy overlay ########################################
        if (marker && !is.na(x$markTab)) {
            ## If marksource not specified, use first marker set by default
            if (marksource == "") {
                marksource <- x$markSource[1]  # Default: Display markers from first source
            } else if (is.character (marksource)
                       && length(marksource)==1) {
                ## Catch cases where marksource doesn't match entries in x$marksource
                if (!any(x$markSource==marksource)) {
                    cat ("gbtools WARNING: marksource doesn't match any entries,
                         defaulting to the first...\n")
                    marksource <- x$markSource[1]
                }
            } else if (is.character(marksource)
                       && length(marksource) > 1) {
                cat ("gbtools WARNING: Only one marker source can be plotted
                     as overlay (marksource parameter). Defaulting to first supplied\n")
                marksource <- marksource[1]
            } else {
                cat ("gbtools ERROR: Please check marksource argument\n")
            }
            ## Overlay markers ###############################################
            markTabTemp <- subset(x$markTab,source==marksource)
            mark.stats <- generatePlotColors2(X,markTabTemp,taxon,consensus,markCutoff,markCustomPalette)
            if (highlightTaxon != "") { ## Highlight specific taxon ##############
                # Subset mark.stats to keep only the specified taxon
                mark.stats <- mark.stats[which(mark.stats$taxon==highlightTaxon),]
                if (dim(mark.stats)[1] == 0) { # If mark.stats contains no rows, warning
                    cat ("gbtools WARNING: Taxon specified by highlightTaxon matches none!\n")
                }
            }
            points(x=mark.stats$xVals,
                   y=mark.stats$yVals,
                   pch=20,cex=cexScaling(mark.stats$Length,
                                         type=symbolScale,
                                         const=symbolScaleParam
                                         ),
                   col=as.character(mark.stats$colors))
            if (legend) {
                # Get those taxa which are not below cutoff
                colorframe <- unique(data.frame(mark.stats$taxon,
                                                mark.stats$colors,
                                                stringsAsFactors=FALSE # Else get error with newrow
                                                )
                                     )
                names(colorframe) <- c("taxon","colors")
                # Check if any taxa are below cutoff
                if (dim(subset(colorframe,colors=="grey50"))[1] == 0) {
                    nogray <- 1
                } else {
                    nogray <- 0
                }
                colorframe <- subset(colorframe,colors!="grey50")
                if (highlightTaxon == "" && nogray == 0) { # Don't add "below cutoff" to legend if only highlighting one taxon
                    newrow <- c("(below cutoff)","grey50")
                    colorframe <-rbind (colorframe,newrow)
                }
                legend("topright",
                       legend=colorframe$taxon,
                       cex=0.6,
                       fill=as.character(colorframe$colors))
            }
        }
    }
    ## Add GC coloring overlay ################################################
    else if (gc) {
        ## Catch exception ####################################################
        if (marker && gc) {
            cat ("gbtools ERROR: Cannot color by both marker and GC values \n")
        } else {
            ## Define GC palette colors (from Albertsen scripts) ##############
            gbr <- colorRampPalette(c("green","blue","orange","red"))
            palette (adjustcolor(gbr(70)))
            ## Draw plot colored by GC ########################################
            plot(x=X$xVals,
                 y=X$yVals,
                 pch=20,
                 cex=cexScaling(X$Length, type=symbolScale, const=symbolScaleParam),
                 col=X$Ref_GC*100,  # This colors points by GC values
                 main=main,
                 log=log,
                 xlab=xlab,
                 ylab=ylab,
                 ...)
            if (legend) {
                legendcolors <- c("20","30","40","50","60","70","80")
                legend ("topright",
                        legend=as.character(legendcolors),
                        fill=as.numeric(legendcolors))
            }
        }
    }
    ## Add SSU marker overlay #################################################
    if (ssu && !is.na(x$ssuTab)) {
        ssu.stats <- mergeScaffMarker(X,
                                      x$ssuTab,
                                      taxon,
                                      consensus=FALSE)
        points(x=ssu.stats$xVals,
               y=ssu.stats$yVals,
               pch=10,
               cex=2,
               col="black")
        ## Text label for SSU markers #########################################
        if (textlabel==TRUE) {
            text(x=ssu.stats$xVals,
                 y=ssu.stats$yVals,
                 as.character(ssu.stats$taxon), # Taxonomic level to get label
                 pos=3,
                 offset=0.2,
                 font=2)
        }
    }
    ## Add tRNA marker overlay ################################################
    if (trna && !is.na(x$trnaTab)) {
        trna.stats <- merge(X,
                            x$trnaTab,
                            by.x="ID",
                            by.y="scaffold")
        points(x=trna.stats$xVals,
               y=trna.stats$yVals,
               pch=4,
               cex=1,
               col="black")
    }
}
