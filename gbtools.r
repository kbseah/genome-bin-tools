add <- function(x1, x2) UseMethod("add")

add.gbtbin <- function(x1,x2) {
## Merge two bins; i.e. take their union
    result <- setOperation(x1=x1,x2=x2,shortlist="all")
    result$call[[length(result$call)+1]] <- match.call()  # Record function call 
    return(result)
}

choosebin <- function(x, ... ) UseMethod ("choosebin")  # Defines generic for choosebin function

choosebin.gbt <- function(x,  # Object of class gbt
                          slice,  # Which slices used for the plot from which points to be chosen?
                          taxon="Class",  # Deprecated - user don't change this
                          num.points=6,  # Number of points to pick in polygon
                          draw.polygon=TRUE,  # Add polygon overlay to plot?
                          save=FALSE,  # Save list of contigs in bin to external file?
                          file="interactive_bin.list"  # Name of file to save list of contigs in bin
                          ) {
    require(sp)
## Wrapper for picking bin interactively from GC-cov or diff-cov plot
    if (!is.numeric(slice) || length(slice) > 2) {
        cat ("gbtools ERROR: Please specify the library(-ies) used to make the plot in focus\n")
    } else {
        bin <- pickBinPoints(num.points=num.points,
                             draw.polygon=draw.polygon)
        if (length(slice)==1) {  # Pick bin from GC-coverage plot
            X <- merge(data.frame(ID=x$scaff$ID,
                                  Ref_GC=x$scaff$Ref_GC),
                       data.frame(ID=x$covs$ID,
                                  Avg_fold=x$covs[slice[1]+1]),
                       by="ID")
            names(X) <- c("ID","Ref_GC","Avg_fold")
            inpolygon <- sp::point.in.polygon(X$Ref_GC,
                                              X$Avg_fold,
                                              bin$x,
                                              bin$y)
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
                                              bin$x,
                                              bin$y)
        }
        X.subset <- X[which(inpolygon==1),]
        X.shortlist <- as.character(X.subset$ID)
        result <- gbtbin(shortlist=X.shortlist,
                         x=x,
                         slice=slice,
                         taxon=taxon,
                         points=bin,
                         save=save,
                         file=file)
        result$call[[length(result$call)+1]] <- match.call()  # Record choosebin() function call
        return(result)
    }
}

countSingleFromTable <- function(x) {
    x.tab <- table(x)
    uniq <- length(which(x.tab==1))
    return(uniq)
}


fastgFishing <- function(x, bin, fastg.file, ... ) UseMethod ("fastgFishing") 

fastgFishing.gbtbin <- function(x,  # Object of class gbt (parent object of the gbtbin)
                                bin,  # Object of class gbtbin defined from x
                                fastg.file,  # Fastg file for assembly of x
                                taxon="Class",  # Deprecated - user pls ignore
                                save=FALSE,  # Save list of contigs to external file?
                                file="fished_bin.list"  # File to save contig list
                                ) {
    command <- "perl"
## REPLACE THIS PATH WITH YOUR OWN PATH !! #########################################
    script.path <- "/home/kbseah/tools/my_scripts/genome-bin-tools/fastg_parser.pl" 
####################################################################################
    command.params <- paste(script.path,"-i",fastg.file,"-o /tmp/tmp.fishing_output -b - -r")
        # By default throws away fastg_parser.pl output to /tmp/
    fished.contigs.list <- system2(command,
                                   command.params,
                                   input=as.character(bin$scaff$ID),
                                   stderr=NULL,
                                   stdout=TRUE)
    newbin <- gbtbin(shortlist=fished.contigs.list,
                     x=x,
                     slice=NA,
                     taxon=taxon,
                     save=save,
                     file=file)
    newbin$call[[length(newbin$call)+1]] <- match.call() # Record function call
    return(newbin)
}
gbtbin <- function(shortlist,x,slice,taxon,points,save,file) UseMethod("gbtbin")
gbtbin.default <- function(shortlist,  # Character vector, contigs to extract from gbt object
                           x,  # Object of class gbt
                           slice,  # Which slice used by choosebin()?
                           taxon,  # Deprecated - user don't change
                           points=NA,  # Number of points in polygon for choosebin()
                           save=FALSE,  # Save contig list to external file?
                           file="interactive_bin.list"  # File name to save contig list
                           ) {
    scaff.subset <- subset(x$scaff, ID%in%shortlist)
    covs.subset <- subset(x$covs, ID%in%shortlist)
    markTab.subset <- NA
    marksource <- x$markSource
    ssuTab.subset <- NA
    trnaTab.subset <- NA
    thecall <- x$call
    userTab <- list()
    userSource <- x$userSource
    ## Summary statistics initialize ##########################################
    bin.nummarkers <- NA
    bin.uniqmarkers <- NA
    bin.numtRNAs <- NA
    bin.uniqtRNAs <- NA
    bin.numSSUs <- NA
    bin.singlemarkers <- NA
    marker.tab <- list()
    tRNAs.tab <- NA
    ## Take subset of markTab #################################################
    if (is.data.frame(x$markTab)) {
        markTab.subset <- subset(x$markTab,scaffold%in% shortlist)
        for (j in 1:length(marksource)) {
            markTab.subsubset <- subset(markTab.subset,source==marksource[j])
            bin.nummarkers[j] <- dim(markTab.subsubset)[1]
            marker.tab[[j]] <- table(markTab.subsubset$gene)
            bin.uniqmarkers[j] <- length(which(marker.tab[[j]] > 0))
            bin.singlemarkers[j] <- length(which(marker.tab[[j]] ==1))
        }
    }
    ## Take subset of ssuTab ##################################################
    if (is.data.frame(x$ssuTab)) {
        ssuTab.subset <- subset(x$ssuTab,scaffold%in%shortlist)
        bin.numSSUs <- dim(ssuTab.subset)[1]
    }
    ## Take subset of trnaTab #################################################
    if (is.data.frame(x$trnaTab)) {
        trnaTab.subset <- subset(x$trnaTab,scaffold%in% shortlist)
        bin.numtRNAs <- dim(trnaTab.subset)[1]
        tRNAs.tab <- table (trnaTab.subset$tRNA_type)
        bin.uniqtRNAs <- length(which(tRNAs.tab > 0))
    }
    ## Take subset of userTab #################################################
    if (length(x$userTab) > 0) {
        for (i in 1:length(x$userTab)) {
            userTabSubset <- subset(x$userTab[[i]],scaffold %in% shortlist)
            userTab[[i]] <- userTabSubset
        }
    }
    ## Summary statistics #####################################################
    bin.length <- sum(scaff.subset$Length)
    bin.numscaffolds <- dim(scaff.subset)[1]
    bin.summary <- list(Total_length=bin.length,
                        Num_scaffolds=bin.numscaffolds,
                        Scaff_length_max=max(scaff.subset$Length),
                        Scaff_length_min=min(scaff.subset$Length),
                        Scaff_length_median=median(scaff.subset$Length),
                        Scaff_length_N50=getN50(scaff.subset$Length),
                        Marker_sources=marksource,
                        Num_markers=bin.nummarkers,
                        Num_unique_markers=bin.uniqmarkers,
                        Num_singlecopy_markers=bin.singlemarkers,
                        Num_SSUs=bin.numSSUs,
                        Num_tRNAs=bin.numtRNAs,
                        Num_tRNAs_types=bin.uniqtRNAs)
    
    ## Write to file, if save option is used ##################################
    if (save) {
        write(as.vector(scaff.subset$ID),file=file)
    }
    
    ## Package and return result ##############################################
    result <- list(scaff=scaff.subset,
                   covs=covs.subset,
                   markTab=markTab.subset,
                   markSource=marksource,
                   ssuTab=ssuTab.subset,
                   trnaTab=trnaTab.subset,
                   userTab=userTab,
                   userSource=userSource,
                   summary=bin.summary,
                   marker.table=marker.tab,
                   tRNA.table=tRNAs.tab,
                   points=points,
                   slice=slice,
                   call=thecall)
    class(result) <- "gbtbin"
    return(result)
}

gbt <- function(covstats,mark,marksource,ssu,trna) UseMethod ("gbt")

gbt.default <- function (covstats,  # Vector of filenames for coverage tables 
                         mark=NA,   # Vector of filenames for marker gene taxonomy tables
                         marksource=NA,  # Vector of source names for each marker gene table
                         ssu=NA,    # Filename for SSU annotation table
                         trna=NA    # Filename for tRNA annotation table
                         ) {
## Create new gbt objects
    if ( class(covstats)!="character" || length(covstats)==0 ) {  # Check that covstats argument is character class
        cat ("gbtools ERROR: covstats argument must be a list of file names!\n")
    }
    else {
        ## Read Coverage tables ##############################################
        if (length(covstats)==1) {
            scaff <- read.table(file=as.character(covstats),
                                sep="\t",
                                header=T)
            covs <- data.frame(ID=scaff$ID,
                               scaff$Avg_fold)
        }
        else {
            scaff <- read.table(file=as.character(covstats[1]),
                                sep="\t",
                                header=T)  # Contains all other data associated per contig
            covs <- data.frame(ID=scaff$ID,
                               scaff$Avg_fold)  # Contains scaffold ID and coverage data
            for (i in 2:length(covstats)) {  # Read the other covstats files
                scafftemp <- read.table(file=as.character(covstats[i]),
                                        sep="\t",
                                        header=T)
                covs <- merge(covs,
                              data.frame(ID=scafftemp$ID,
                                         scafftemp$Avg_fold),
                              by="ID")
            }
        }
        
        ## Read taxonomic marker table ##########################################
        if ( !is.na(marksource[1]) ) { # Read list of marker sources
            markTab <- NA
            if ( !is.na(mark[1]) ) {  # Read marker table
                for (i in 1:length(marksource)) {
                    markTabTemp <- read.table(file=as.character(mark[i]),
                                              sep="\t",  # Tab separated input only
                                              fill=TRUE,  # Fill ragged lines 
                                              header=T)
                    namesvec <- names(markTabTemp)
                    if (length(which(!markTabTemp$scaffold %in% scaff$ID)) > 0) {
                        # Catch marker tables where scaffold IDs dont match
                        # what's in the scaff table
                        cat ("gbtools WARNING: Scaffold IDs in marker table ")
                        print (as.character(mark[i]))
                        cat (" doesn't match contig coverage tables \n")
                    }
                    numMarksTemp <- dim(markTabTemp)[1]
                    sourcevector <- rep(as.character(marksource[i]),numMarksTemp)
                    markTabTemp <- cbind(markTabTemp,sourcevector)
                    names(markTabTemp) <- c(namesvec,"source")
                    if ( length(markTab)==1 && is.na (markTab)) {
                        markTab <- markTabTemp
                    } else {
                        markTab <- rbind(markTab,markTabTemp)
                    }
                }
            } else {
                cat ("gbtools WARNING: Marker tables not supplied. marksource parameter ignored \n")
                marksource <- NA
                numMarks <- NA
                markTab <- NA
            }
        } else {
            cat ("gbtools WARNING: marksource not supplied. Any marker tables supplied ignored \n")
            marksource <- NA
            numMarks <- NA
            markTab <- NA
        }
        
        ## Read SSU marker table ##############################################
        if ( !is.na(ssu[1]) ) {  
            ssuTab <- read.table(file=as.character(ssu),sep="\t",header=T)
            numSsu <- dim(ssuTab)[1]
        } else {
            ssu <- NA
            ssuTab <- NA
            numSsu <- NA
        }
        
        ## Read tRNA marker table #############################################
        if ( !is.na(trna[1]) ) {  # Read tRNA marker table
            trnaTab <- read.table(file=as.character(trna),sep="\t",skip=3,header=F)
            names(trnaTab) <- c("scaffold",
                                "tRNA_no",
                                "tRNA_begin",
                                "tRNA_end",
                                "tRNA_type",
                                "Anticodon",
                                "Intron_begin",
                                "Intron_end",
                                "Cove_score")
            numTrna <- dim(trnaTab)[1]
        } else {
            trna <- NA
            trnaTab <- NA
            numTrna <- NA
        }
        
        ## Generate summary statistics #######################################
        if (!is.na(markTab)[1]) {
            Num_markers <- table(markTab$source)
        } else { Num_markers <- NA }
        summarystats <- list(Total_length=sum(scaff$Length),
                             Num_scaffolds=length(scaff$ID),
                             Scaff_length_max=max(scaff$Length),
                             Scaff_length_min=min(scaff$Length),
                             Scaff_length_median=median(scaff$Length),
                             Scaff_length_N50=getN50(scaff$Length),
                             Num_markers=Num_markers,
                             Num_SSU=numSsu,
                             Num_tRNAs=numTrna)
        
        ## Package and return result #########################################
        userTab <- list()  # Create userTab as an empty list to hold user-custom data
        userSource <- ""  # Create userSource vector to hold names of user-custom data
        result <- list(scaff=scaff,
                       covs=covs,
                       markTab=markTab,
                       markSource=marksource,
                       ssuTab=ssuTab,
                       trnaTab=trnaTab,
                       userTab=userTab,
                       userSource=userSource,
                       summary=summarystats)
        result$call <- list()
        result$call[[1]] <- match.call()  # Record function call that produces this gbt object
        class(result) <- "gbt"
        result
    }    
}

generateLegendColors <- function(scaffold.stats,  # Same params as generatePlotColors()
                                 marker.list,
                                 taxon,
                                 consensus) {
## Generates colors for plot legends when coloring by markers
    ## Some table merging to have points to plot for the markers ###################
    marker.stats <- mergeScaffMarker(scaffold.stats,
                                     marker.list,
                                     taxon,
                                     consensus) 
    marker.list[,"taxon"] <- marker.list[,which(names(marker.list)==taxon)]
    ## Count how many taxa are only supported by one marker gene ###################
    singleton.taxa <- names(table(marker.list$taxon)[which(table(marker.list$taxon)==1)])
    ## Which taxon has the most marker genes? ######################################
    top.taxon <- names(table(marker.list$taxon)[which.max(table(marker.list$taxon))]) 
    taxnames <- names(table(marker.stats$taxon))  # Names of taxa
    ## Vector to hold color names ##################################################
    taxcolors <- rep("",length(names(table(marker.stats$taxon))))
    taxcolors[which(names(table(marker.stats$taxon))  # Give singleton taxa color "grey50" 
                    %in% singleton.taxa)] <- "grey50"
    numsingletons <- length(taxcolors[which(names(table(marker.stats$taxon))
                                            %in% singleton.taxa)]) # Count singleton taxa
    taxcolors[which(names(table(marker.stats$taxon))
                    ==top.taxon)] <- "red"  # Give taxon with most marker genes color "red"
    ## How many other colors do we need, given that all singletons have same color? #
    numcolors <- length(table(marker.stats$taxon)) - 1 - numsingletons
    ## Generate needed colors, from yellow to magenta, giving red a wide berth ######
    thecolors <- rainbow(numcolors,
                         start=1/6,
                         end=5/6)
    taxcolors[which(!(names(table(marker.stats$taxon)) %in% singleton.taxa)
                    & names(table(marker.stats$taxon))!=top.taxon)] <- thecolors
    ## Data frame containing which colors correspond to which taxa ##################
    colorframe <- data.frame(taxon=taxnames,
                             colors=taxcolors) 
    ## Merge this by taxon into the marker.stats table for plotting (this works #####
    ## even when consensus option is called) ########################################
    marker.stats <- merge(marker.stats,
                          colorframe,
                          by="taxon")
    return(colorframe)
}

generatePlotColors <- function(scaffold.stats,  # scaff table from gbt object
                               marker.list,  # markTab table from gbt object
                               taxon,  # Taxonomic level to do the coloring
                               consensus  # Logical-if taxon assgs conflict, take consens?
                               ) {           # This took a very long time to get it right
## Generates colors for marker gene phylotypes in plot
    ## Merge tables to have points to plot for the markers #########################
    marker.stats <- mergeScaffMarker(scaffold.stats,
                                     marker.list,
                                     taxon,
                                     consensus) 
    marker.list[,"taxon"] <- marker.list[,which(names(marker.list)==taxon)]
    ## Count taxa supported only by one marker gene ################################
    singleton.taxa <- names(table(marker.list$taxon)[which(table(marker.list$taxon)==1)])
    ## Which taxon has the most marker genes? ######################################
    top.taxon <- names(table(marker.list$taxon)[which.max(table(marker.list$taxon))])
        # Use which.max() because it breaks ties. Otherwise all genomes with same
        # number of marker genes will have same color!
        # Important: Identification of singleton taxa uses the original marker.list
        #            because after "consensus", each scaffold has only one taxon
        #            assignment and scaffolds with >1 marker will be undercounted
    ## For plot colors - identify singleton taxa and the taxon with the highest ####
    ## marker counts, and assign them special colors ###############################
    taxnames <- names(table(marker.stats$taxon))  # Names of taxa
    taxcolors <- rep("",length(names(table(marker.stats$taxon))))  # Vector to hold color names
    taxcolors[which(names(table(marker.stats$taxon))
                    %in% singleton.taxa)] <- "grey50"  # Give singletons color "grey50"
    numsingletons <- length(taxcolors[which(names(table(marker.stats$taxon))
                                            %in% singleton.taxa)]) # Count singleton taxa
    taxcolors[which(names(table(marker.stats$taxon))
                    ==top.taxon)] <- "red"  # Give taxon with most marker genes color "red"
    ## How many other colors do we need, given that all singletons have same color?#
    numcolors <- length(table(marker.stats$taxon)) - 1 - numsingletons
    ## Generate needed colors, from yellow to magenta, giving red a wide berth #####
    thecolors <- rainbow(numcolors,
                         start=1/6,
                         end=5/6) 
    taxcolors[which(!(names(table(marker.stats$taxon))%in% singleton.taxa)
                    & names(table(marker.stats$taxon))!=top.taxon)] <- thecolors
    ## Data frame containing which colors correspond to which taxa #################
    colorframe <- data.frame(taxon=taxnames,
                             colors=taxcolors)
    ## Merge this by taxon into the marker.stats table for plotting (this works ####
    ## even when consensus option is called) #######################################
    marker.stats <- merge(marker.stats,
                          colorframe,
                          by="taxon")
    return(marker.stats)
}

getN50 <- function(x) {
    # Adapted from R-bloggers:
    # http://www.r-bloggers.com/calculating-an-n50-from-velvet-output/
    x.sort <- sort(x,decreasing=TRUE)
    n50 <- x[cumsum(x.sort) >= sum(x.sort)/2][1]
    return(n50)
}

lej <- function(x1, x2) UseMethod ("lej")

lej.gbtbin <- function(x1,x2) {
## Take difference between two bins - non commutative! i.e. left exclusive join
    shortlist <- x1$scaff$ID[which(!x1$scaff$ID %in% x2$scaff$ID)]
    result <- setOperation.gbtbin(x1=x1,
                                  x2=x2,
                                  shortlist=shortlist)
    result$call[[length(result$call)+1]] <- match.call()  # Record function call 
    return(result)
}

mergeScaffMarker <- function(scaffold.stats,marker.list,taxon,consensus=TRUE) {
## Merge table of scaffold statistics (output from pileup.sh in BBMap package)
## and table of marker statistics parsed by parse_phylotype_result.pl
    require(plyr)
    marker.list[,"taxon"] <- marker.list[,which(names(marker.list)==taxon)]
    marker.stats <- merge(scaffold.stats,
                          marker.list,
                          by.x="ID",
                          by.y="scaffold")
    ## For scaffolds with multiple marker genes, take majority consensus of marker #
    ## taxon assignment ############################################################
    if (consensus) {    
        #scaffs.with.multi <- as.vector(names(table(marker.stats$ID)[which(table(marker.stats$ID)>1)]))
        #consensus.list <- ddply(marker.list,
        #                       .(scaffold),
        #                       function(x) levels(x$taxon)[which.max(tabulate(x$taxon))])
        consensus.list <- plyr::ddply(marker.list,
                                      .(scaffold),
                                      summarize,
                                      taxon=levels(taxon)[which.max(tabulate(taxon))])
        marker.stats <- merge(scaffold.stats,
                              consensus.list,
                              by.x="ID",
                              by.y="scaffold")
    }
    return(marker.stats)
}

pickBinPoints <- function(num.points=6,  # How many points in polygon?
                          draw.polygon=TRUE  # Overlay polygon on plot?
                          ) {
## Wrapper for locator() and polygon() to perform interactive binning on the current
## plot. Returns the polygon vertices which can be used in get.bin.stats()
    thepoints <- locator(num.points,pch=20,type="p")
    if (draw.polygon) { polygon(thepoints) }
    return(thepoints)
}

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
                     userAxis=NA,  # Use userTab variables for axis
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
                                          xVals=x$scaff$Ref_GC)
                } else if (userAxis[1]=="cov") {
                    userX.x <- merge(data.frame(ID=x$scaff$ID,
                                                Ref_GC=x$scaff$Ref_GC,
                                                Length=x$scaff$Length),
                                     data.frame(ID=x$covs$ID,
                                                xVals=x$covs[slice[1]+1]),
                                     by="ID")
                } else if (any(userAxis[1]==x$userSource)) {
                    sourcenum <- which(x$userSource==userAxis[1])
                    userX.x <- merge(data.frame(ID=x$scaff$ID,
                                                Ref_GC=x$scaff$Ref_GC,
                                                Length=x$scaff$Length),
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
                                          yVals=x$covs[slice[1]+1])
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
                names(X) <- c("ID","Ref_GC","Length","xVals","yVals")
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
                              xVals=x$scaff$Ref_GC),
                   data.frame(ID=x$covs$ID,
                              yVals=x$covs[slice[1]+1]),
                   by="ID")
        names(X) <- c("ID","Ref_GC","Length","xVals","yVals")
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
                                Length=x$scaff$Length),
                    data.frame (ID=x$covs$ID,
                                xVals=x$covs[slice[1]+1],
                                yVals=x$covs[slice[2]+1]),
                    by="ID")
        names(X) <- c("ID","Ref_GC","Length","xVals","yVals")
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
             pch=20,cex=sqrt(X$Length)/100,
             col=col,
             log=log,
             main=main,
             xlab=xlab,
             ylab=ylab,
             ...)
        ## Add marker taxonomy overlay ########################################
        if (marker && !is.na(x$markTab)) {
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
            markTabTemp <- subset(x$markTab,source==marksource)
            mark.stats <- generatePlotColors(X,markTabTemp,taxon,consensus)
            points(x=mark.stats$xVals,
                   y=mark.stats$yVals,
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
                 cex=sqrt(X$Length)/100,
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

plot.gbtbin <- plot.gbt # Inherit plot.gbt function

points.gbtbin <- function(x,  # Object of class gbtbin
                          col="black",  # Overlay plot points color
                          slice="default",  # Which slice to use for plotting?
                          cutoff=0,  # Min contig length to plot
                          pch=20, ...) {
    ## Defaults to the same slice used to choose the bin ###########################
    if (slice == "default") { 
        slice <- x$slice
    }
    ## Catch invalid slice values ##################################################
    if (is.na(slice) || !is.numeric(slice) || length(slice) > 2) { 
        cat ("gbtools ERROR: Please specify valid value for slice option for this bin\n")
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
               cex=sqrt(as.numeric(X$Length))/100,
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
               cex=sqrt(as.numeric(X$Length))/100,
               col=col, ...)
    }
    ## Throw error message for invalid slice options ################################
    else { cat ("gbtools ERROR: Please specify valid value for slice option for this bin\n")}
}

print.gbtbin <- function(x) {
    cat("*** Object of class gbtbin ***\n")
    cat ("\n*** Scaffolds ***\n")
    lengthdf <- data.frame(x$summary$Total_length,
                           x$summary$Num_scaffolds,
                           x$summary$Scaff_length_max,
                           x$summary$Scaff_length_min,
                           x$summary$Scaff_length_median,
                           x$summary$Scaff_length_N50)
    names(lengthdf) <- c("Total", "Scaffolds", "Max", "Min", "Median", "N50")
    print(lengthdf)
    cat("\n*** Markers ***\n")
    markerdf <- data.frame(x$summary$Marker_sources,
                           x$summary$Num_markers,
                           x$summary$Num_unique_markers,
                           x$summary$Num_singlecopy_markers)
    names(markerdf) <- c("Source","Total","Unique","Singlecopy")
    print(markerdf)
    cat("\n*** SSU markers ***\n")
    print(x$summary$Num_SSUs)
    cat("\n*** tRNA_markers ***\n")
    print(x$summary$Num_tRNAs)
    cat("\n*** User-supplied variables ***\n")
    print(x$userSource)
    cat("\n*** Function call history ***\n")
    print(x$call)
}

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

setOperation <- function(x1, x2, shortlist) UseMethod("setOperation")

setOperation.gbtbin <- function(x1,x2,shortlist) {
## Generic operation for merging and subsetting two gbtbin objects
    ## Merge the two scaff tables ##################################################
    scaff.add <- unique (rbind(x1$scaff,x2$scaff))
    ## Catch special value for shortlist ###########################################
    if (shortlist[1] == "all") {
        shortlist <- as.character(scaff.add$ID)
    }
    ## Subset scaff and coverage tables ############################################
    scaff.add <- subset (scaff.add, ID %in% shortlist)
    covs.add <- unique (rbind(x1$covs,x2$covs))
    covs.add <- subset(covs.add, ID %in% shortlist)
    ## Initialize marker, ssu, and tRNA tables #####################################
    markTab.add <- NA
    ssuTab.add <- NA
    trnaTab.add <- NA
    marksource <- x1$markSource  # markSource should be identical in x1 and x2 -- TODO: Implement a check?
    ## Initialize variables for summary stats ######################################
    bin.nummarkers <- NA    # In case the marker.list is not supplied
    bin.uniqmarkers <- NA
    bin.numtRNAs <- NA      # Likewise for number of tRNAs
    bin.uniqtRNAs <- NA
    bin.numSSUs <- NA
    bin.singlemarkers <- NA
    marker.tab <- list()  # Make empty list object (NA will throw error later)
    tRNAs.tab <- NA
    ## Combine marker tables #######################################################
    #if (!is.na(x1$markTab) || !is.na(x2$markTab)) {
    if (is.data.frame(x1$markTab) && is.data.frame(x2$markTab)) {
        if (dim(x1$markTab)[1] > 0 || dim(x2$markTab)[1] > 0) { 
            markTab.add <- unique (rbind(x1$markTab,x2$markTab))
            markTab.add <- subset(markTab.add, scaffold %in% shortlist)
            for (j in 1:length(marksource)) {
                markTab.subsubset <- subset(markTab.add,source==marksource[j])
                bin.nummarkers[j] <- dim(markTab.subsubset)[1]
                marker.tab[[j]] <- table(markTab.subsubset$gene)
                bin.uniqmarkers[j] <- length(which(marker.tab[[j]] > 0))
                bin.singlemarkers[j] <- length(which(marker.tab[[j]] ==1))
            }
        }
    }
    ## Combine SSU Tables ##########################################################
    #if (!is.na(x1$ssuTab) || !is.na(x2$ssuTab)) {
    if (is.data.frame(x1$ssuTab) && is.data.frame(x2$ssuTab)) {
        if ( dim(x1$ssuTab)[1]>0 || dim(x2$ssuTab)[1]>0 ) {
            ssuTab.add <- unique (rbind(x1$ssuTab,x2$ssuTab))
            ssuTab.add <-subset(ssuTab.add, scaffold %in% shortlist)
            bin.numSSUs <- dim(ssuTab.add)[1]
        }
    }
    ## Combine trna Tables #########################################################
    #if (!is.na(x1$trnaTab) || !is.na(x2$trnaTab)) {
    if (is.data.frame(x1$trnaTab) && is.data.frame(x2$trnaTab)) {
        if ( dim(x1$trnaTab)[1] > 0 || dim(x2$trnaTab)[1] > 0 ) {
            trnaTab.add <- unique (rbind(x1$trnaTab,x2$trnaTab))
            trnaTab.add <- subset(trnaTab.add, scaffold %in% shortlist)
            bin.numtRNAs <- dim(trnaTab.add)[1]
            tRNAs.tab <- table(trnaTab.add$tRNA_type)
            bin.uniqtRNAs <- length(which(tRNAs.tab > 0))
        }
    }
    ## Summary statistics ##########################################################
    bin.length <- sum(scaff.add$Length)
    bin.numscaffolds <- dim(scaff.add)[1]
    bin.summary <- list (Total_length=bin.length,
                         Num_scaffolds=bin.numscaffolds,
                         Scaff_length_max=max(scaff.add$Length),
                         Scaff_length_min=min(scaff.add$Length),
                         Scaff_length_median=median(scaff.add$Length),
                         Scaff_length_N50=getN50(scaff.add$Length),
                         Marker_sources=marksource,
                         Num_markers=bin.nummarkers,
                         Num_unique_markers=bin.uniqmarkers,
                         Num_singlecopy_markers=bin.singlemarkers,
                         Num_SSUs=bin.numSSUs,
                         Num_tRNAs=bin.numtRNAs,
                         Num_tRNAs_types=bin.uniqtRNAs)
    ## Package and return results ##################################################
    result <- list (scaff=scaff.add,
                    covs=covs.add,
                    markTab=markTab.add,
                    markSource=marksource,
                    ssuTab=ssuTab.add,
                    trnaTab=trnaTab.add,
                    summary=bin.summary,
                    marker.table=marker.tab,
                    tRNA.table=tRNAs.tab,
                    points=NA,
                    slice=NA)
    class(result) <- "gbtbin"
    return(result)
}

summary.gbtbin <- function (x) {
    print(x)  # Print the standard summary
    ## Show the marker tables
    cat ("\n*** Polygon for choosebin (if applicable) ***\n")
    print(x$points)
    cat ("\n*** Table(s) of marker genes ***\n")
    print(x$marker.table)
    cat ("\n*** Table of tRNA genes ***\n")
    print(x$tRNA.table)
}

summary.gbt <- print.gbt  # Identical to "print" behavior

userAdd <- function(x,userTab,userSource) UseMethod("userAdd")

userAdd.gbt <- function(x,
                        userTab,
                        userSource=NA
                        ) {
    ## Check that userTab is data.frame with col "scaffold" ###################
    if (!is.data.frame(userTab) ||
        length(which(names(userTab)=="scaffold"))==0 ||
        is.na(userSource) ) {
        cat("gbtools ERROR: Please check inputs. See help(userAdd) \n")
    } else {
        ## Check that userTab scaffold IDs match x scaffold IDs ###############
        if (length(which(userTab$scaffold %in% x$scaff$ID))==0) {
            cat ("gbtools ERROR: Scaffold IDs in userTab don't match gbt object\n")
        } else {
            x$userTab[[length(x$userTab)+1]] <- userTab  # Append userTab
            x$userSource[length(x$userTab)] <- userSource # Append userSource
            # NB: Using c() will create discrepancy between userTab and userSource
            # because c() on an empty vector will create first element ""
            x$call[[length(x$call)+1]] <- match.call()  # Record function call
            return(x)  # Return result
        }
    }
}

winnow <- function (x, gc, len, covmin, covmax, slice, save, file) UseMethod ("winnow")

winnow.gbt <- function (x,
                        gc=c(0,1),
                        len=c(0,Inf),
                        covmin=NA,
                        covmax=NA,
                        slice=NA,
                        save=FALSE,
                        file="bin_scaffolds.list") {
## "Winnow" a gbt object by GC%, Length, and/or coverage cutoffs
    ## Shortlist scaffolds that match GC and Length criteria #######################
    scafflist <- as.character(x$scaff$ID[which(x$scaff$Ref_GC > gc[1]
                                               & x$scaff$Ref_GC < gc[2]
                                               & x$scaff$Length > len[1]
                                               & x$scaff$Length < len[2])])
    ## Check if coverage cutoffs given, with relevant slices #######################
    if ( is.numeric(slice) ) {
        if (length(covmin)==length(covmax)
            && length(covmin) == length(slice)) {  # Check that values match
            covslist <- as.character(x$covs$ID[which(x$covs[slice[1]+1] > covmin[1]
                                                     & x$covs[slice[1]+1] < covmax[1])])
            if (length(slice)> 1) {
                for (i in 2:length(slice)) {
                    covslist2 <- as.character(x$covs$ID[which(x$covs[slice[i]+1]> covmin[i]
                                                              & x$covs[slice[i]+1]< covmax[i])])
                    covslist <- intersect(covslist,covslist2)
                }
            }
            scafflist <- intersect(scafflist, covslist)  # Update scaffolds shortlist
        } else {
            cat("gbtools ERROR: Lengths of covmin, covmax, and slice parameters
                do not match\n")
        }
    }
    ## Package and return result ###################################################
    bin <- gbtbin(shortlist=scafflist,
                  x=x,
                  slice=NA,
                  taxon=taxon,
                  points=NA,
                  save=save,
                  file=file)
    bin$call[[length(bin$call)+1]] <- match.call()  # Record function call 
    return(bin)
}

winnow.gbtbin <- winnow.gbt  # Inherit behavior of gbt method

winnowMark <- function(x,marksource,param,value,save,file) UseMethod("winnowMark")

winnowMark.gbt <- function(x,  # Object of class gbt
                           marksource=NA, # Which marker set to use
                           param="Class",  # Which taxonomic level to choose?
                           value="Gammaproteobacteria",  # Which taxon to choose?
                           save=FALSE,  # Save list of contigs to external file?
                           file="bin_scaffolds.list"  # File to save list of contigs
                           ) {
## Winnow a gbt object by its marker table values
    if (is.na(marksource)) {
        cat("gbtools WARNING: No marksource supplied, using the first marker set by default...")
        marksource=levels(x$markTab$source)[1]
    }
    markTab.subset <- subset(x$markTab,source==marksource)
    scafflist <- as.character(markTab.subset$scaffold[which(markTab.subset[,which(names(markTab.subset)==param)]==value)])
    bin <- gbtbin (shortlist=scafflist,
                   x=x,
                   points=NA,
                   slice=NA,
                   save=save,
                   file=file)
    bin$call[[length(bin$call)+1]] <- match.call()  # Record function call
    return(bin)
}

winnowMark.gbtbin <- winnowMark.gbt

write.gbt <- function(x,  # Object of class gbt or gbtbin
                      file="gbtools_export.list" # File to export list of contigs
                      ) {
    write(as.character(x$scaff$ID),file=file)
}

write.gbtbin <- write.gbt

importBins <- function(x, # Object of class gbt
                       file, # File with table of bins (1st col) and contig names (2nd col)
                       to.list=T # Create list of gbtbin objects, rather than individual imports
                       ) {
    thetab <- read.table(file=file, sep="\t",header=F)
    names(thetab) <- c("bin","contig")
    binsvector <- as.vector(levels(thetab$bin))
    if (to.list) {
        output <- NULL # Initialize the list to return
        for (i in 1:length(binsvector)) {
            output[[i]] <- gbtbin.default(shortlist=as.vector(thetab$contig[which(thetab$bin==binsvector[i])]), x=x, slice=1)
        }
        names(output) <- binsvector
        return(output)
    }
    else if (!to.list) {
        # Create a dummy function assignfunc() because assign() works by side-effect!
        assignfunc <- function(y) assign(as.character(binsvector[y]),
                                         gbtbin.default(shortlist=as.vector(thetab$contig[which(thetab$bin==binsvector[y])]),
                                                        x=x,
                                                        slice=1
                                                        ),
                                         env=.GlobalEnv
                                         )
        for (i in 1:length(binsvector)) {
            assignfunc(i)
        }
    }
}

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
    # Check that gbt object is really a gbt or gbtbin object, else abandon plotting
    if (class(x) != "gbt" && class(x) != "gbtbin") {
        cat ("\nObject ")
        print (deparse(substitute(x)))
        cat (" is not a gbt or gbtbin object!\n")
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

tabOverlapBins <- function(x, # First list of gbtbin objects 
                           y, # Second list of gbtbin objects
                           binNames.x="", # list of names for x. Default -- number sequentially
                           binNames.y="", # list of names for y
                           weight=TRUE, # weight overlaps by number of bases if true (default), else by number of contigs
                           by="x" # Output type: "x" (default), "y", or "raw" -- by "x" means 'as fraction of rows', and by "y" means 'as fraction of columns'
                           ) {
    # Check that inputs are lists of gbtbin objects
    if (!is.list(x) | !is.list(y)) {
        cat ("gbtools ERROR: The input variables x and y must be character vectors giving names of gbtbin objects \n")
    }
    else {
        if (binNames.x=="" || binNames.y=="") {
            binNames.x <- 1:length(x)
            binNames.y <- 1:length(y)
            binNames.x <- sapply(binNames.x,
                                 function(X) paste("x",X,sep=""))
            binNames.y <- sapply(binNames.y,
                                 function(X) paste("y",X,sep=""))
        }
        if (!weight) {
            counts <- NULL # Initialize the vector of counts x in y
            countsA <- NULL # Initialize vector of fractional counts x in y divided by total of x
            countsB <- NULL # Initialize vector of fractional counts x in y divided by total of y
            for (i in 1:length(x)) {
                for (j in 1:length(y)) {
                    curr.count <- length(which(x[[i]]$scaff$ID %in% y[[j]]$scaff$ID))
                    counts <- c(counts, curr.count) # Append count of overlaps
                    countsA <- c(countsA, curr.count / length(x[[i]]$scaff$ID)) # As fraction of bin x
                    countsB <- c(countsB, curr.count / length(y[[j]]$scaff$ID)) # As fraction of bin y
                    
                }
            }
            counts.matrix <- matrix(counts, byrow=TRUE, nrow=length(x), dimnames=list(binNames.x,binNames.y))
            countsA.matrix <- matrix(countsA, byrow=TRUE, nrow=length(x), dimnames=list(binNames.x,binNames.y))
            countsB.matrix <- matrix(countsB, byrow=TRUE, nrow=length(x), dimnames=list(binNames.x,binNames.y))
            if (by=="x") {
                return (countsA.matrix)
            } else if (by=="y") {
                return (countsB.matrix)
            } else if (by=="raw") {
                return (counts.matrix)
            } else {
                cat ("gbtools ERROR: Invalid value for \"by\" -- should be \"x\", \"y\", or \"raw\" \n" )
            }
        } else {
            lens <- NULL # Initialize the vector of bases x in y
            lensA <- NULL # Initialize the vector of fractional bases x in y divided by total of x
            lensB <- NULL
            for (i in 1:length(x)) {
                for (j in 1:length(y)) {
                    curr.len <- sum(x[[i]]$scaff$Length[which(x[[i]]$scaff$ID %in% y[[j]]$scaff$ID)])
                    lens <- c(lens, curr.len)
                    lensA <- c(lensA, curr.len / sum(x[[i]]$scaff$Length))
                    lensB <- c(lensB, curr.len / sum(y[[j]]$scaff$Length))
                    
                }
            }
            lens.matrix <- matrix(lens,byrow=TRUE,nrow=length(x),dimnames=list(binNames.x,binNames.y))
            lensA.matrix <- matrix(lensA, byrow=TRUE,nrow=length(x),dimnames=list(binNames.x,binNames.y))
            lensB.matrix <- matrix(lensB, byrow=TRUE,nrow=length(x),dimnames=list(binNames.x,binNames.y))
            if (by=="x") {
                return (lensA.matrix)
            } else if (by=="y") {
                return (lensB.matrix) 
            } else if (by=="raw") {
                return(lens.matrix)
            } else {
                cat ("gbtools ERROR: Invalid value for \"by\" -- should be \"x\", \"y\", or \"raw\" \n" )
            }
        }
    }
}

mergeOverlapBins <- function(x, # List of gbtbin objects
                             y, # Second list of gbtbin objects
                             binNames.x="", # Names for gbtbin objects
                             binNames.y="", # Names of gbtbin objects
                             weight=TRUE, # Count overlaps by number of bases (default), else number of contigs
                             by="y", # Condition by fraction in y (default) or x
                             mergeto="x", # Merge y into x (default) or x into y; usually the opposite to "by"
                             threshold="0.8", # Threshold of shared fraction for merging bins
                             out="mergedBin" # Name prefix for output bins
                             ) {
    # Check that inputs are lists of gbtbin objects
    if (!is.list(x) || !is.list(y)) {
        cat ("gbtools ERROR: Input should be lists of gbtbin objects \n")
    } else {
        if (binNames.x=="" || binNames.y=="") {
            binNames.x <- 1:length(x)
            binNames.y <- 1:length(y)
            binNames.x <- sapply(binNames.x,
                                 function(X) paste("x",X,sep=""))
            binNames.y <- sapply(binNames.y,
                                 function(X) paste("y",X,sep=""))
        }
        # Calculate the overlap table
        overlapTab <- tabOverlapBins(x=x,y=y,
                                     binNames.x=binNames.x,
                                     binNames.y=binNames.y,
                                     weight=weight,
                                     by=by)
        # Initialize the list of merged bins
        outList <- NULL
        counter <- 1 # Counter to name the merged bins
        # If merging by "x"
        if (mergeto=="x") {
            for (i in 1:length(x)) {
                hits <- which(overlapTab[i,] > threshold)
                if ( length(hits) > 0 ) {
                    # merge the first y object into the x object
                    newbin <- add(x[[i]],y[[hits[1]]])
                    # If more than one y object exceeds threshold, merge those in sequentially
                    if ( length(hits) > 1) {
                        for (j in 2:length(hits)) {
                            newbin <- add(newbin,y[[hits[j]]])
                        }
                    }
                    # Report which bins were merged
                    cat ("Merged:",
                         as.character(binNames.x[i]),
                         as.character(binNames.y[hits]),
                         "\n",
                         sep=" ")
                    # Add the new bin to list of merged bins
                    outList[[paste(as.character(out),as.character(counter),sep="")]] <- newbin
                    counter <- counter+1
                }
            }
        } else if (mergeto=="y") {   # If merging by "y"
            for (i in 1:length(y)) {
                hits <- which(overlapTab[,i] > threshold)
                if ( length(hits) > 0 ) {
                    # merge the first y object into the x object
                    newbin <- add(y[[i]],x[[hits[1]]])
                    # If more than one y object exceeds threshold, merge those in sequentially
                    if ( length(hits) > 1) {
                        for (j in 2:length(hits)) {
                            newbin <- add(newbin,x[[hits[j]]])
                        }
                    }
                    # Report which bins were merged
                    cat ("Merged:",
                         as.character(binNames.y[i]),
                         as.character(binNames.x[hits]),
                         "\n",
                         sep=" ")
                    # Add the new bin to list of merged bins
                    outList[[paste(as.character(out),as.character(counter),sep="")]] <- newbin
                    counter <- counter+1
                }
            }
        }
        # Return list of gbtbin objects
        return(outList)
    }
}