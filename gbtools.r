#' Take union of two gbtbin objects
#'
#' Take union of two gbtbin objects. Equivalent to the R union function
#'
#' Self explanatory...
#' 
#' @param x1 Object of class gbtbin
#' @param x2 Object of class gbtbin
#'
#' @return Object of class gbtbin
#'
#' @seealso \code{\link{lej}}
#' @export
add <- function(x1, x2) UseMethod("add")
add.gbtbin <- function(x1,x2) {
## Merge two bins; i.e. take their union
    result <- setOperation(x1=x1,x2=x2,shortlist="all")
    result$call <- match.call()  # Record function call that returned this merged bin
    return(result)
}#' Choose bin interactively from plot of gbt object
#'
#' Choose genome bin from GC-coverage or differential coverage plot of a
#' gbt object
#'
#' @param x Object of class gbt, used to generate the plot
#' @param slice The same slice parameter used to generate the plot
#' @inheritParams pickBinPoints
#'
#' @importFrom sp point.in.polygon
#' @return Object of class gbtbin
#' @seealso \code{\link{plot.gbt}}
#' @export
#'
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
        result$call <- match.call()  # Record the choosebin() call used to produce this bin
        return(result)
    }
}
#' Tabulate objects and count how many singletons
#'
#' @param x Object of class data.frame or vector
#' @return Numeric vector of length 1
#' @keywords internal
countSingleFromTable <- function(x) {
    x.tab <- table(x)
    uniq <- length(which(x.tab==1))
    return(uniq)
}
#' Perform connectivity fishing with Fastg file
#'
#' Find contigs connected to existing gbtbin object, using connectivity from
#' an external Fastg file
#'
#' @param x Object of class gbt, from which the bin object was defined
#' @param bin Object of class gbtbin, defined from the gbt object x
#' @param fastg.file Path to Fastg file for the metagenome assembly of x
#' @param save Logical - Save list of contigs of new bin to external file?
#' @param file File to save contig list 
#' @return Object of class gbtbin
#' @export
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
    newbin <- gbtbin(shortlist=fished.contigs.list,x=x,slice=NA,taxon=taxon,save=save,file=file)
    return(newbin)
}

#' Create gbtbin from gbt object
#'
#' Creates object of class gbtbin from existing gbt object, representing a
#' genome bin defined from a metagenome assembly.
#'
#' This function is called by \code{\link{choosebin}}, but can be called
#' directly by the user, given a shortlist of contigs in the bin.
#'
#' @param shortlist List of contigs in the new bin (required)
#' @param x Object of class gbt
#' @param slice Numeric vector describing which coverage values were used to
#'               draw the plot from which the bin was selected by choosebin.
#'               Use NA if not calling this function from choosebin.
#' @param taxon For internal use
#' @param points Numeric vector describing points used to select the bin from
#'                a plot by choosebin. Use NA if not calling this function
#'                from choosebin.
#' @param save Save list of contigs to external file? (logical, default FALSE)
#' @param file File name to save contigs, if save=TRUE.
#'
#' @return Object of class gbtbin
#' @seealso \code{\link{gbt}}, \code{\link{choosebin}}
#' @export
gbtbin <- function(shortlist,x,slice,taxon,points,save,file) UseMethod("gbtbin")
gbtbin.default <- function(shortlist,  # Character vector, contigs to extract from gbt object
                           x,  # Object of class gbt
                           slice,  # Which slice used by choosebin()?
                           taxon,  # Deprecated - user don't change
                           points=NA,  # Number of points in polygon for choosebin()
                           save=FALSE,  # Save contig list to external file?
                           file="interactive_bin.list"  # File name to save contig list
                           ) {
    scaff.subset <- subset(x$scaff, ID%in% shortlist)
    covs.subset <- subset(x$covs, ID%in%shortlist)
    markTab.subset <- NA
    marksource <- x$markSource
    ssuTab.subset <- NA
    trnaTab.subset <- NA
    # Summary statistics initialize
    bin.nummarkers <- NA
    bin.uniqmarkers <- NA
    bin.numtRNAs <- NA
    bin.uniqtRNAs <- NA
    bin.numSSUs <- NA
    bin.singlemarkers <- NA
    marker.tab <- list()
    tRNAs.tab <- NA
    ## Take subset of markTab ##############################
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
    ## Take subset of ssuTab ##############################
    if (is.data.frame(x$ssuTab)) {
        ssuTab.subset <- subset(x$ssuTab,scaffold%in%shortlist)
        bin.numSSUs <- dim(ssuTab.subset)[1]
    }
    ## Take subset of trnaTab ##############################
    if (is.data.frame(x$trnaTab)) {
        trnaTab.subset <- subset(x$trnaTab,scaffold%in% shortlist)
        bin.numtRNAs <- dim(trnaTab.subset)[1]
        tRNAs.tab <- table (trnaTab.subset$tRNA_type)
        bin.uniqtRNAs <- length(which(tRNAs.tab > 0))
    }
    ## Summary statistics ######################################
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
    
    ## Write to file, if save option is used ######################
    if (save) {
        write(as.vector(scaff.subset$ID),file=file)
    }
    
    ## Package and return result ######################################
    result <- list(scaff=scaff.subset,
                   covs=covs.subset,
                   markTab=markTab.subset,
                   markSource=marksource,
                   ssuTab=ssuTab.subset,
                   trnaTab=trnaTab.subset,
                   summary=bin.summary,
                   marker.table=marker.tab,
                   tRNA.table=tRNAs.tab,
                   points=points,
                   slice=slice)
    class(result) <- "gbtbin"
    return(result)
}
#' Interactive visualization of metagenome assemblies in R
#'
#' gbtools lets you visualize metagenome assemblies by GC-coverage or 
#' differential coverage plots, and interactively bin them. 
#'
#' See website for more details.
#'
#' @references \url{https://github.com/kbseah/genome-bin-tools}
#' @author Brandon Seah, \email{kbseah@@mpi-bremen.de}
#' @docType package
#' @name gbtools
NULL#' Create new gbt object
#'
#' Creates new gbt object from coverage data, taxonomic markers, and other
#' data
#'
#' See documentation online https://github.com/kbseah/genome-bin-tools for
#' fuller instructions on generating the input files required.
#'
#' @param covstats File(s) with coverage statistics of a metagenomic assembly;
#'                  output of pileup.sh in BBTools package (required). More
#'                  than one file can be specified with c().
#' @param mark Table of scaffolds with marker genes and taxonomic information
#'              for each marker. E.g. AMPHORA2 or Phyla-AMPHORA marker sets
#'              and output parsed by parse_phylotype_result.pl. (optional)
#' @param ssu Table of scaffolds with SSU rRNA genes, and taxonomic info
#'             for each SSU. E.g. use get_ssu_for_genome_bin_tools.pl.
#'             (optional)
#' @param tnra Table of tRNA genes found in assembly. Can use the output from
#'              tRNAscan-SE directly. (optional)
#'
#' @return Object of class gbt
#'
#' @seealso \code{\link{gbtbin}}, \code{link{choosebin}}
#'
#' @export
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
            cat ("gbtools WARNING: marksource not supplied.\n")
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
        summarystats <- list(Total_length=sum(scaff$Length),
                             Num_scaffolds=length(scaff$ID),
                             Scaff_length_max=max(scaff$Length),
                             Scaff_length_min=min(scaff$Length),
                             Scaff_length_median=median(scaff$Length),
                             Scaff_length_N50=getN50(scaff$Length),
                             Num_markers=table(markTab$source),
                             Num_SSU=numSsu,
                             Num_tRNAs=numTrna)
        
        ## Package and return result #########################################
        result <- list(scaff=scaff,
                       covs=covs,
                       markTab=markTab,
                       markSource=marksource,
                       ssuTab=ssuTab,
                       trnaTab=trnaTab,
                       summary=summarystats)
        result$call <- match.call()  # Record function call that produces this gbt object
        class(result) <- "gbt"
        result
    }    
}#' Generates colors for plot legends when coloring by markers
#'
#' @param scaffold.stats scaff data.frame from a gbt object
#' @param marker.list List of markers from markTab data.frame of gbt object
#' @param taxon Taxon names for which to generate colors
#' @param consensus If scaffold has more than one marker, take majority-rule
#'                    consensus taxonomic assignment? (logical)
#' @return data.frame with color assignments for each scaffold
#' @keywords internal 
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
}#' Generates colors for marker gene phylotypes in plot
#'
#' @inheritParams generateLegendColors
#' @return data.frame with color assignments for each scaffold
#' @keywords internal
#'
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
}#' Calculate N50 from contig lengths
#'
#' @param x Numeric vector
#' @return N50 value - min contig length that contains half of total bases
#' @keywords internal
getN50 <- function(x) {
    # Adapted from R-bloggers:
    # http://www.r-bloggers.com/calculating-an-n50-from-velvet-output/
    x.sort <- sort(x,decreasing=TRUE)
    n50 <- x[cumsum(x.sort) >= sum(x.sort)/2][1]
    return(n50)
}#' Take difference between two gbtbin objects
#'
#' Takes the reverse complement of two gbtbin objects. Equivalent to setdiff
#' in R, or left-exclusive-join in SQL. Non commutative!
#'
#' Self explanatory...
#'
#' @inheritParams add
#'
#' @seealso \code{\link{add}}
#'
#' @export
lej <- function(x1, x2) UseMethod ("lej")
lej.gbtbin <- function(x1,x2) {
## Take difference between two bins - non commutative! i.e. left exclusive join
    shortlist <- x1$scaff$ID[which(!x1$scaff$ID %in% x2$scaff$ID)]
    result <- setOperation.gbtbin(x1=x1,
                                  x2=x2,
                                  shortlist=shortlist)
    result$call <- match.call()  # Record function call that returned this merged bin
    return(result)
}
#' Merge scaffolds table with marker table
#'
#' Merge table of scaffold statistics (output from pileup.sh in BBMap package)
#' and table of marker statistics parsed by parse_phylotype_result.pl
#'
#' @inheritParams generateLegendColors
#'
#' @return data.frame with marker data and scaffold statistics
#' @keywords internal
#' @importFrom plyr ddply
#'
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
}#' Define polygon for genome bin on plot
#' 
#' Wrapper for locator() and polygon() to perform interactive binning on the
#' current plot. Returns the polygon vertices which can be used in
#' get.bin.stats()
#'
#' @param num.points Number of points in polygon (integer)
#' @param draw.polygon Draw polygon as overlay on plot (logical)
#'
#' @return numerical vector containing polygon vertices coordinates
#' @keywords internal
#'
pickBinPoints <- function(num.points=6,  # How many points in polygon?
                          draw.polygon=TRUE  # Overlay polygon on plot?
                          ) {
## Wrapper for locator() and polygon() to perform interactive binning on the current
## plot. Returns the polygon vertices which can be used in get.bin.stats()
    thepoints <- locator(num.points,pch=20,type="p")
    if (draw.polygon) { polygon(thepoints) }
    return(thepoints)
}#' Plot object of class gbtbin
#'
#' Plot GC-coverage or differential coverage plots from a gbtbin object
#'
#' See documentation on \code{\link{plot.gbt}} for more details on plotting
#' @inheritParams plot.gbt
#' @seealso \code{\link{plot.gbt}}, \code{\link{points.gbtbin}}
#' @return New graphics window and a plot
#' @export
#'
plot.gbtbin <- function(x, slice="default", ...) {
    ## Defaults to same slice used to choose the bin ################################
    if (slice == "default") { 
        slice <- x$slice
    }
    ## Catch invalid slice values ###################################################
    if (is.na(x$slice) || !is.numeric(x$slice) || length(x$slice) > 2) { 
        cat ("gbtools ERROR: Please specify valid value for slice option for this bin\n")
    }
    ## Else inherit same plot method as gbt class, for simplicty's sake #############
    else {
        plot.gbt (x=x, slice=slice, ...) 
    }
}#' Plot object of class gbt
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
}#' Add points to plot of gbt or gbtbin object
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
}#' Print summary of gbt objectbin
#'
#' @param x Object of class gbtbin
#' @return data.frame summarizing contents of x
#' @seealso \code{\link{gbtbin}}, \code{\link{summary.gbtbin}}
#' @export
print.gbtbin <- function(x) {
    cat("Object of class gbtbin\n")
    cat ("\nScaffolds:\n")
    lengthdf <- data.frame(x$summary$Total_length,
                           x$summary$Num_scaffolds,
                           x$summary$Scaff_length_max,
                           x$summary$Scaff_length_min,
                           x$summary$Scaff_length_median,
                           x$summary$Scaff_length_N50)
    names(lengthdf) <- c("Total", "Scaffolds", "Max", "Min", "Median", "N50")
    print(lengthdf)
    cat("\nMarkers:\n")
    markerdf <- data.frame(x$summary$Marker_sources,
                           x$summary$Num_markers,
                           x$summary$Num_unique_markers,
                           x$summary$Num_singlecopy_markers)
    names(markerdf) <- c("Source","Total","Unique","Singlecopy")
    print(markerdf)
    cat("\nSSU markers:\n")
    print(x$summary$Num_SSUs)
    cat("\ntRNA_markers:\n")
    print(x$summary$Num_tRNAs)
    cat("\nCall:\t")
    print(x$call)
}
#' Print summary of gbt object
#'
#' @param x Object of class gbt
#' @return data.frame summarizing contents of x
#' @seealso \code{\link{gbt}}, \code{\link{summary.gbt}}
#' @export
print.gbt <- function(x) {
    cat("Object of class gbt\n")
    cat("\nScaffolds:\n")
    lengthdf <- data.frame(x$summary$Total_length,
                           x$summary$Num_scaffolds,
                           x$summary$Scaff_length_max,
                           x$summary$Scaff_length_min,
                           x$summary$Scaff_length_median,
                           x$summary$Scaff_length_N50)
    names(lengthdf) <- c("Total", "Scaffolds", "Max", "Min", "Median", "N50")
    print(lengthdf)
    cat("\nMarker counts by source:")
    print(x$summary$Num_markers)
    cat("\nSSU markers:\n")
    print(x$summary$Num_SSU)
    cat("\ntRNA markers:\n")
    print(x$summary$Num_tRNAs)
    cat("\nCall:\t")
    print(x$call)
}#' Generic operation for merging and subsetting two gbtbin objects
#' 
#' @param x1 Object of class gbtbin
#' @param x2 Object of class gbtbin
#' @param shortlist Vector of contig IDs to make new bin
#' @return Object of class gbtbin
#' @keywords internal

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
    if (!is.na(x1$markTab) || !is.na(x2$markTab)) {
        markTab.add <- unique (rbind(x1$markTab,x2$markTab))
        for (j in 1:length(marksource)) {
            markTab.subsubset <- subset(markTab.add,source==marksource[j])
            bin.nummarkers[j] <- dim(markTab.subsubset)[1]
            marker.tab[[j]] <- table(markTab.subsubset$gene)
            bin.uniqmarkers[j] <- length(which(marker.tab[[j]] > 0))
            bin.singlemarkers[j] <- length(which(marker.tab[[j]] ==1))
        }
    }
    ## Combine SSU Tables ##########################################################
    if (!is.na(x1$ssuTab) || !is.na(x2$ssuTab)) {
        ssuTab.add <- unique (rbind(x1$ssuTab,x2$ssuTab))
        bin.numSSUs <- dim(ssuTab.add)[1]
    }
    ## Combine trna Tables #########################################################
    if (!is.na(x1$trnaTab) || !is.na(x2$trnaTab)) {
        trnaTab.add <- unique (rbind(x1$trnaTab,x2$trnaTab))
        bin.numtRNAs <- dim(trnaTab.add)[1]
        tRNAs.tab <- table(trnaTab.add$tRNA_type)
        bin.uniqtRNAs <- length(which(tRNAs.tab > 0))
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
}#' Print summary of gbtbin object
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
}#' Print summary of gbt object
#'
#' @param x Object of class gbt
#' @return data.frame summarizing contents of x
#' @seealso \code{\link{gbt}}, \code{\link{print.gbt}}
#' @export
summary.gbt <- print.gbt  # Identical to "print" behavior
#' Subset a gbt or gbtbin object by marker gene taxonomy
#'
#' @param param Taxonomic level to make subset (default "Class")
#' @param value Value of the taxon to make subset (default
#'               "Gammaproteobacteria")
#' @inheritParams winnow
#' @return Object of class gbtbin
#' @seealso \code{\link{winnow}}, \code{\link{gbt}}
#' @export
#'
winnowMark <- function(x,param,value,save,file) UseMethod("winnowMark")
winnowMark.gbt <- function(x,  # Object of class gbt
                           param="Class",  # Which taxonomic level to choose?
                           value="Gammaproteobacteria",  # Which taxon to choose?
                           save=FALSE,  # Save list of contigs to external file?
                           file="bin_scaffolds.list"  # File to save list of contigs
                           ) {
## Winnow a gbt object by its marker table values
    scafflist <- as.character(x$markTab$scaffold[which(x$markTab[,which(names(x$markTab)
                                                                        ==param)]
                                                       ==value)])
    bin <- gbtbin (shortlist=scafflist,
                   x=x,
                   points=NA,
                   slice=NA,
                   save=save,
                   file=file)
    bin$call <- match.call()  # Record function call that produced this winnow bin
    return(bin)
}
winnowMark.gbtbin <- winnowMark.gbt#' Subset a gbt or gbtbin object by GC, length, or coverage cutoffs
#'
#' @param x Object of class gbt or gbtbin
#' @param gc Vector giving min and max GC values (default c(0,Inf))
#' @param len Vector giving min and max contig lengths (default c(0,Inf))
#' @param covmin Vector giving minimum coverage values per sample (default NA)
#' @param covmax Vector giving maximum coverage values per sample (default NA)
#' @param slice Vector of sample numbers for the coverage cutoffs (default NA)
#' @param save Save list of contigs to external file? (logical, defautl FALSE)
#' @param file File name for export of contig list.
#' @return Object of class gbtbin
#' @seealso \code{\link{winnowMark}}, \code{\link{gbt}}
#' @export
#'
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
    bin$call <- match.call()  # Record function call that produced this winnow bin
    return(bin)
}
winnow.gbtbin <- winnow.gbt  # Inherit behavior of gbt method