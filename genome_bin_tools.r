#!/usr/bin/env R

## Tools for interactive genome binning in R

## Version 3 - 2014-10-28 -- Import tRNA predictions from tRNAscan-SE
## Version 2 - 2014-10-27 -- Imports full taxon string from phylotyping and SSU results. 
## Version 1 - 2014-10-22
## Contact: kbseah@mpi-bremen.de

## Required packages:
##  sp
##  plyr

## Required input files:
##  Output from pileup.sh tool (BBMap package) that was generated with a reference Fasta file (important!)
##  Parsed phylotyping output from AMPHORA or Phyla-AMPHORA

########################################################################################################################
## Single plot tools:

merge.scaff.marker <- function(scaffold.stats,marker.list,taxon,consensus=TRUE) {
## Merge table of scaffold statistics (output from pileup.sh in BBMap package) and table of marker statistics parsed by parse_phylotype_result.pl
## This function needed by other functions in this file
    marker.list[,"taxon"] <- marker.list[,which(names(marker.list)==taxon)]
    marker.stats <- merge(scaffold.stats,marker.list,by.x="ID",by.y="scaffold")
    if (consensus) {    # For scaffolds with multiple marker genes, take majority consensus of marker taxon assignment
        require(plyr)
        #scaffs.with.multi <- as.vector(names(table(marker.stats$ID)[which(table(marker.stats$ID)>1)]))
        #consensus.list <- ddply(marker.list, .(scaffold), function(x) levels(x$taxon)[which.max(tabulate(x$taxon))])
        consensus.list <- ddply(marker.list, .(scaffold), summarize, taxon=levels(taxon)[which.max(tabulate(taxon))])
        marker.stats <- merge(scaffold.stats,consensus.list,by.x="ID",by.y="scaffold")
    }
    return(marker.stats)
}

gc.cov.plot <- function(scaffold.stats,marker.list="",cutoff=0,taxon="Class",consensus=TRUE,legend=FALSE,xlim=c(min(scaffold.stats$Ref_GC),max(scaffold.stats$Ref_GC)),
                        ylim=c(1,max(scaffold.stats$Avg_fold)),log="y",main="",xlab="GC",ylab="Coverage") {
## Given pileup.sh output and parsed phylotyping result, generate GC-cov plot with colored markers if supplied
## Input:
##   scaffold.stats -- Table of scaffold coverage, length, and GC values imported from pileup.sh output
##   marker.list    -- Table of marker statistics parsed by parse_phylotype_result.pl
##   cutoff         -- Don't plot contigs < cutoff (default 0); but if they contain markers, they will be plotted
## Output:
##   Plot of coverage vs. GC, optionally with scaffolds containing markers colored
    if (cutoff>0) {                                                             # If a minimum contig length cutoff is specified, subset the table
        scaffold.stats <- subset(scaffold.stats,Length>=cutoff)
    }
    plot(scaffold.stats$Ref_GC,scaffold.stats$Avg_fold,pch=20,cex=sqrt(scaffold.stats$Length)/100,xlim=xlim,ylim=ylim,log=log,main=main,xlab=xlab,ylab=ylab,col="grey")
    if (marker.list != "" ) {                                                             # If marker list is specified, add markers to plot. TO DO: Choose consensus taxon for scaffolds with >1 marker
        marker.stats <- generate.plot.colors(scaffold.stats,marker.list,taxon, consensus)
        points(marker.stats$Ref_GC,marker.stats$Avg_fold,pch=20,cex=sqrt(marker.stats$Length)/100,col=as.character(marker.stats$colors)) # Add points for scaffolds with marker genes, colored by their taxon
        if (legend) {       # If requested to add a legend to plot
            colorframe <- generate.legend.colors(scaffold.stats,marker.list,taxon,consensus)
            new.colorframe <- subset(colorframe,colors!="grey50")
            newrow <- c("singletons","grey50")
            new.colorframe <- rbind (new.colorframe,newrow)
            legend ("topright",legend=new.colorframe$taxon,cex=0.6,fill=as.character(new.colorframe$colors))
        }
    }
}

generate.plot.colors <- function(scaffold.stats, marker.list, taxon, consensus) {
    marker.stats <- merge.scaff.marker(scaffold.stats,marker.list,taxon, consensus)          # Some table merging to have points to plot for the markers
    marker.list[,"taxon"] <- marker.list[,which(names(marker.list)==taxon)]
    singleton.taxa <- names(table(marker.list$taxon)[which(table(marker.list$taxon)==1)])       # Count how many taxa are only supported by one marker gene
    top.taxon <- names(table(marker.list$taxon)[which.max(table(marker.list$taxon))])           # Which taxon has the most marker genes?
        # Use which.max() because it breaks ties. Otherwise all genomes with same number of marker genes will have same color!
        # Important: Identification of singleton taxa uses the original marker.list because after "consensus",
        #  each scaffold has only one taxon assignment and scaffolds with >1 marker will be undercounted
    ### For plot colors - identify singleton taxa and the taxon with the highest marker counts, and assign them special colors
    taxnames <- names(table(marker.stats$taxon))                                                    # Names of taxa
    taxcolors <- rep("",length(names(table(marker.stats$taxon))))                                   # Create vector to hold color names
    taxcolors[which(names(table(marker.stats$taxon)) %in% singleton.taxa)] <- "grey50"              # Which taxa are singletons? Give them the color "grey50"
    numsingletons <- length(taxcolors[which(names(table(marker.stats$taxon)) %in% singleton.taxa)]) # Count how many singleton taxa
    taxcolors[which(names(table(marker.stats$taxon))==top.taxon)] <- "red"                          # Which taxon has the most marker genes? Give it the color "red"
    numcolors <- length(table(marker.stats$taxon)) - 1 - numsingletons                              # How many other colors do we need, given that all singletons have same color?
    thecolors <- rainbow(numcolors,start=1/6,end=5/6)                                               # Generate needed colors, from yellow to magenta, giving red a wide berth
    taxcolors[which(!(names(table(marker.stats$taxon)) %in% singleton.taxa)  & names(table(marker.stats$taxon))!=top.taxon)] <- thecolors
    colorframe <- data.frame(taxon=taxnames,colors=taxcolors)                                       # Data frame containing which colors correspond to which taxa
    marker.stats <- merge(marker.stats,colorframe,by="taxon")                                       # Merge this by taxon into the marker.stats table for plotting (this works even when consensus option is called)
    return(marker.stats)
}

generate.legend.colors <- function(scaffold.stats, marker.list,taxon, consensus) {
    marker.stats <- merge.scaff.marker(scaffold.stats,marker.list,taxon, consensus)          # Some table merging to have points to plot for the markers
    marker.list[,"taxon"] <- marker.list[,which(names(marker.list)==taxon)]
    singleton.taxa <- names(table(marker.list$taxon)[which(table(marker.list$taxon)==1)])       # Count how many taxa are only supported by one marker gene
    top.taxon <- names(table(marker.list$taxon)[which.max(table(marker.list$taxon))])           # Which taxon has the most marker genes?
    taxnames <- names(table(marker.stats$taxon))                                                    # Names of taxa
    taxcolors <- rep("",length(names(table(marker.stats$taxon))))                                   # Create vector to hold color names
    taxcolors[which(names(table(marker.stats$taxon)) %in% singleton.taxa)] <- "grey50"              # Which taxa are singletons? Give them the color "grey50"
    numsingletons <- length(taxcolors[which(names(table(marker.stats$taxon)) %in% singleton.taxa)]) # Count how many singleton taxa
    taxcolors[which(names(table(marker.stats$taxon))==top.taxon)] <- "red"                          # Which taxon has the most marker genes? Give it the color "red"
    numcolors <- length(table(marker.stats$taxon)) - 1 - numsingletons                              # How many other colors do we need, given that all singletons have same color?
    thecolors <- rainbow(numcolors,start=1/6,end=5/6)                                               # Generate needed colors, from yellow to magenta, giving red a wide berth
    taxcolors[which(!(names(table(marker.stats$taxon)) %in% singleton.taxa)  & names(table(marker.stats$taxon))!=top.taxon)] <- thecolors
    colorframe <- data.frame(taxon=taxnames,colors=taxcolors)                                       # Data frame containing which colors correspond to which taxa
    marker.stats <- merge(marker.stats,colorframe,by="taxon")                                       # Merge this by taxon into the marker.stats table for plotting (this works even when consensus option is called)
    return(colorframe)
}

add.ssu.gc.cov.plot <- function(scaffold.stats,ssu.list,taxon="Order",textlabel=FALSE,pch=10,cex=2,font=2) {
## Mark positions of scaffolds containing SSU rRNA genes on GC-coverage plot, with optional labeling of taxons
    ssu.stats <- merge.scaff.marker(scaffold.stats,ssu.list,taxon,consensus=FALSE)
    points(ssu.stats$Ref_GC,ssu.stats$Avg_fold,pch=pch,cex=cex,col="black")
    if(textlabel==TRUE) {
        text(ssu.stats$Ref_GC,ssu.stats$Avg_fold,as.character(ssu.stats$taxon),pos=3,offset=0.2,font=font)
    }
}

add.tRNA.gc.cov.plot <- function(scaffold.stats,trna.list,pch=4,cex=1) {
## Mark positions of scaffolds containing tRNA genes on GC-coverage plot
    trna.stats <- merge(scaffold.stats,trna.list,by.x="ID",by.y="scaffold")
    points(trna.stats$Ref_GC,trna.stats$Avg_fold,pch=pch,cex=cex,col="black")
}

pick.bin.points <- function(num.points=6,draw.polygon=TRUE) {
## Wrapper for locator() and polygon() to perform interactive binning on the current plot. Returns the polygon vertices which can be used in get.bin.stats()
    thepoints <- locator(num.points,pch=20,type="p")
    if (draw.polygon) { polygon(thepoints) }
    return(thepoints)
}

get.bin.stats <- function(locator.points,scaffold.stats,marker.list="",trna.list="",taxon,save=FALSE,file="interactive_bin.list") {
## Get bin from locator polygon and scaffold statistics
## Input:
##  locator.points -- Points from locator() function used to define genome bin interactively in coverage vs. GC plot,
##  scaffold.stats -- Table of scaffold coverage, length, and GC values imported from pileup.sh output
##  marker.list   -- Table of marker statistics parsed by parse_phylotype_result.pl
##  trna.list -- Table of tRNAs predicted by tRNAscan-SE, imported and headers renamed by import.tRNAscan.results() [optional]
    
    require(sp)
    inpolygon <- point.in.polygon(scaffold.stats$Ref_GC,scaffold.stats$Avg_fold,locator.points$x,locator.points$y)      # Check which scaffolds are contained in that polygon
    scaffold.stats.subset <- scaffold.stats[which(inpolygon==1),]
    
    bin.nummarkers <- NA    # Initialize value of bin.nummarkers for summary, in case the marker.list is not supplied
    bin.uniqmarkers <- NA
    bin.numtRNAs <- NA      # Likewise for number of tRNAs
    bin.uniqtRNAs <- NA
    marker.tab <- NA
    marker.stats.subset <- NA
    tRNAs.tab <- NA
    trna.stats.subset <- NA
    
    if (marker.list != "") {
        marker.stats <- merge.scaff.marker(scaffold.stats,marker.list,taxon,consensus=FALSE)  # Important! "consensus" should be set to FALSE so that all the markers are counted properly
        inpolygon2 <- point.in.polygon(marker.stats$Ref_GC,marker.stats$Avg_fold,locator.points$x,locator.points$y)         # Check which marker-bearing-scaffolds are contained in that polygon
        marker.stats.subset <- marker.stats[which(inpolygon2==1),]
        bin.nummarkers <- dim(marker.stats.subset)[1]                                                                       # Total number of markers in the bin
        marker.tab <- table(marker.stats.subset$gene)                                                                       # Table of counts of each marker that is present (zeroes not shown)
        bin.uniqmarkers <- length(which(marker.tab > 0))                                                                    # Count total number of unique markers
        bin.singlemarkers <- length(which(marker.tab == 1))                                                                 # Count total number of single-copy markers
    }
    
    if (trna.list != "") {
        trna.stats <- merge(scaffold.stats,trna.list,by.x="ID",by.y="scaffold")
        inpolygon3 <- point.in.polygon(trna.stats$Ref_GC,trna.stats$Avg_fold,locator.points$x,locator.points$y)
        trna.stats.subset <- trna.stats[which(inpolygon3==1),]
        bin.numtRNAs <- dim(trna.stats.subset)[1]
        tRNAs.tab <- table(trna.stats.subset$tRNA_type)
        bin.uniqtRNAs <- length(which(tRNAs.tab > 0))                                                                       # Count total number of unique tRNAs
    }
    
    bin.length <- sum(scaffold.stats.subset$Length)                                                                     # Total length of all scaffolds in the bin
    bin.numscaffolds <- dim(scaffold.stats.subset)[1]                                                                   # Total number of scaffolds in the bin
    bin.summary <- data.frame(Total_length=bin.length,Num_scaffolds=bin.numscaffolds,Num_markers=bin.nummarkers,Num_unique_markers=bin.uniqmarkers,Num_singlecopy_markers=bin.singlemarkers,
                              Num_tRNAs=bin.numtRNAs,Num_tRNAs_types=bin.uniqtRNAs)
    if (save) {                                                                                                         # Option to export list of scaffolds contained in this bin, e.g. for reassembly
        write(as.vector(scaffold.stats.subset$ID),file=file)
    }
    
    result <- list(summary=bin.summary,marker.table=marker.tab,tRNA.table=tRNAs.tab,scaff=scaffold.stats.subset,mark=marker.stats.subset,tRNA=trna.stats.subset,points=locator.points)
    return(result)
}

choose.gc.cov.bin <- function(scaffold.stats,marker.list="",trna.list="",taxon="Class",num.points=6,draw.polygon=TRUE,save=FALSE,file="interactive_bin.list") {
## Wrapper for picking out bin on plot interactively and immediately reporting the statistics on scaffolds contained in the bin
    thebin <- pick.bin.points(num.points=num.points,draw.polygon=draw.polygon)
    theresult <- get.bin.stats(locator.points=thebin,scaffold.stats=scaffold.stats,marker.list=marker.list,trna.list=trna.list,taxon=taxon,save=save,file=file)
    return(theresult)
}

import.tRNAscan.results <- function(file) {
    thetable <- read.table(file=file,header=F,sep="\t",skip=3)
    names(thetable) <- c("scaffold","tRNA_no","tRNA_begin","tRNA_end","tRNA_type","Anticodon","Intron_begin","Intron_end","Cove_score")
    return(thetable)
}

########################################################################################################################
## Differential coverage double plot tools:


diffcov.plot <- function(scaffold.stats1,scaffold.stats2,marker.list="",taxon="Class",gc=TRUE,legend=FALSE,cutoff=0,assembly="",
                         consensus=TRUE,xlim=c(min(subset(scaffold.stats1,Avg_fold>0)$Avg_fold),max(scaffold.stats1$Avg_fold)),
                         ylim=c(min(subset(scaffold.stats2,Avg_fold>0)$Avg_fold),max(scaffold.stats2$Avg_fold)),
                         xlab="Coverage 1",ylab="Coverage 2") {
## Generate diffcov plot from pair of pileup.sh outputs, colored by GC content and (optionally) by marker genes taxonomic affiliation
## Input:
##   pileup.sh output of mapping two different read sets onto the SAME reference assembly (scaffold.stats1 and scaffold.stats2 respectively)
##   cutoff -- minimum contig size to plot (default 0)
##   assembly -- Name of assembly to put in plot title
##   gc -- plot the GC plot (default TRUE), if FALSE and the marker.list is supplied, only plot the marker.list. If FALSE and marker.list is not supplied, return something rude
    if (cutoff > 0) {
        scaffold.stats1 <- subset(scaffold.stats1,Length>=cutoff)
        scaffold.stats2 <- subset(scaffold.stats2,Length>=cutoff)
    }
    if (assembly != ""){                                        # If an assembly name is supplied, make it pretty for the plot title
        titlename <- paste("for assembly",assembly,sep=" ")
    }
    else {
        titlename <- ""
    }
    cov1 <- data.frame(ID=scaffold.stats1$ID,cov1=scaffold.stats1$Avg_fold,gc=scaffold.stats1$Ref_GC,len=scaffold.stats1$Length)    # Reformat and merge the coverage information from tables
    cov2 <- data.frame(ID=scaffold.stats2$ID,cov2=scaffold.stats2$Avg_fold)
    diffcov <- merge(cov1,cov2,by="ID")
    gbr <- colorRampPalette(c("green","blue","orange","red"))   # Define colors for GC palette -- from Albertsen script
    if (marker.list != "" & gc==TRUE) {                         # Plot BOTH diffcov with marker genes colored and with GC colored
        #marker.stats <- merge.scaff.marker(diffcov,marker.list,consensus)
        marker.stats <- generate.plot.colors(diffcov,marker.list,taxon,consensus)
        taxa = as.factor(marker.stats$taxon)
        col.taxa = rainbow(length(levels(marker.stats$taxon)))
        palette(col.taxa)
        par(mfrow=c(1,2))
        plot(diffcov$cov1,diffcov$cov2,log="xy",pch=20,cex=sqrt(diffcov$len)/100,col="grey",
             main=paste("Differential coverage plot",titlename,"colored by marker genes",sep=" "),xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim)
        points(marker.stats$cov1,marker.stats$cov2,pch=20,cex=sqrt(marker.stats$len)/100,col=as.character(marker.stats$colors))
        palette(adjustcolor(gbr(70)))
        plot(diffcov$cov1,diffcov$cov2,log="xy",pch=20,cex=sqrt(diffcov$len)/100,col=diffcov$gc*100,
             main=paste("Differential coverage plot",titlename,"colored by GC",sep=" "),xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim)
    }
    else if (marker.list == "" & gc==TRUE) {                    # Plot ONLY GC colored plot
        palette(adjustcolor(gbr(70)))
        plot(diffcov$cov1,diffcov$cov2,log="xy",pch=20,cex=sqrt(diffcov$len)/100,col=diffcov$gc*100,
             main=paste("Differential coverage plot",titlename,"colored by GC",sep=" "),xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim)
    }
    else if (marker.list != "" & gc==FALSE) {                   # Plot ONLY marker genes colored plot
        #marker.stats <- merge.scaff.marker(diffcov,marker.list,consensus)
        marker.stats <- generate.plot.colors(diffcov,marker.list,taxon,consensus)
        taxa = as.factor(marker.stats$taxon)
        col.taxa = rainbow(length(levels(marker.stats$taxon)))
        palette(col.taxa)
        plot(diffcov$cov1,diffcov$cov2,log="xy",pch=20,cex=sqrt(diffcov$len)/100,col="grey",
             main=paste("Differential coverage plot",titlename,"colored by marker genes",sep=" "),xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim)
        points(marker.stats$cov1,marker.stats$cov2,pch=20,cex=sqrt(marker.stats$len)/100,col=as.character(marker.stats$colors))
        palette(adjustcolor(gbr(70)))
        if (legend) {
            colorframe <- generate.legend.colors(diffcov,marker.list,taxon,consensus)
            new.colorframe <- subset(colorframe,colors!="grey50")
            newrow <- c("singletons","grey50")
            new.colorframe <- rbind (new.colorframe,newrow)
            legend ("topright",legend=new.colorframe$taxon,cex=0.6,fill=as.character(new.colorframe$colors))
        }
    }
    else if (marker.list == "" & gc==FALSE) {                   # Plot without any coloring
        plot(diffcov$cov1,diffcov$cov2,log="xy",pch=20,cex=sqrt(diffcov$len)/100,col="grey",
             main=paste("Differential coverage plot",titlename,sep=" "),xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim)
    }
}

add.ssu.diffcov.plot <- function(scaffold.stats1,scaffold.stats2,ssu.list,taxon="Class",textlabel=FALSE,pch=10,cex=2,font=2) {
    cov1 <- data.frame(ID=scaffold.stats1$ID,cov1=scaffold.stats1$Avg_fold,gc=scaffold.stats1$Ref_GC,len=scaffold.stats1$Length)    # Reformat and merge the coverage information from tables
    cov2 <- data.frame(ID=scaffold.stats2$ID,cov2=scaffold.stats2$Avg_fold)
    diffcov <- merge(cov1,cov2,by="ID")
    ssu.stats <- merge.scaff.marker(diffcov,ssu.list,taxon,consensus=FALSE)
    points(ssu.stats$cov1,ssu.stats$cov2,pch=pch,cex=cex,col="black")
    if(textlabel==TRUE) {
        text(ssu.stats$cov1,ssu.stats$cov2,ssu.stats$taxon,cex=0.8,pos=3,offset=0.2,font=font)
    }
}

add.tRNA.diffcov.plot <- function(scaffold.stats1,scaffold.stats2,trna.list,pch=4,cex=1) {
## Mark positions of scaffolds containing tRNA genes on GC-coverage plot
    cov1 <- data.frame(ID=scaffold.stats1$ID,cov1=scaffold.stats1$Avg_fold,gc=scaffold.stats1$Ref_GC,len=scaffold.stats1$Length)    # Reformat and merge the coverage information from tables
    cov2 <- data.frame(ID=scaffold.stats2$ID,cov2=scaffold.stats2$Avg_fold)
    diffcov <- merge(cov1,cov2,by="ID")
    trna.stats <- merge(diffcov,trna.list,by.x="ID",by.y="scaffold")
    points(trna.stats$cov1,trna.stats$cov2,pch=pch,cex=cex,col="black")
}

get.diffcov.bin.stats <- function(locator.points,scaffold.stats1,scaffold.stats2,marker.list="",trna.list="",taxon="Class",save=FALSE,file="interactive_bin.list") {
## Get bin from locator polygon of diffcov plot. Be advised that this has not been tested on a double plot -- use locator on single plot only!
## Input:
##  locator.points -- Points from locator() function used to define genome bin interactively in coverage vs. GC plot,
##  scaffold.stats1 -- Table of scaffold coverage, length, and GC values imported from pileup.sh output for the first sample (X axis)
##  scaffold.stats2 -- Table of scaffold coverage, length, and GC values imported from pileup.sh output for second sample (Y axis)
##  marker.list   -- Table of marker statistics parsed by parse_phylotype_result.pl
##  trna.list     -- Table of tRNAs predicted by tRNAscan-SE, parsed by import.tRNAscan.results()
    cov1 <- data.frame(ID=scaffold.stats1$ID,cov1=scaffold.stats1$Avg_fold,gc=scaffold.stats1$Ref_GC,len=scaffold.stats1$Length)
    cov2 <- data.frame(ID=scaffold.stats2$ID,cov2=scaffold.stats2$Avg_fold)
    diffcov <- merge(cov1,cov2,by="ID")
    require(sp)
    inpolygon <- point.in.polygon(diffcov$cov1,diffcov$cov2,locator.points$x,locator.points$y)               # Check which scaffolds are contained in that polygon
    scaffold.stats.subset <- diffcov[which(inpolygon==1),]
    
    bin.nummarkers <- NA    # Initialize value of bin.nummarkers for summary, in case the marker.list is not supplied
    bin.uniqmarkers <- NA
    bin.numtRNAs <- NA      # Likewise for number of tRNAs
    bin.uniqtRNAs <- NA
    marker.tab <- NA
    marker.stats.subset <- NA
    tRNAs.tab <- NA
    trna.stats.subset <- NA
    
    if (marker.list != "") {
        marker.stats <- merge.scaff.marker(diffcov,marker.list,taxon,consensus=FALSE)           # Important! "consensus" should be set to FALSE so that all the markers are counted properly
        inpolygon2 <- point.in.polygon(marker.stats$cov1,marker.stats$cov2,locator.points$x,locator.points$y)         # Check which marker-bearing-scaffolds are contained in that polygon
        marker.stats.subset <- marker.stats[which(inpolygon2==1),]
        bin.nummarkers <- dim(marker.stats.subset)[1]                                                                       # Total number of markers in the bin
        marker.tab <- table(marker.stats.subset$gene)                                                                       # Table of counts of each marker that is present (zeroes not shown)
        bin.uniqmarkers <- length(which(marker.tab > 0))                                                                    # Count total number of unique markers
        bin.singlemarkers <- length(which(marker.tab == 1))                                                                 # Count total number of single-copy markers
    }
    
    if (trna.list != "") {
        trna.stats <- merge(diffcov,trna.list,by.x="ID",by.y="scaffold")
        inpolygon3 <- point.in.polygon(trna.stats$cov1,trna.stats$cov2,locator.points$x,locator.points$y)
        trna.stats.subset <- trna.stats[which(inpolygon3==1),]
        bin.numtRNAs <- dim(trna.stats.subset)[1]
        tRNAs.tab <- table(trna.stats.subset$tRNA_type)
        bin.uniqtRNAs <- length(which(tRNAs.tab > 0))                                                                       # Count total number of unique tRNAs
    }
    
    bin.length <- sum(scaffold.stats.subset$len)                                                             # Total length of all scaffolds in the bin
    bin.numscaffolds <- dim(scaffold.stats.subset)[1]                                                        # Total number of scaffolds in the bin

    bin.summary <- data.frame(Total_length=bin.length,Num_scaffolds=bin.numscaffolds,Num_markers=bin.nummarkers,Num_unique_markers=bin.uniqmarkers,Num_singlecopy_markers=bin.singlemarkers,
                              Num_tRNAs=bin.numtRNAs,Num_tRNAs_types=bin.uniqtRNAs)
    if (save) {                                                                                              # Option to export list of scaffolds contained in this bin, e.g. for reassembly
        write(as.vector(scaffold.stats.subset$ID),file=file)
    }
    result <- list(summary=bin.summary,marker.table=marker.tab,tRNA.table=tRNAs.tab,scaff=scaffold.stats.subset,mark=marker.stats.subset,tRNA=trna.stats.subset,points=locator.points)
    return(result)
}

choose.diffcov.bin <- function(scaffold.stats1,scaffold.stats2,marker.list="",trna.list="",num.points=6,draw.polygon=TRUE,save=FALSE,file="interactive_bin.list") {
## Wrapper for picking out bin on plot interactively and immediately reporting the statistics on scaffolds contained in the bin
    thebin <- pick.bin.points(num.points=num.points,draw.polygon=draw.polygon)
    theresult <- get.diffcov.bin.stats(locator.points=thebin,scaffold.stats1=scaffold.stats1,scaffold.stats2=scaffold.stats2,marker.list=marker.list,trna.list=trna.list,save=save,file=file)
    return(theresult)
}