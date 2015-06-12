#!/usr/bin/env R

## Version 1 - 2014-10-22
## Contact: kbseah@mpi-bremen.de

## R Script to parse contig coverage info from output of pileup.sh (BBMap package)
## and phylotyping results parsed from AMPHORA2 or Phyla-AMPHORA pipelines
## and output pretty plots of coverage vs. GC, colored by marker genes
## !! This script is designed to be called by assembly_parse_plot.pl, which does preliminary formatting !!

## Usage:
## $ Rscript assembly_parse_plot.r <coverage_table> <parsed_phylotyping_result> <assembly_name>

## Read commandline arguments
args <- commandArgs(trailingOnly=TRUE)
cov.file <- args[2]     # File containing coverage, GC% and length info for each contig
phy.file <- args[3]     # File containing parsed phylotyping results
assem.name <- args[4]   # Name of the genomic assembly / project number
rm (args)

## Import the tables
cov.tab <- read.table(cov.file,header=T,sep="\t")
phy.tab <- read.table(phy.file,header=T,sep="\t")

## Subset the scaffold coverage table for plotting
cov.tab.reduced <- subset(cov.tab,Length>200)       # Taking only contigs with length > 200
cov.tab.withmarkers <- merge(cov.tab,phy.tab,by.x="ID",by.y="scaffold")

# Plain coverage vs. GC plot
pdf(paste(assem.name,"gccov_plot","plain","pdf",sep="."),width=10,height=7)
plot(cov.tab.reduced$Ref_GC, cov.tab.reduced$Avg_fold, pch=20, cex=sqrt(cov.tab.reduced$Length)/100, col="grey", log="y", main=paste("Plot of assembly",assem.name,sep=" "), xlab="GC", ylab="Coverage")
dev.off()

png(paste(assem.name,"gccov_plot","plain","png",sep="."),width=1000,height=700)
plot(cov.tab.reduced$Ref_GC, cov.tab.reduced$Avg_fold, pch=20, cex=sqrt(cov.tab.reduced$Length)/100, col="grey", log="y", main=paste("Plot of assembly",assem.name,sep=" "), xlab="GC", ylab="Coverage")
dev.off()

# Scaffolds containing marker genes colored by taxonomic affiliation
taxa <- as.factor(phy.tab$taxon)                        # Convert taxa to factors for coloring
col.taxa <- rainbow(length(levels(phy.tab$taxon)))      # Define the palette
palette(col.taxa)
pdf(paste(assem.name,"gccov_plot","markers","pdf",sep="."),width=10,height=7)
plot(cov.tab.reduced$Ref_GC, cov.tab.reduced$Avg_fold, pch=20, cex=sqrt(cov.tab.reduced$Length)/100, col="grey", log="y", main=paste("Plot of assembly",assem.name,sep=" "), xlab="GC", ylab="Coverage")
points(cov.tab.withmarkers$Ref_GC,cov.tab.withmarkers$Avg_fold, pch=20, cex=sqrt(cov.tab.withmarkers$Length)/100, col=taxa)
dev.off()

png(paste(assem.name,"gccov_plot","markers","png",sep="."),width=1000,height=700)
plot(cov.tab.reduced$Ref_GC, cov.tab.reduced$Avg_fold, pch=20, cex=sqrt(cov.tab.reduced$Length)/100, col="grey", log="y", main=paste("Plot of assembly",assem.name,sep=" "), xlab="GC", ylab="Coverage")
points(cov.tab.withmarkers$Ref_GC,cov.tab.withmarkers$Avg_fold, pch=20, cex=sqrt(cov.tab.withmarkers$Length)/100, col=taxa)
dev.off()