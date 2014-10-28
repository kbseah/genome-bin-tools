# Interactive metagenome binning tools

Contact: Brandon Seah (kbseah@mpi-bremen.de)

## Introduction

Various tools and approaches exist for metagenomic binning - the process of defining individual genomes in a metagenomic assembly. These tools are designed for interactive exploration and binning of low-diversity microbial metagenomes in R.

A useful way to visualize a metagenomic assembly is to plot the coverage (depth) and GC% of the assembled scaffolds. Scaffolds coming from the same genome would tend to have similar coverage and GC%, and so form clusters in the plots. To aid in distinguishing the clusters, the taxonomic affiliation of each scaffold can be evaluated either by searching the entire scaffold sequence against a database like NCBI nr, or by searching specific marker genes. 

Examples of tools that use GC-coverage plots and taxonomic annotation:
 * Blobology (https://github.com/blaxterlab/blobology)
 * Metawatt (http://sourceforge.net/projects/metawatt/)

Another visualization or binning method relies on the variation in coverage for different genomes between different samples. If the coverage of a metagenome assembly in one sample is plotted against the coverage in another sample, individual genomes would again tend to cluster together.

Examples of tools that use differential coverage binning:
 * Multi-metagenome (http://madsalbertsen.github.io/multi-metagenome/)
 * GroopM (http://minillinim.github.io/GroopM/)

Genome-bin-tools offers:
 * Higher-level functions for plotting - Save time spent on typing and copy-pasting commands
 * Designed to work with free assembly and annotation tools (BBMap for mapping, barrnap for finding rRNAs, tRNAscan-SE for finding tRNAs, AMPHORA2 for finding marker genes)
 * Needing minimal software installation - start R, import some tables, load R functions, and go!
 * Interactive - select bins, see summary statistics for bins, save scaffold lists for later processing

 
 ## Preparing your data
 
 ### Assemble the metagenome
 
 ### Identify marker genes and find phylogenetic affiliation
 
 ### Identify rRNA and tRNA genes (optional)
 
 ### Do a quick preliminary plot
  
 ## Import data to R
 
 ## Explore GC-coverage plots
 
 ## Explore differential coverage plots