# Interactive metagenome binning tools

Contact: Brandon Seah (kbseah@mpi-bremen.de)

Cite: Brandon Seah (2014), genome-bin-tools, Online: https://github.com/kbseah/genome-bin-tools

## 0. Introduction

Various tools and approaches exist for metagenomic binning - the process of defining individual genomes in a metagenomic assembly. These tools are designed for interactive exploration and binning of low-diversity microbial metagenomes in R.

A useful way to visualize a metagenomic assembly is to plot the coverage (depth) and GC% of the assembled scaffolds. Scaffolds coming from the same genome would tend to have similar coverage and GC%, and so form clusters in the plots. To aid in distinguishing the clusters, the taxonomic affiliation of each scaffold can be evaluated either by searching the entire scaffold sequence against a database like NCBI nr, or by searching specific marker genes. 

Examples of tools that use GC-coverage plots and taxonomic annotation include Blobology (https://github.com/blaxterlab/blobology) and Metawatt (http://sourceforge.net/projects/metawatt/).

Another visualization or binning method relies on the variation in coverage for different genomes between different samples. If the coverage of a metagenome assembly in one sample is plotted against the coverage in another sample, individual genomes would again tend to cluster together.

Examples of tools that use differential coverage binning: Multi-metagenome (http://madsalbertsen.github.io/multi-metagenome/), GroopM (http://minillinim.github.io/GroopM/).

Genome-bin-tools builds on concepts from Multi-metagenome, but it offers more:
 * Higher-level functions for plotting - Save time spent on typing and copy-pasting commands
 * Designed to work with free assembly and annotation tools (BBMap for mapping, barrnap for finding rRNAs, tRNAscan-SE for finding tRNAs, AMPHORA2 for finding marker genes)
 * Needing minimal software installation - start R, import some tables, load R functions, and go!
 * Interactive - select bins, see summary statistics for bins, save scaffold lists for later processing

## 1. Produce and annotate metagenomic assembly

If you want to follow along, you can use the data from Albertsen et al. 2013, which is in the folder `example_data`. The file names in the following examples use the file names of the example data. (NB: The original header names of the assembly Fasta file were edited to be shorter.)

### 1a. Assemble the metagenome and calculate coverage

Use your favorite assemblier, like IDBA-UD (http://i.cs.hku.hk/~alse/hkubrg/projects/idba/) or SPAdes (http://bioinf.spbau.ru/spades), to assemble your metagenome. The assembly should be in a Fasta file. Calculate coverage with bbmap.sh (http://sourceforge.net/projects/bbmap/) by mapping the original reads used to assemble the metagenome back onto the assembly:

```
 $ bbmap.sh ref=HPminus_assembly.fasta nodisk in=HPminus_reads.fq.gz covstats=HPminus.coverage
```

BBmap is convenient because it calculates length and GC along with the coverage per scaffold, and outputs them all to the same file. If you use a different mapping tool, you will have to aggregate all these data together in a tab-separated table with the same header names as the bbmap output.

For differential coverage binning, you will need a second read library from a different sample where at least some of the genomes present in the first sample should also be present. Map the second read library onto the same assembly and save the coverage values as a different file:

```
 $ bbmap.sh ref=HPminus_assembly.fasta nodisk in=HPplus_reads.fq.gz covstats=HPplus.coverage
```

### 1b. Identify marker genes and find phylogenetic affiliation (optional)

Use AMPHORA2 (https://github.com/martinwu/AMPHORA2) or Phyla-AMPHORA (https://github.com/martinwu/Phyla_AMPHORA) to identify conserved marker genes in your assembly, and to assign a taxonomic position. Parse the output of their script Phylotyping.pl (here called `phylotype.result`) for import into R:

```
 $ perl parse_phylotype_result.pl -p phylotype.result > phylotype.result.parsed
```

This generates a file called `phylotype.result.parsed` which will be imported into R.

### 1c. Identify rRNA genes (optional)

Use barrnap (http://www.vicbioinformatics.com/software.barrnap.shtml) to detect SSU rRNA genes in the assembly, and assign phylotype by extacting sequences using fastaFrombed (http://bedtools.readthedocs.org/en/latest/) and then using Usearch (http://www.drive5.com/usearch/) against a curated SILVA (www.arb-silva.de/) database. The database has to be prepared in a specific way (instructions to come) but is identical to the Usearch-indexed database used by PhyloFlash (https://bitbucket.org/HGV/phyloflash.git). PhyloFlash is also a great tool, why not check it out? (Disclosure: I helped to develop PhyloFlash).

The rRNA extraction and output parsing is done with a wrapper script:
```
 $ perl get_ssu_for_genome_bin_tools.pl -d <path/to/ssu/database> -c <number_CPUs> -a HPminus_assembly.fasta -o <output_prefix> 
```

This generates a file called <output_prefix>.ssu.tab (in the example data, HPminus.ssu.tab) which will be imported into R.

### 1d. Identify tRNA genes (optional)

Use tRNAscan-SE version 1.23 (http://selab.janelia.org/tRNAscan-SE/) to find tRNA genes. 

```
 $ tRNAscan-SE -G -o HPminus.trna.tab HPminus_assembly.fasta
```

The output HPminus.trna.tab is directly imported into R.

### 1e. Do a quick preliminary plot

(To be updated)

## 2. Load functions into R

Start R. Required packages are `sp` and `plyr`, which can be installed like so:

```R
 > install.packages("sp")
 > install.packages("plyr")
```

Install the genome.bin.tools package in R:

```R
 > install.packages("/path/to/genome.bin.tools_1.0.tar.gz",repos=NULL,type="source")
```

This is recommended because you can call help() to read the documentation for each function within the R environment.

Alternatively, load the R functions with source (recommended if you want to tweak them or use experimental features):

```R
 > source("genome_bin_tools.r")
```

## 3. Explore GC-coverage plots

GC-coverage plots are generated from single samples (i.e. coverage statistics from mapping a single read library onto a single assembly).

The data are imported as objects of class `genomestats` using the function of the same name:

```R
 > d <- genomestats(covstats="HPminus.coverage",marker.list="phylotype.result.parsed",ssu.list="HPminus.ssu.tab",trna.list="HPminus.trna.tab")
```

These import and parse the files generated earlier. Only the coverage file `assembly.coverage` is required; the rest are optional (though having them will make more informative and useful plots).

Type the name of the object to see a summary:

```R
 > d
```

### 3a. Plotting

```R
 > plot(d) # basic plot, if marker.list was imported, then colored automatically by marker taxonomy
 > plot(d,cutoff=2000) # Do not show scaffolds shorter than 2000 bp
 > plot(d,taxon="Phylum") # Color marker genes at the taxonomic level of "Phylum"
 > plot(d,ssu=TRUE) # Overlay crosshairs marking scaffolds that contain SSU rRNA genes
 > plot(d,ssu=TRUE,textlabel=TRUE) # Add text labels beside crosshairs showing phylotype assigned to SSU rRNA genes
 > plot(d,trna=TRUE) # Add cross marks showing scaffolds that contain tRNA rRNA genes
 > plot(d,legend=TRUE) # Add a legend showing which color corresponds to which taxonomic group; grey points are singleton taxa
```

Zoom into specific areas of the plot by altering the `xlim` and `ylim` parameters, as with the basic `plot` function in R.

### 3b. Interactively choosing genomic bins

If you see a cluster of scaffolds which you would like to save as a bin, you can choose it interactively by picking the points that draw a polygon surrounding the scaffolds you want:

```R
 > bin1 <- choosebin.genomestats(d) # basic function
 > bin1 <- choosebin.genomestats(d, save=TRUE, file="bin1.scaffolds.list") # Save the names of scaffolds in this bin to an external file called bin1.scaffolds.list
 > bin1 <- choosebin.genomestats(d,num.points=10) # Change the number of corners of the polygon
```

`bin` is now an object of class genomestatsbin. Type the name of the bin object to see a summary of the bin. If you imported the marker, SSU, and/or tRNA data, a summary of how many of each are contained in the bin will be reported (this is useful if the marker genes are typical single-copy genes, for example).

### 3c. Fishing for connected contigs using Fastg files (experimental)

Fastg files are generated by newer versions of the SPAdes assembler, and contain contig connectivity information generated during the assembly process. These can be useful, e.g. to "fish" scaffolds that are from the same genome but which were inadvertently left out of the interactively chosen bin.

This needs a genomestatsbin object and produces another genomestatsbin object.

For now you will have to manually edit the `genome_bin_tools.r` file to specify the location of the script `fastg_parser.pl` (included with this package) and also a directory to store temporary files (default: /tmp/). 

Fish out a new bin from an old bin using a Fastg file:
```R
 > bin2 <- fastg_fishing.genomestatsbin(d,bin1,"path/to/fastg/file.fastg")
 > bin2 <- fastg_fishing.genomestatsbin(d,bin1,"path/to/fastg/file.fastg", save=TRUE,file="bin2.scaffolds.list") # Save list of scaffolds in new bin to external file
```

You can compare the two bins by plotting them overlaid:
```R
 > plot(d,marker=FALSE) # The underlying GC-coverage plot, with coloring turned off
 > points(bin2,col="blue") # The new bin in blue
 > points(bin1,col="black") # The original bin in black
```

## 4. Explore differential coverage plots

Differential coverage plots are generated from two separate coverage files (in this example: `assembly.coverage` and `assembly.coverage2`). The tools use analogous object classes `diffcovstats` and `diffcovstatsbin` (for bins defined from differential coverage plots).

Import data:

```R
 > D <- diffcovstats(covstats1="assembly.coverage",covstats2="assembly.coverage2",marker.list="phylotype.results.parsed",ssu.list="assembly.ssu.tab",trna.list="trnascan.results")
```

Type the object name to see a summary. 

### 4a. Plotting

Analogous to plotting for genomestats objects, but with option of coloring by markers or by GC value:

```R
 > plot(D) # Basic plot. Defaults to coloring by marker genes, if data imported
 > plot(D,marker=FALSE,gc=FALSE) # Uncolored plot
 > plot(D,marker=TRUE,legend=TRUE) # Add legend
 > plot(D,gc=TRUE,marker=FALSE) # Color by GC
 > plot(D,gc=TRUE,marker=FALSE,legend=TRUE) # Add color scale for GC values
 > plot(D,ssu=TRUE) # Mark scaffolds containing SSU rRNA genes with crosshairs
 > plot(D,trna=TRUE) # Mark scaffolds containing tRNA genes with crosses
```

### 4b. Interactively choosing genomic bins

Identical syntax to `choosebin.genomestats`:

```R
 > Bin1 <- choosebin.diffcovstats(D)
```

### 4c. Fishing for connected contigs with Fastg data (experimental)

Identical syntax to `fastg_fishing.genomestats`:

```R
 > Bin2 <- fastg_fishing.genomestats(D,Bin1,"/path/to/file.fastg")
```

