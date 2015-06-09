# Interactive metagenome binning tools

Contact: Brandon Seah (kbseah@mpi-bremen.de)

Cite: Brandon Seah (2015), genome-bin-tools, Online: https://github.com/kbseah/genome-bin-tools
or via Zenodo: [![DOI](https://zenodo.org/badge/10602/kbseah/genome-bin-tools.svg)](http://dx.doi.org/10.5281/zenodo.15812)

Cite dependencies if you use them:
* **R** -  R Core Team. 2014. R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. (http://www.R-project.org/)
* **BBMap** - Bushnell B. 2015. BBMap (http://sourceforge.net/projects/bbmap/)
* **AMPHORA2** - Wu M, Scott AJ. 2012. Bioinformatics 28 (7) : 1033-1034.
* **barrnap** - Seemann T. 2014. barrnap (http://www.vicbioinformatics.com/software.barrnap.shtml)
* **Usearch** - Edgar RC 2010. Bioinformatics 26 (19) : 2460-2461.
* **ARB-SILVA** - Quast C et al. 2013. Nucleic Acids Research 41 (D1) : D590-D596.
* **tRNAscan-SE** - Lowe T, Eddy S. 1997. Nucleic Acids Research 25 : 955-964.
* **Blobology** - Kumar S et al. 2013. Frontiers in Genetics 4 : 237

## 0. Introduction

Various tools and approaches exist for metagenomic binning - the process of defining individual genomes in a metagenomic assembly. These tools are designed for interactive exploration and binning of low-diversity microbial metagenomes in R.

A useful way to visualize a metagenomic assembly is to plot the coverage (depth) and GC% of the assembled scaffolds. Scaffolds coming from the same genome would tend to have similar coverage and GC%, and so form clusters in the plots. To aid in distinguishing the clusters, the taxonomic affiliation of each scaffold can be evaluated either by searching the entire scaffold sequence against a database like NCBI nr, or by searching specific marker genes. 

Examples of tools that use GC-coverage plots and taxonomic annotation include [Blobology] (https://github.com/blaxterlab/blobology) and [Metawatt] (http://sourceforge.net/projects/metawatt/).

Another visualization or binning method relies on the variation in coverage for different genomes between different samples. If the coverage of a metagenome assembly in one sample is plotted against the coverage in another sample, individual genomes would again tend to cluster together.

Examples of tools that use differential coverage binning: [Multi-metagenome] (http://madsalbertsen.github.io/multi-metagenome/), [GroopM] (http://minillinim.github.io/GroopM/).

Genome-bin-tools builds on concepts from Multi-metagenome, but it offers more:
 * Higher-level functions for plotting - Save time spent on typing and copy-pasting commands
 * Designed to work with free assembly and annotation tools (BBMap for mapping, barrnap for finding rRNAs, tRNAscan-SE for finding tRNAs, AMPHORA2 for finding marker genes)
 * Needing minimal software installation - start R, import some tables, load R functions, and go!
 * Interactive - select bins, see summary statistics for bins, save scaffold lists for later processing

## 1. Produce and annotate metagenomic assembly

If you want to follow along, you can use the data from Albertsen et al. 2013, which is in the folder `example_data`. The file names in the following examples use the file names of the example data. (NB: The original header names of the assembly Fasta file were edited to be shorter.)

### 1a. Assemble the metagenome and calculate coverage

Use your favorite assemblier, like [IDBA-UD] (http://i.cs.hku.hk/~alse/hkubrg/projects/idba/) or [SPAdes] (http://bioinf.spbau.ru/spades), to assemble your metagenome. The assembly should be in a Fasta file. Calculate coverage with [bbmap.sh] (http://sourceforge.net/projects/bbmap/) by mapping the original reads used to assemble the metagenome back onto the assembly:

```
 $ bbmap.sh ref=HPminus_assembly.fasta nodisk in=HPminus_reads.fq.gz covstats=HPminus.coverage
```

If you already have a SAM file with the mapping, you can use `pileup.sh` (also from the BBmap suite of tools) directly:

```
 $ pileup.sh in=mapping.sam out=coverage.table
```

BBmap is convenient because it calculates length and GC along with the coverage per scaffold, and outputs them all to the same file. If you use a different mapping tool, you will have to aggregate all these data together in a tab-separated table with the same header names as the bbmap output.

NOTE: Newer versions of BBmap (after 33.x) insert a comment character `#` in the header line of the covstats file. Remove this character before running the scripts.

For differential coverage binning, you will need a second read library from a different sample where at least some of the genomes present in the first sample should also be present. Map the second read library onto the same assembly and save the coverage values as a different file:

```
 $ bbmap.sh ref=HPminus_assembly.fasta nodisk in=HPplus_reads.fq.gz covstats=HPplus.coverage
```

### 1b. Identify marker genes and find phylogenetic affiliation (optional)

Use [AMPHORA2] (https://github.com/martinwu/AMPHORA2) or [Phyla-AMPHORA] (https://github.com/martinwu/Phyla_AMPHORA) to identify conserved marker genes in your assembly, and to assign a taxonomic position. Parse the output of their script Phylotyping.pl (here called `phylotype.result`) for import into R:

```
 $ perl parse_phylotype_result.pl -p phylotype.result > phylotype.result.parsed
```

This generates a file called `phylotype.result.parsed` which will be imported into R.

An alternative is to Blast contigs directly against a database like NCBI nt and use the NCBI taxon IDs to assign a taxon to each contig. This is the approach used by [Blobology] (https://github.com/blaxterlab/blobology); one of their scripts has been modified to produce an output table compatible with `genome-bin-tools`.
They require the NCBI nt Blast database and NCBI taxonomy dump, both of which are available from the NCBI FTP site (see the Blobology manual for instructions). Requires Blast+.

```
 $ perl blob_annotate_mod.pl --assembly HPminus_assembly.fasta --blastdb /path/to/ncbi/nt --out contig_taxonomy.tab --num_threads 8 --taxdump /path/to/ncbi/taxdump/
```

The output file `contig_taxonomy.tab` can be substituted in the subsequent workflow for `phylotype.result.parsed`. 

### 1c. Identify rRNA genes (optional)

Use [barrnap] (http://www.vicbioinformatics.com/software.barrnap.shtml) to detect SSU rRNA genes in the assembly, and assign phylotype by extacting sequences using [fastaFrombed] (http://bedtools.readthedocs.org/en/latest/) and then using [Usearch] (http://www.drive5.com/usearch/) against a curated [SILVA] (www.arb-silva.de/) database. The database has to be prepared in a specific way (instructions to come) but is identical to the Usearch-indexed database used by [phyloFlash] (https://github.com/HRGV/phyloFlash). PhyloFlash is also a great tool, why not check it out? (Disclosure: I helped to develop phyloFlash).

The rRNA extraction and output parsing is done with a wrapper script:
```
 $ perl get_ssu_for_genome_bin_tools.pl -d <path/to/ssu/database> -c <number_CPUs> -a HPminus_assembly.fasta -o <output_prefix> 
```

This generates a file called <output_prefix>.ssu.tab (in the example data, HPminus.ssu.tab) which will be imported into R.

### 1d. Identify tRNA genes (optional)

Use [tRNAscan-SE version 1.23] (http://selab.janelia.org/tRNAscan-SE/) to find tRNA genes. 

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

Install the gbtools package in R:

```R
 > install.packages("/PATH/TO/gbtools_2.0.tar.gz",repos=NULL,type="source")
```

where `/PATH/TO/` is replaced with the path to wherever you have the R source package.

This is recommended because you can call `help()` or `?` to read the documentation for each function within the R environment.

Alternatively, load the R functions with `source` (recommended if you want to tweak them or use experimental features):

```R
 > source("gbtools.r")
```

## 3. Import your data into R

`gbtools` treats each assembly and the associated coverage data as a single object of class `gbt`. If you have several coverage tables for the same assembly, you can import them together with the taxon marker, SSU rRNA, and tRNA annotations with the `gbt` function. The examples in this README can be replicated with the files in the `example_data` folder, adapted from the [Multi-metagenome] (http://madsalbertsen.github.io/multi-metagenome/) data set.

```R
 # If you have a single coverage table, give the filename at the covstats= option
 > d <- gbt(covstats="HPminus.coverage",mark="phylotype.result.parsed",ssu="HPminus.ssu.tab",trna="HPminus.trna.tab")
 # With two or more coverage tables, you can import them together with the c() function
 > d <- gbt(covstats=c("HPminus.coverage","HPplus.coverage"),mark="phylotype.result.parsed",ssu="HPminus.ssu.tab",trna="HPminus.trna.tab")
```
Only the coverage table is required. The rest are optional (though having them will provide more information and prettier plots).

Type the name of the object to see summary statistics:

```R
 > d
 > summary(d) # Does the same thing
```

## 4. Explore GC-coverage plots

### 4a. Basic plotting

GC-coverage plots are generated from single samples (i.e. coverage statistics from mapping a single read library onto a single assembly). Remember that you could import coverage data from several samples with the `gbt()` function into a single `gbt` object. You can specify which sample you want to plot using the `slice=` option:

```R
 > plot(d) # basic plot, defaults to plotting GC-coverage plot with first sample's coverage
 >         # if marker.list was imported, then colored automatically by marker taxonomy
 > plot(d,slice=2) # Plot the GC-coverage plot of the second sample instead.
 > plot(d,cutoff=2000) # Do not show scaffolds shorter than 2000 bp
 > plot(d,taxon="Phylum") # Color marker genes at the taxonomic level of "Phylum"
 > plot(d,ssu=TRUE) # Overlay crosshairs marking scaffolds that contain SSU rRNA genes
 > plot(d,ssu=TRUE,textlabel=TRUE) # Add text labels beside crosshairs showing phylotype assigned to SSU rRNA genes
 > plot(d,trna=TRUE) # Add cross marks showing scaffolds that contain tRNA rRNA genes
 > plot(d,legend=TRUE) # Add a legend showing which color corresponds to which taxonomic group; grey points are singleton taxa
```

Zoom into specific areas of the plot by altering the `xlim` and `ylim` parameters, as with the basic `plot` function in R.

### 4b. Interactively choosing genomic bins

If you see a cluster of scaffolds which you would like to save as a bin, you can choose it interactively by picking the points that draw a polygon surrounding the scaffolds you want. You must specify the sample that was used to draw the plot that you are now interacting with, using the `slice=` option. Specifying the wrong sample will give meaningless results.

```R
 > bin1 <- choosebin(d, slice=1) # basic function
 > bin1 <- choosebin(d, slice=1, save=TRUE, file="bin1.scaffolds.list") # Save the names of scaffolds in this bin to an external file called bin1.scaffolds.list
 > bin1 <- choosebin(d, slice=1, num.points=10) # Change the number of corners of the polygon
```

`bin` is now an object of class `gbtbin`. Type the name of the bin object to see a summary of the bin. If you imported the marker, SSU, and/or tRNA data, a summary of how many of each are contained in the bin will be reported (this is useful if the marker genes are typical single-copy genes, for example).

### 4c. Visualize bin for different samples

Let's say you've defined a bin based on the GC-coverage plot of sample 1. You want to see if they are still clustering together in sample 2. All you have to do is to draw a new plot and draw in the points of your bin, but change the `slice=` parameter to show sample 2 instead.

```R
 > plot(d,slice=2) # Plot the GC-coverage plot for the assembly using coverage data from sample 2
 > points(bin1, slice=2) # Overlay with your bin, using coverage data from sample 2
```

### 4d. Taking subsets of genomic bins

If you want to get a subset of scaffolds based on GC, coverage, or length, use the function `winnow()`

```R
 > bin2 <- winnow(d,slice=1, covmin=200, covmax=200) # Return the subset of contigs in d which have coverage above 200 in sample 1
 > bin2 <- winnow(d,gc=c(0.2,0.5),len=c(1000,Inf)) # Return the subset of scaffolds in d that have GC between 20-50%, and length above 1000
 > bin2 <- winnow(d,slice=c(1,2),covmin=c(200,300),covmax=c(Inf,5000)) # Subset of contigs that have coverage above 200 in sample 1, and between 300-5000 in sample 2.
```

If you want to get a subset of scaffolds containing marker genes that belong to a particular taxon, use the function `winnowMark()`

```R
 > bin3 <- winnowMark(d,param="Class",value="Gammaproteobacteria") # Returns all scaffolds containing a marker gene whose value for "Class" is "Gammaproteobacteria
```
The functions `winnow()` and `winnowMark()` work for both `gbt` and `gbtbin` objects

### 4e. Adding and subtracting bins

You can also perform set operations on bins. The `add` function takes the union, while the `lej` function takes the set difference. `lej` is non-commutative. I.e. `lej(bin1,bin2)` is not equivalent to `lej(bin2,bin1)`. It returns the members of the first object that are not in the second object.

```R
 > mergedbin <- add(bin1, bin2)  # Returns bin with scaffolds present in both bins
 > bin1not2 <- lej(bin1, bin2)   # Returns scaffolds in bin1 that are not in bin2
```
These functions work only on bin objects (class `gbtbin`).

### 4f. Fishing for connected contigs using Fastg files (experimental)

Fastg files are generated by newer versions of the SPAdes assembler, and contain contig connectivity information generated during the assembly process. These can be useful, e.g. to "fish" scaffolds that are from the same genome but which were inadvertently left out of the interactively chosen bin.

This needs a genomestatsbin object and produces another genomestatsbin object.

For now you will have to manually edit the `genome_bin_tools.r` file to specify the location of the script `fastg_parser.pl` (included with this package) and also a directory to store temporary files (default: `/tmp/`). 

Fish out a new bin from an old bin using a Fastg file:
```R
 > bin4 <- fastgFishing(d,bin1,"path/to/fastg/file.fastg")
 > bin4 <- fastgFishing(d,bin1,"path/to/fastg/file.fastg", save=TRUE,file="bin2.scaffolds.list") # Save list of scaffolds in new bin to external file
```

You can compare the two bins by plotting them overlaid:
```R
 > plot(d,marker=FALSE,slice=1) # The underlying GC-coverage plot, with coloring turned off
 > points(bin2,col="blue",slice=1) # The new bin in blue
 > points(bin1,col="black",slice=1) # The original bin in black
```

## 5. Explore differential coverage plots

Differential coverage plots are generated from two separate coverage files (in this example: `HPminus.coverage` and `HPplus.coverage`). You've already seen how to import them into a `gbt` object earlier in step 3. 

### 5a. Plotting

To make differential coverage plots instead of GC-coverage plots, simply specify a pair of samples to the `slice=` parameter. For example, to have sample 1 coverage on the x-axis and sample 2 coverage on the y-axis, use the parameter `slice=c(1,2)`. Differential coverage plots are automatically generated with log scales for both axes, but you can override this by specifying your own option to the `log=` parameter.

Like with GC-coverage plots, you can overlay taxonomic marker genes, SSU rRNA markers, and tRNA markers, as well as add legend for the color scheme. In addition, you can choose to color the points by contig GC% instead of marker taxonomy (but you can't do both).

```R
 > plot(d,slice=c(1,2)) # Basic plot. Defaults to coloring by marker genes, if data imported
 > plot(d,slice=c(1,2),marker=FALSE,gc=FALSE) # Uncolored plot
 > plot(d,slice=c(1,2),marker=TRUE,legend=TRUE) # Add legend
 > plot(d,slice=c(1,2),gc=TRUE,marker=FALSE) # Color by GC
 > plot(d,slice=c(1,2),gc=TRUE,marker=FALSE,legend=TRUE) # Add color scale for GC values
 > plot(d,slice=c(1,2),ssu=TRUE) # Mark scaffolds containing SSU rRNA genes with crosshairs
 > plot(d,slice=c(1,2),trna=TRUE) # Mark scaffolds containing tRNA genes with crosses
```

### 5b. Interactively choosing genomic bins

The same `choosebin` function is used to choose a bin from a differential coverage plot, except that the `slice=` parameter must be the same as the one used to draw the plot. Adding and subtracting bins works in the same way, as is fishing with Fastg information. 

## 6. Mapping and reassembly

Once you have your final bin, either the shortlist of scaffolds in that bin was exported with the save parameter when that bin was defined, or you can use the native R function write():

```R
 > write(as.character(bin1$scaff$ID),file="shortlist.file") 
```

This gives you a list of scaffold names in your bin. Use a tool like faSomeRecords to retrieve the corresponding Fasta sequence records. Use a mapper like BBmap to get reads that map to those scaffolds, and then reassemble with your favorite assembler.

```
 $ faSomeRecords original_assembly.fasta shortlist.file shortlist.fasta 
 $ bbmap.sh ref=shortlist.fasta outputunmapped=f outm=reads_for_reassembly.fastq.gz in=original_reads.fastq
```

Good luck!
