# gbtools: Interactive tools for metagenome visualization and binning in R

Refer to the [`gbtools` manual](https://github.com/kbseah/genome-bin-tools/wiki) for detailed instructions.

## Quick start

Assuming you already have metagenomes assembled, coverage calculated by mapping, marker genes identified and taxonomically classified. What next?

 * Download the latest release of gbtools
 * Install package dependencies in R (`sp` and `plyr` packages)
 * Install `gbtools` package in R
 * Format input data, and check with `input_validator.pl` (in the `accessory_scripts` folder)
 * Load your data into R: `A <- gbt(covstats=c("coverage_file1.tab","coverage_file2.tab"),mark="taxonomic_markers.tab",marksource="source",ssu="ssu_markers.tab",trna="trna_markers.tab")`
 * Plot GC-coverage plots: `plot(A, slice=1)`
 * Select contigs by coverage/gc/length cutoffs: `A.bin1 <- winnow(A,gc=c(0.25,0.35),len=c(10000,Inf),covmin=100,covmax=Inf,slice=1)`
 * Overlay bin on plot: `points(A.bin1,slice=1)`
 * Plot differential-coverage plot: `plot(A,slice=c(1,2))`
 * Overlay bin on differential coverage plot: `points(A.bin1, slice=c(1,2))`
 * Choose bin by interactively selecting a cloud of points on plot: `A.bin2 <- choosebin(A, slice=c(1,2))`
 * View summary statistics for a bin: `summary(A.bin2)`

## Getting help

Problems with using `gbtools`? Create a new issue using the GitHub issue-tracker on the right. Or send me an email, with "gbtools help" in the subject line.

Problems with input file formats? Read [the appendix](Appendix.\ Input\ file\ formats) and use the `input_validator.pl` script to check your input files.

## Citations 

Cite: Brandon Seah (2015), genome-bin-tools, Online: https://github.com/kbseah/genome-bin-tools

or via Zenodo: [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.18530.svg)](http://dx.doi.org/10.5281/zenodo.18530)

Cite dependencies if you use them:
* **R** -  R Core Team. 2014. R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. (http://www.R-project.org/)
* **BBMap** - Bushnell B. 2015. BBMap (http://sourceforge.net/projects/bbmap/)
* **AMPHORA2** - Wu M, Scott AJ. 2012. Bioinformatics 28 (7) : 1033-1034.
* **barrnap** - Seemann T. 2014. barrnap (http://www.vicbioinformatics.com/software.barrnap.shtml)
* **Usearch** - Edgar RC 2010. Bioinformatics 26 (19) : 2460-2461.
* **ARB-SILVA** - Quast C et al. 2013. Nucleic Acids Research 41 (D1) : D590-D596.
* **tRNAscan-SE** - Lowe T, Eddy S. 1997. Nucleic Acids Research 25 : 955-964.
* **Blobology** - Kumar S et al. 2013. Frontiers in Genetics 4 : 237

Contact: Brandon Seah (`kbseah@mpi-bremen.de`)
