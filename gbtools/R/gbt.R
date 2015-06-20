#' Create new gbt object
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
#' @seealso \code{\link{gbtbin}}, \code{\link{choosebin}}
#'
#' @export
gbt <- function(covstats,mark,marksource,ssu,trna) UseMethod ("gbt")
