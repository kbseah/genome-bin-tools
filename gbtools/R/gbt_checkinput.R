#' Check input files for gbtools
#'
#' Check input files for the \code{\link{gbt}} function before importing them.
#' This calls a Perl script, located in inst/Perl/input_validator.pl within the
#' R package source. The input validator script can also be called outside of
#' the R environment, e.g. if you wish to check the files before starting an R
#' session. Calling input_validator.pl without arguments will show usage help.
#'
#' See documentation online https://github.com/kbseah/genome-bin-tools for
#' fuller instructions on formats of the input files required.
#'
#' @param covstats File(s) with coverage statistics of a metagenomic assembly;
#'        output of pileup.sh in BBTools package (required). More  than one file
#'        can be specified with c().
#' @param mark Table of scaffolds with marker genes and taxonomic information
#'        for each marker. E.g. AMPHORA2 or Phyla-AMPHORA marker sets and output
#'        parsed by parse_phylotype_result.pl. More than one file can be
#'        specified with c() (optional)
#' @param ssu Table of scaffolds with SSU rRNA genes, and taxonomic info
#'        for each SSU. E.g. use get_ssu_for_genome_bin_tools.pl. (optional)
#' @param trna Table of tRNA genes found in assembly. Can use the output from
#'        tRNAscan-SE directly. (optional)
#' @param outdir Folder to write fixed versions of the input file, where this
#'        is possible.
#' @param log Path to write log file of the input validator script, which you
#'        can read for more details (Default: current working folder)
#' @param script.path Path to input validator script. Default is to look in the
#'        location where gbtools package is installed.
#'
#' @return data.frame of file names and errors found per file
#' @seealso \code{\link{gbt}} - Import data and create gbt object for analysis
#' @export

gbt_checkinput <- function(covstats,
                           mark=NULL,
                           ssu=NULL,
                           trna=NULL,
                           outdir=NULL,
                           log=NULL,
                           script.path=system.file("Perl","input_validator.pl",package="gbtools")
                           ) {
    command <- "perl"
    # Required parameters
    command.params <- paste(script.path,
                            "--rflag",
                            "--covstats", paste(covstats,collapse=","))
    # Optional parameters
    if (!is.null(mark)) {
        command.params <- paste(command.params,
                                "--mark", paste(mark,collapse=","))
    }
    if (!is.null(ssu)) {
        command.params <- paste(command.params,
                                "--ssu", ssu)
    }
    if (!is.null(trna)) {
        command.params <- paste(command.params,
                                "--trna", trna)
    }
    if (!is.null(outdir)) {
        command.params <- paste(command.params,
                                "--outdir", outdir)
    }
    result.vec <- system2 (command,
                           command.params,
                           stderr=NULL,
                           stdout=TRUE)
    result.df <- as.data.frame(matrix(result.vec, ncol=2, byrow=TRUE))
    names(result.df) <- c("file","errors")
    result.df$errors <- as.numeric(levels(result.df$errors))[result.df$errors]
    
    # Check for errors and return counts
    if (sum(result.df$errors)>0) {
        message ("Errors found in input files. Check validator log for details")
    } else {
        message ("No errors found in input files")
    }
    return (result.df)
}