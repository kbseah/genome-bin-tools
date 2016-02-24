#' Perform connectivity fishing with Fastg and paths files from SPAdes 3.6.2+
#'
#' Find contigs connected to a given gbtbin and report a new gbtbin object
#' Requires Fastg and paths files produced by SPAdes 3.6.2+; old version of
#' SPAdes produces different format.
#'
#' @param x Object of class gbt (parent object of the gbtbin object)
#' @param bin Object of class gbtbin, derived from x above
#' @param fastg.file Fastg formatted assembly graph from SPAdes
#' @param paths Paths file mapping assembly graph edge names to scaffold/contig names
#' @param save Logical: Save list of fished contigs to external file? (Default: No)
#' @param file File name to save list of fished contigs, if save=TRUE
#' @return Object of class gbtbin
#' @export

fastgFish.gbt <- function (x, # Object of class gbt (parent object of the bin)
                              bin, # Object of class gbtbin
                              fastg.file, # Fastg assembly graph
                              paths.file, # scaffolds.paths or contigs.paths file
                              save=FALSE, # Save result to external file? Default no
                              file="fished_bin.list" # File name to save result
                              ) {
    ## REPLACE THIS PATH WITH SCRIPT PATH ON YOUR LOCAL SYSTEM #############################################
    script.path <- "/home/kbseah/tools/my_scripts/genome-bin-tools/accessory_scripts/fastg_paths_fishing.pl"
    ########################################################################################################
    command <- "perl"
    command.params <- paste(script.path,
                            "-g", fastg.file,
                            "-p", paths.file,
                            "-o /tmp/tmp.fishing_output",
                            "-b -",
                            "-r")
    fished.contigs.list <- system2 (command,
                                    command.params,
                                    input=as.character(bin$scaff$ID),
                                    stderr=NULL,
                                    stdout=TRUE)
    newbin <- gbtbin (shortlist=fished.contigs.list,
                      x=x,
                      slice=NA,
                      taxon="Class",
                      save=save,
                      file=file)
    newbin$call[[length(newbin$call)+1]] <- match.call() # Record function call
    return(newbin)
}