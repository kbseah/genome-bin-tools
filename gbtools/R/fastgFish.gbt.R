#' @export

fastgFish.gbt <- function (x, # Object of class gbt (parent object of the bin)
                           bin, # Object of class gbtbin
                           fastg.file, # Fastg assembly graph
                           paths.file, # scaffolds.paths or contigs.paths file
                           fasta.file, # scaffolds.fasta or contigs.fasta file
                           depth=0, # Number of iterations to fish (Default 0 - exhaustive)
                           save=FALSE, # Save result to external file? Default no
                           file="fished_bin.list", # File name to save result
                           script.path=system.file("Perl","fastg_paths_fishing.pl",package="gbtools")
                           ) {
    script.path=script.path
    command <- "perl"
    command.params <- paste(script.path,
                            "-g", fastg.file,
                            "-p", paths.file,
                            "-s", fasta.file,
                            "-o /tmp/tmp.fishing_output",
                            "-b -",
                            "-i", depth,
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