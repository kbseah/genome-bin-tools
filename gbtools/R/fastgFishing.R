#' Perform connectivity fishing with Fastg file
#'
#' Find contigs connected to existing gbtbin object, using connectivity from
#' an external Fastg file
#'
#' @param x Object of class gbt, from which the bin object was defined
#' @param bin Object of class gbtbin, defined from the gbt object x
#' @param fastg.file Path to Fastg file for the metagenome assembly of x
#' @param save Logical - Save list of contigs of new bin to external file?
#' @param file File to save contig list 
#' @return Object of class gbtbin
#' @export
fastgFishing <- function(x, bin, fastg.file, ... ) UseMethod ("fastgFishing") 
fastgFishing.gbtbin <- function(x,  # Object of class gbt (parent object of the gbtbin)
                                bin,  # Object of class gbtbin defined from x
                                fastg.file,  # Fastg file for assembly of x
                                taxon="Class",  # Deprecated - user pls ignore
                                save=FALSE,  # Save list of contigs to external file?
                                file="fished_bin.list"  # File to save contig list
                                ) {
    command <- "perl"
## REPLACE THIS PATH WITH YOUR OWN PATH !! #########################################
    script.path <- "/home/kbseah/tools/my_scripts/genome-bin-tools/fastg_parser.pl" 
####################################################################################
    command.params <- paste(script.path,"-i",fastg.file,"-o /tmp/tmp.fishing_output -b - -r")
        # By default throws away fastg_parser.pl output to /tmp/
    fished.contigs.list <- system2(command,
                                   command.params,
                                   input=as.character(bin$scaff$ID),
                                   stderr=NULL,
                                   stdout=TRUE)
    newbin <- gbtbin(shortlist=fished.contigs.list,x=x,slice=NA,taxon=taxon,save=save,file=file)
    return(newbin)
}

