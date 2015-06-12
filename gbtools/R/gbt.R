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
#' @seealso \code{\link{gbtbin}}, \code{link{choosebin}}
#'
#' @export
gbt <- function(covstats,mark,marksource,ssu,trna) UseMethod ("gbt")
gbt.default <- function (covstats,  # Vector of filenames for coverage tables 
                         mark=NA,   # Vector of filenames for marker gene taxonomy tables
                         marksource=NA,  # Vector of source names for each marker gene table
                         ssu=NA,    # Filename for SSU annotation table
                         trna=NA    # Filename for tRNA annotation table
                         ) {
## Create new gbt objects
    if ( class(covstats)!="character" || length(covstats)==0 ) {  # Check that covstats argument is character class
        cat ("gbtools ERROR: covstats argument must be a list of file names!\n")
    }
    else {
        ## Read Coverage tables ##############################################
        if (length(covstats)==1) {
            scaff <- read.table(file=as.character(covstats),
                                sep="\t",
                                header=T)
            covs <- data.frame(ID=scaff$ID,
                               scaff$Avg_fold)
        }
        else {
            scaff <- read.table(file=as.character(covstats[1]),
                                sep="\t",
                                header=T)  # Contains all other data associated per contig
            covs <- data.frame(ID=scaff$ID,
                               scaff$Avg_fold)  # Contains scaffold ID and coverage data
            for (i in 2:length(covstats)) {  # Read the other covstats files
                scafftemp <- read.table(file=as.character(covstats[i]),
                                        sep="\t",
                                        header=T)
                covs <- merge(covs,
                              data.frame(ID=scafftemp$ID,
                                         scafftemp$Avg_fold),
                              by="ID")
            }
        }
        
        ## Read taxonomic marker table ##########################################
        if ( !is.na(marksource[1]) ) { # Read list of marker sources
            markTab <- NA
            if ( !is.na(mark[1]) ) {  # Read marker table
                for (i in 1:length(marksource)) {
                    markTabTemp <- read.table(file=as.character(mark[i]),
                                              sep="\t",  # Tab separated input only
                                              fill=TRUE,  # Fill ragged lines 
                                              header=T)
                    namesvec <- names(markTabTemp)
                    if (length(which(!markTabTemp$scaffold %in% scaff$ID)) > 0) {
                        # Catch marker tables where scaffold IDs dont match
                        # what's in the scaff table
                        cat ("gbtools WARNING: Scaffold IDs in marker table ")
                        print (as.character(mark[i]))
                        cat (" doesn't match contig coverage tables \n")
                    }
                    numMarksTemp <- dim(markTabTemp)[1]
                    sourcevector <- rep(as.character(marksource[i]),numMarksTemp)
                    markTabTemp <- cbind(markTabTemp,sourcevector)
                    names(markTabTemp) <- c(namesvec,"source")
                    if ( length(markTab)==1 && is.na (markTab)) {
                        markTab <- markTabTemp
                    } else {
                        markTab <- rbind(markTab,markTabTemp)
                    }
                }
            } else {
                cat ("gbtools WARNING: Marker tables not supplied. marksource parameter ignored \n")
                marksource <- NA
                numMarks <- NA
                markTab <- NA
            }
        } else {
            cat ("gbtools WARNING: marksource not supplied.\n")
        }
        
        ## Read SSU marker table ##############################################
        if ( !is.na(ssu[1]) ) {  
            ssuTab <- read.table(file=as.character(ssu),sep="\t",header=T)
            numSsu <- dim(ssuTab)[1]
        } else {
            ssu <- NA
            ssuTab <- NA
            numSsu <- NA
        }
        
        ## Read tRNA marker table #############################################
        if ( !is.na(trna[1]) ) {  # Read tRNA marker table
            trnaTab <- read.table(file=as.character(trna),sep="\t",skip=3,header=F)
            names(trnaTab) <- c("scaffold",
                                "tRNA_no",
                                "tRNA_begin",
                                "tRNA_end",
                                "tRNA_type",
                                "Anticodon",
                                "Intron_begin",
                                "Intron_end",
                                "Cove_score")
            numTrna <- dim(trnaTab)[1]
        } else {
            trna <- NA
            trnaTab <- NA
            numTrna <- NA
        }
        
        ## Generate summary statistics #######################################
        summarystats <- list(Total_length=sum(scaff$Length),
                             Num_scaffolds=length(scaff$ID),
                             Scaff_length_max=max(scaff$Length),
                             Scaff_length_min=min(scaff$Length),
                             Scaff_length_median=median(scaff$Length),
                             Scaff_length_N50=getN50(scaff$Length),
                             Num_markers=table(markTab$source),
                             Num_SSU=numSsu,
                             Num_tRNAs=numTrna)
        
        ## Package and return result #########################################
        userTab <- list()  # Create userTab as an empty list to hold user-custom data
        userSource <- ""  # Create userSource vector to hold names of user-custom data
        result <- list(scaff=scaff,
                       covs=covs,
                       markTab=markTab,
                       markSource=marksource,
                       ssuTab=ssuTab,
                       trnaTab=trnaTab,
                       userTab=userTab,
                       userSource=userSource,
                       summary=summarystats)
        result$call <- list()
        result$call[[1]] <- match.call()  # Record function call that produces this gbt object
        class(result) <- "gbt"
        result
    }    
}