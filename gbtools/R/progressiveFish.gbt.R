progressiveFish.gbt <- function(x, 
                            bin, 
                            depthLimit=5, 
                            fastg.file, 
                            paths.file, 
                            fasta.file, 
                            script.path="~/tools/my_scripts/genome-bin-tools/accessory_scripts/fastg_paths_fishing.pl") {
    # Initialize output variables
    binList <- vector("list", depthLimit+1)
    Fishing_depth <- 0:depthLimit
    Total_length <- vector("numeric",depthLimit+1)
    Num_scaffolds <- vector("numeric",depthLimit+1)
    Num_SSUs <- vector("numeric",depthLimit+1)
    marksource <- bin$summary$Marker_sources
    markList <- vector("list",length(marksource))
    # Summary values from the "seed" bin
    binList[[1]] = bin
    Total_length[1] <- bin$summary$Total_length
    Num_scaffolds[1] <- bin$summary$Num_scaffolds
    Num_SSUs[1] <- bin$summary$Num_SSUs
    # Marker summaries from "seed" bin
    for (i in 1:length(marksource)) {
        Fishing_depth <- 0:depthLimit
        Num_markers <- vector("numeric",depthLimit+1)
        Num_unique_markers <- vector("numeric",depthLimit+1)
        Num_singlecopy_markers <- vector("numeric",depthLimit+1)
        markList[[i]] <- data.frame(Fishing_depth=Fishing_depth,
                                   Num_markers=Num_markers,
                                   Num_unique_markers=Num_unique_markers,
                                   Num_singlecopy_markers=Num_singlecopy_markers)
        markList[[i]]$Fishing_depth[1] <- 0
        markList[[i]]$Num_markers[1] <- bin$summary$Num_markers[i]
        markList[[i]]$Num_unique_markers[1] <- bin$summary$Num_unique_markers[i]
        markList[[i]]$Num_singlecopy_markers[1] <- bin$summary$Num_singlecopy_markers[i]
    }
    # Iterate across depth values and get connected bins
    for (i in 1:depthLimit) {
        # Perform Fastg fishing
        binConnect <- fastgFish(x=x,
                                bin=bin,
                                fastg.file = fastg.file,
                                paths.file = paths.file,
                                fasta.file=fasta.file, 
                                depth = i,
                                script.path = script.path)
        # Add connected bin to list of output bins
        binList[[i+1]] = binConnect
        # Add summary statistics to vectors
        Total_length[i+1] = binConnect$summary$Total_length
        Num_scaffolds[i+1] = binConnect$summary$Num_scaffolds
        Num_SSUs[i+1] = binConnect$summary$Num_SSUs
        # Add marker statistics to lists
        for (j in 1:length(marksource)) {
            markList[[j]]$Num_markers[i+1] <- binConnect$summary$Num_markers[j]
            markList[[j]]$Num_unique_markers[i+1] <- binConnect$summary$Num_unique_markers[j]
            markList[[j]]$Num_singlecopy_markers[i+1] <- binConnect$summary$Num_singlecopy_markers[j]
        }
    }
    # Combine summary stats into data frames
    summary <- data.frame(Fishing_depth=Fishing_depth,
                         Total_length=Total_length,
                         Num_scaffolds=Num_scaffolds,
                         Num_SSUs=Num_SSUs)
    # Return output as list
    output <- list(binList=binList,summary=summary,markList=markList,marksource=marksource)
    class(output) <- "progressiveFish"
    return(output)
}