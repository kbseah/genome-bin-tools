#' @export
userAdd.gbt <- function(x,
                        userTab,
                        userSource=NA
                        ) {
    ## Check that userTab is data.frame with col "scaffold" ###################
    if (!is.data.frame(userTab) ||
        length(which(names(userTab)=="scaffold"))==0 ||
        is.na(userSource) ) {
        cat("gbtools ERROR: Please check inputs. See help(userAdd) \n")
    } else {
        ## Check that userTab scaffold IDs match x scaffold IDs ###############
        if (length(which(userTab$scaffold %in% x$scaff$ID))==0) {
            cat ("gbtools ERROR: Scaffold IDs in userTab don't match gbt object\n")
        } else {
            x$userTab[[length(x$userTab)+1]] <- userTab  # Append userTab
            x$userSource[length(x$userTab)] <- userSource # Append userSource
            # NB: Using c() will create discrepancy between userTab and userSource
            # because c() on an empty vector will create first element ""
            x$call[[length(x$call)+1]] <- match.call()  # Record function call
            return(x)  # Return result
        }
    }
}
