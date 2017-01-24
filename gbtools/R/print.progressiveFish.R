print.progressiveFish <- function(x) {
    cat("*** Object of class progressiveFish ***\n")
    cat("Fishing depth: ")
    print(length(x$binList)[1])
    cat("Summary statistics:\n")
    print(x$summary)
    cat("Number of marker sets: ")
    print (length(x$markList)[1])
}