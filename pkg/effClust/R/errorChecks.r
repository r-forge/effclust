#' Common error checks
#' @noRd
errorChecks <- function(object, cluster, exclude, include.only, data=NULL) {

    if ( inherits(object, "formula") && is.null(data) )
         stop("Please provide data.")

    if (is.character(cluster) && length(cluster) == 1 &&
        !(cluster %in% names(data)))
        stop(paste(cluster, "not found in data."))

    if ( !is.null(exclude) && !is.null(include.only) )
       stop("Only one of exclude or include.only may be specified.")

    # Checks to make sure cluster specification works, when not a formula.
    if ( !inherits(cluster, "formula") ) {

        if ( inherits(object, "formula") && !length(cluster) %in% c(1, NROW(data)) )
            stop("length(cluster) must equal NROW(data) when cluster is a variable.")

        if ( !inherits(object, "formula") && length(cluster) != nobs(object) )
            stop("length(cluster) is different than number of observations.
                 \nIt might help to specify cluster as formula." )

        if ( inherits(object, "formula") && anyNA(cluster) )
            warning("NAs in cluster argument. Do these match missing
                    cases in the X matrix?")
    }

    if (inherits(cluster, "formula") & length(cluster) > 2) {
            stop("Multiway clustering is not allowed;
                 cluster formula should have only one term.")
        }

    return(NULL)
}
