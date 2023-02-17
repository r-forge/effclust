#' Common error checks
#' @noRd
errorChecks <- function(object, cluster, exclude, include.only, data=NULL) {

    if ( inherits(object, "formula") && is.null(data) )
         stop("Please provide data.")

    if ( inherits(object, "formula") && grepl("|", format(object), fixed=TRUE) )
         stop("Two-part formulas not supported when object is a formula.\n
               Partial out fixed effects or use fixest::feols() instead.")

    if (is.character(cluster) && length(cluster) == 1 &&
        !(cluster %in% names(data)))
        stop(paste(cluster, "not found in data."))

    if ( !is.null(exclude) && !is.null(include.only) )
       stop("Only one of exclude or include.only may be specified.")

    # Checks to make sure cluster specification works
    if (!inherits(cluster, "formula")) {

        if ( length(cluster) == 1 && !cluster %in% names(data) )
            stop(paste(cluster, "not found in data."))

        if ( inherits(object, "formula") && !length(cluster) %in% c(1, NROW(data)) )
            stop("length(cluster) must equal NROW(data).")

        if ( !inherits(object, "formula") && length(cluster) != nobs(object) )
            stop("length(cluster) is different than number of observations.
                 \nIt might help to specify cluster as formula." )

        if ( inherits(object, "formula") && anyNA(cluster) )
            warning("NAs in cluster argument. These much match missing
                    cases in the X matrix.")
    }

    if (inherits(cluster, "formula")) {
        tt <- stats::terms(cluster)
        tlabs <- attr(tt, "term.labels")
        if (length(tlabs) > 1) {
            stop("Multiway clustering is not allowed;
                 cluster formula should have only one term.")
        }
    }

    return(NULL)
}
