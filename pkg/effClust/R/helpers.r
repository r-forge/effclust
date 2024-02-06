#' Process the include.only and exclude args and pass back
#' a list of coefs to compute effective number of clusters for.
#' @noRd
incl.excl <- function(include.only, exclude, tags, fixed) {
    if (is.null(include.only) && is.null(exclude)) {
        return(tags)
    }
    if (!is.null(include.only)) {
        # Expand the regex for inclusion.
        include.vars <- unlist(sapply(include.only, grep, tags,
            fixed = fixed, value = TRUE
        ))
        # What should be removed?
        exclude.vars <- setdiff(tags, include.vars)
        exclude.vars <- exclude.vars[!is.na(exclude.vars)]
    }
    if (!is.null(exclude)) {
        # Expand regex for exclusion.
        exclude.vars <- unlist(sapply(exclude, grep, tags,
            fixed = fixed, value = TRUE
        ))
        exclude.vars <- exclude.vars[!is.na(exclude.vars)]
    }
    return(setdiff(tags, exclude.vars))
}
