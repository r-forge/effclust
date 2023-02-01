#' @rdname effClust
#' @export
effClust.lm <- function(object, cluster, data=NA, subset=NA,
                     include.only=NA, exclude=NA,
                     fixed=FALSE, nominal=FALSE, rho=0.999) {

    errorChecks(object, cluster, data, subset, exclude, include.only)

    # One more check below after cluster is a variable for sure.

    if (is.character(cluster) & length(cluster)==1) {
        cluster <- stats::as.formula(paste0)("~", cluster)
    }

    ## First make cluster into a vector, if it's not already.
    if (inherits(cluster, "formula")) {
        tt <- stats::terms(cluster)
        cl.name <- attr(tt, "term.labels")
        # Is it like ~id or ~d$id?
        #no.dollar <- regexpr("$", cl.name, fixed=TRUE) == -1L
        #if (no.dollar) { # It's like ~id, so make it into ~data$id.
            #cl.name <- paste0("data$", cl.name)
            attributes(tt)$term.labels <- cl.name
            cluster <- stats::reformulate(attr(tt, "term.labels"))
        #}
        # Now integrate cluster into the model frame, coping with NAs.
        cl.exp <- stats::expand.model.frame(object, cluster,
                                    na.expand=TRUE)
        cluster <- cl.exp[ , cl.name]
    }

     if (length(cluster) != nobs(object)) {
        stop("length(cluster) is different than number of observations.
              \nIt might help to specify cluster as formula." )
    }

    # If cluster is character or factor, make it integer.
    # This only happens if the argument was a variable.
    if (is.factor(cluster)) cluster <- as.integer(cluster)
    if (is.character(cluster)) cluster <- as.integer(as.factor(cluster))

    all.tags <- names(coef(object))[!is.na(coef(object))]
    X <- model.matrix(object)[ , all.tags, drop=FALSE]
    XpXinv <- chol2inv(qr.R(qr(X)))
    rownames(XpXinv) <- all.tags
    colnames(XpXinv) <- all.tags

   tags <- incl.excl(include.only, exclude, tags=all.tags, fixed=fixed)

   calcGstar(X, XpXinv, cluster, tags, rho, nominal)
}