#' @rdname effClust
#' @export
effClust.lm <- function(object, cluster,
                     include.only=NULL, exclude=NULL,
                     fixed=FALSE, nominal=FALSE, rho=0.999, ...) {

    errorChecks(object, cluster, exclude, include.only)
    # One more possible warning below after cluster is a variable for sure.

    if (is.character(cluster) & length(cluster)==1) {
        cluster <- stats::as.formula(paste0)("~", cluster)
    }

    ## First make cluster into a vector, if it's not already.
    if (inherits(cluster, "formula")) {
        cl.name <- labels(stats::terms(cluster))
        if (cl.name %in% labels(object)) { 
            cluster <- stats::model.frame(object)[ , cl.name]
        } else {
            # Integrate cluster into the model frame, which copes with NAs.
                # Note: expand.model.frame throws a "lengths differ" error if 
                # cluster is a formula like ~d2$id where d2 has a different 
                # number of rows than the data used for object.
            cl.exp <- stats::expand.model.frame(object, cl.name, na.expand=TRUE)
            cluster <- cl.exp[ , cl.name]
        }
    }

    if (anyNA(cluster)) 
        warning("cluster variable includes NA for some complete cases of regressors.")

    # # If cluster is character or factor, make it integer.
    # # This only happens if the argument was a variable.
    # if (is.factor(cluster)) cluster <- as.integer(cluster)
    # if (is.character(cluster)) cluster <- as.integer(as.factor(cluster))

    #all.tags <- names(coef(object))[!is.na(coef(object))]
    #X <- model.matrix(object)[ , all.tags, drop=FALSE]
    X <- stats::model.matrix(object)
    if (anyNA(coef(object))) {
        # A variable was dropped to avoid linear dependency.
        # Need to remove that variable from the model matrix.
        # It was already removed from the qr, so that can be 
        # used directly..
        b <- names(coef(object))[!is.na(coef(object))] 
        X <- X[ , b]
    }
    all.tags <- colnames(X)

   tags <- incl.excl(include.only, exclude, tags=all.tags, fixed=fixed)

   effClust.default(X, cluster, tags, rho, nominal, XpXinv=NULL)
}