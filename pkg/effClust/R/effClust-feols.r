#' @rdname effClust
#' @export
effClust.fixest <- function(object, cluster, data=NA, subset=NA,
                     include.only=NA, exclude=NA,
                     fixed=FALSE, nominal=FALSE, rho=0.999) {


    # Error checks

    if (object$method != "feols")
        stop("fixest method must be feols.")

    errorChecks(object, cluster, data, subset, exclude, include.only)

    # One more check below after X matrix is obtained.

    ## Make cluster into a numeric vector, if it's not already.
    cl.name <- NA

    if (inherits(cluster, "formula")) {
        tt <- stats::terms(cluster)
        cl.name <- attr(tt, "term.labels")
    }

    if (is.character(cluster) & length(cluster)==1) {
        cl.name <- cluster
    }

    # If cl.name is still NA, then cluster is the full vector of
    # cluster identifiers, so we skip this block.
    if (!is.na(cl.name)) {
        # Add the cluster identifier to object formula then
        # build a subsetted model frame with cluster in it.
        # Then extract cluster.
        # Can't just convert that to a model matrix because
        # it is the untransformed data, so X is constructed
        # separately below.
        tt <- stats::terms(object)
        f2 <- stats::reformulate(c(attr(tt, "term.labels"), cl.name),
                                 intercept=0)
        # Need to recover the original model frame first so that the
        # cluster identifier can be added to it.
        # A bunch of language-level rigamarole here because "fixest"
        # objects don't have a proper model.frame function.
        mf <- object$call
        mf[[1]] <- as.name("model.frame")
        mf$fml <- f2   # mf$fml was actually object$fml_all
        names(mf)[c(1,2)] <- c("model.frame", "formula")
        mf <- mf[c("model.frame","formula","data","subset")]
        mdata <- eval(mf, parent.frame())
        cluster <- mdata[ , cl.name]
    }  # Otherwise, cluster is a vector, assumed to match.

    # If cluster is character or factor, make it integer.
    # This only happens if the argument was a variable.
    if (is.factor(cluster)) cluster <- as.integer(cluster)
    if (is.character(cluster)) cluster <- as.integer(as.factor(cluster))

    all.tags <- names(coef(object)) # feols doesn't include dropped variables.

    # Unlike plm::model.matrix, fixest::model.matrix does not return the
    # transformed data, so need to de-mean it to get X.
    X <- fixest::demean(model.matrix(object, type="rhs"),
                        model.matrix(object, type="fixef"))
    # The following substitute is simpler, but generates a 
    # data.table error if the data for object was a data.table.
    #X <- fixest::demean(object)

    if (length(cluster) != nobs(object)) {
        stop("length(cluster) is different than number of observations.
              \nIt might help to specify cluster as formula.")
    }

    XpXinv <- chol2inv(qr.R(qr(X)))
    rownames(XpXinv) <- names(coef(object))
    colnames(XpXinv) <- rownames(XpXinv)

   tags <- incl.excl(include.only, exclude, tags=all.tags, fixed=fixed)

   calcGstar(X, XpXinv, cluster, tags, rho, nominal)

}