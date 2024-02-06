#' @rdname effClust
#' @export
effClust.fixest <- function(object, cluster,
                     include.only=NULL, exclude=NULL,
                     fixed=FALSE, nominal=FALSE, rho=0.999, ...) {


    # Error checks

    if (object$method != "feols")
        stop("fixest method must be feols.")

    errorChecks(object, cluster, exclude, include.only)

    # One more check below after X matrix is obtained.

    ## Make cluster into a numeric vector, if it's not already.
    cl.name <- NA

    if (inherits(cluster, "formula")) {
        tt <- stats::terms(cluster)
        cl.name <- attr(tt, "term.labels")
    }

    if (is.character(cluster) && length(cluster)==1) {
        cl.name <- cluster
    }

    tt <- stats::terms(object)

    # If cl.name is still NA, then cluster is the full vector of
    # cluster identifiers, so we skip this block.
    if (!is.na(cl.name)) {
        # Add the cluster identifier to object formula then
        # build a subsetted model frame with cluster in it.
        # Then extract cluster.
        # Can't just convert that to a model matrix because
        # it is the untransformed data, so X is constructed
        # separately below.
        # model.frame(object) would return the transformed data.
        tt.lab <- labels(tt)  # Does not include terms after | in formula.
        # Two cases:
        #   1. cl.name is in the "main" formula --> just go
        #   2. cl.name is not --> need to tack it on.
        f2 <- stats::formula(object, type="linear") # the "main" part
        if (!cl.name %in% tt.lab) {
            # Is there a response variable in f2?
            # Otherwise with a one-sided formula, update gives back a formula
            # like .~X1
            yesY <- ifelse(length(f2) == 3, ".", "")
            f2 <- stats::update(f2, new=str2lang(paste0(yesY, "~.+", cl.name)))
        }

        # Need to omit subset from the call to model.frame because with
        # fixest it is either a vector or a formula.  The vector will be
        # in parent.frame(), not here, and model.frame doesn't accept
        # a formula in subset. So, we  need to subset it ourselves.
        # ALSO (grrr), need to set na.action = na.pass. Otherwise mdata
        # can end up with different rows than s.set below.  Then use
        # na.omit() after s.set is integrated into mdata.
        mf <- object$call
        mf$na.action = "na.pass"
        mf[[1]] <- as.name("model.frame")
        mf$fml <- f2   # mf$fml is actually object$fml_all, not formula(object)
        names(mf)[c(1,2)] <- c("model.frame", "formula")
        mf <- mf[c("model.frame","formula","data","na.action")]
        mdata <- eval(mf, parent.frame())  # Not subsetted.

        ssCall <- object$call$subset
        if (!is.null(ssCall)) {
            if (inherits(eval(ssCall), "formula")) {
                mf <- object$call
                mf$na.action = "na.pass"
                mf[[1]] <- as.name("model.frame")
                mf$fml <- mf$subset
                names(mf)[c(1,2)] <- c("model.frame", "formula")
                mf <- mf[c("model.frame","formula","data","na.action")]
                s.set <- eval(mf, parent.frame())
                s.set <- s.set[ ,1] # need a vector, not a data.frame
            } else { # It's a vector in parent.frame()
                s.set <- eval(object$call$subset, parent.frame())
            }
        mdata <- mdata[s.set, ]
        }
        mdata <- na.omit(mdata)
        cluster <- mdata[ , cl.name]
    }

    # Unlike plm::model.matrix, fixest::model.matrix does not return the
    # transformed data, so need to de-mean it to get X.
    X <- fixest::demean(model.matrix(object, type="rhs"),
                        model.matrix(object, type="fixef"))
    # The following substitute is simpler, but generates a
    # data.table error if the data for object was a data.table.
    #X <- fixest::demean(object)

    all.tags <- colnames(X)

    if (anyNA(cluster))
        warning("cluster variable includes NA for some complete cases of regressors.")

   tags <- incl.excl(include.only, exclude, tags=all.tags, fixed=fixed)

   effClust.default(X, cluster, tags, rho, nominal, XpXinv=NULL)

}
