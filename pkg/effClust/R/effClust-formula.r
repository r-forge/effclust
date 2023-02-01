#' @rdname effClust
#' @export
effClust.formula <- function(object, cluster, data=NA, subset=NA,
                     include.only=NA, exclude=NA,
                     fixed=FALSE, nominal=FALSE, rho=0.999) {

    cl.name <- NA

    errorChecks(object, cluster, data, subset, exclude, include.only)

    # One more check below after X matrix is obtained.

    # Subset data if needed.
    if (!isNA(subset)) {
        mf <- match.call()                     # copy of call to effClust
        mf[[1]] <- as.name("model.frame")      # Replace effClust with model.frame.
        names(mf)[1] <- "model.frame"          # Name that entry.
        mf <- mf[c("model.frame","data","subset")] # Keep only these args.
        data <- eval(mf, parent.frame())       # Evaluate call.
    }


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
        # Add the cluster identifier to object formula so that
        # model.matrix will do its magic about NAs.
        # Then remove it from X.
        tt <- stats::terms(object)
        f2 <- stats::reformulate(c(attr(tt, "term.labels"), cl.name),
                            intercept=attr(tt, "intercept"))
        X <- model.matrix(f2, stats::model.frame(f2, data))
        cluster <- X[ , cl.name]
        X <- X[ , -which(colnames(X) == cl.name), drop=FALSE]
    } else { # cluster is a vector, assumed to match properly with X.
        X <- model.matrix(object, stats::model.frame(object, data))
    }

    if (is.factor(cluster)) cluster <- as.integer(cluster)
    if (is.character(cluster)) cluster <- as.integer(as.factor(cluster))

    if (length(cluster) != NROW(X))
        stop("length(cluster) is different than number of observations.
        \nIt might help to specify cluster as a formula.")


    all.tags <- colnames(X)
    qrX <- qr(X)
    if (qrX$rank < ncol(X)) stop("Formula includes linear dependencies.")
    XpXinv <- chol2inv(qr.R(qrX))
    rownames(XpXinv) <- colnames(X)
    colnames(XpXinv) <- colnames(X)

   tags <- incl.excl(include.only, exclude, tags=all.tags, fixed=fixed)

   calcGstar(X, XpXinv, cluster, tags, rho, nominal)

}