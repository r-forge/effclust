#' @rdname effClust
#' @export
effClust.formula <- function(object, cluster, data, subset=NULL,
                     include.only=NULL, exclude=NULL,
                     fixed=FALSE, nominal=FALSE, rho=0.999, ...) {

    cl.name <- NA

    errorChecks(object, cluster, exclude, include.only, data)

    # One more check below after X matrix is obtained.

    # # Subset data if needed.
    # if (!is.null(subset)) {
    #     mf <- match.call()                     # copy of call to effClust
    #     mf[[1]] <- as.name("model.frame")      # Replace effClust with model.frame.
    #     names(mf)[1] <- "model.frame"          # Name that entry.
    #     mf <- mf[c("model.frame","data","subset")] # Keep only these args.
    #     data <- eval(mf, parent.frame())       # Evaluate call.
    #}


    ## Make cluster into a numeric vector, if it's not already.
    cl.name <- NA

    if (inherits(cluster, "formula")) {
        tt <- stats::terms(cluster)
        cl.name <- attr(tt, "term.labels")
    }

    if (is.character(cluster) && length(cluster)==1) {
        cl.name <- cluster
    }

    f2 <- object
    # If cl.name is still NA, then cluster is the full vector of
    # cluster identifiers, so we skip this block.
    if (!is.na(cl.name)) {
        # If cl.name not a term in object, tack it onto the formula.
        tt <- stats::terms(object)
        tt.lab <- labels(tt)
    
        if (!cl.name %in% tt.lab) { 
            # Is there a response variable in f2?
            # Otherwise with one-sided formula update gives back a formula
            # like .~X1
            yesY <- ifelse(length(f2) == 3, ".", "")
            f2 <- stats::update(f2, new=str2lang(paste0(yesY, "~.+", cl.name)))
        }
    }

    mf <- match.call()
    mf[[1]] <- as.name("model.frame")
    mf$object <- f2
    names(mf)[c(1,2)] <- c("model.frame","formula")
    mf <- mf[c("model.frame","formula","data","subset")]
    X <- eval(mf, parent.frame())

    if (!is.na(cl.name)) {
        # Otherwise, cluster is a full vector and we don't 
        # need to extract it from X or worry about removing it.
        cluster <- X[, cl.name]
        if (!cl.name %in% tt.lab) X[ , cl.name] <- NULL
    }

    X <- stats::model.matrix(object, X)
        
    if (length(cluster) != NROW(X))
        stop("length(cluster) is different than number of observations.
        \nIt might help to specify cluster as a formula.")

    all.tags <- colnames(X)
    # qrX <- qr(X)
    # if (qrX$rank < ncol(X)) stop("Formula includes linear dependencies.")
    # XpXinv <- chol2inv(qr.R(qrX))
    # dimnames(XpXinv) <- list(all.tags, all.tags)

   tags <- incl.excl(include.only, exclude, tags=all.tags, fixed=fixed)

   effClust.default(X, cluster, tags, rho, nominal, XpXinv=NULL)

}