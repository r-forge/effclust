    #' @rdname effClust
#' @export
effClust.formula <- function(object, cluster, data, subset=NULL,
                     include.only=NULL, exclude=NULL,
                     fixed=FALSE, nominal=FALSE, rho=0.999, ...) {

    errorChecks(object, cluster, exclude, include.only, data)

    cl.name <- NA

# 1. What kind of thing is cluster?

    if (inherits(cluster, "formula")) {
        cl.name <- all.vars(cluster)
    }
    if (is.character(cluster) && length(cluster)==1) {
        cl.name <- cluster
    } # else cluster is an "outside" vector and cl.name stays NA.

# 2. If cluster is an "outside" vector, cbind it to data.
#    errorChecks() has already verified it will fit.
    if (is.character(cluster) && length(cluster) > 1) {
        data$cl4effClust <- cluster
        cl.name <- "cl4effClust"
    }

# 3. Apply all.vars to main formula.  Is cluster there?
#    Yes -> cluster fixed effects.
#    No -> we need to integrate it into data to make model.frame.
    reg.vars <- all.vars(object)  # This includes fixed effect vars after a |.
    if ( !cl.name %in% reg.vars ) {
      var.combo <- c(reg.vars, cl.name)
    } else {
      var.combo <- reg.vars
    }

# 4. Make a model.frame with all variables in the formula + cluster.
    f2 <- stats::as.formula( paste("~", paste(var.combo, collapse=" + ")) )
    mf <- str2lang("model.frame(f2, data=data[ , var.combo], subset=subset)")
    X <- eval(mf)

# 5. Separate cluster into separate vector and remove it from data.
    cluster <- X[ , cl.name]
    X <- X[ , reg.vars]

# 6. If there is a dependent variable in object, remove it from data.
#    Leave it in until now so that model.frame can deal with possible
#    NAs in it.  Also make a formula that doesn't have it for making
#    the model matrix.
#    If length(object) == 3, there is a LHS variable.
#    In that is true, the LHS is the second element.
    f3 <- object
    if (length(object) == 3) {
        dv <- all.vars(object[[2]])
        X[ , dv] <- NULL
        f3 <- stats::delete.response(stats::terms(object))
    }

# 7. If main formula has |, send the data off to be demeaned.
#    Otherwise, X is ready to be passed to model.matrix..
    if ( any(grepl("|", all.names(object), fixed=TRUE)) ) {
        fd <- f3[[2]]  # the whole RHS
        fd[[1]] <- as.name("~")   # Change | to ~.
        # Next line because fixest doesn't play well with data.table.
        if (any(class(data) == "data.table")) class(data) <- "data.frame"
        X <- fixest::demean(stats::as.formula(fd), data=data, as.matrix=TRUE)
    } else {
        X <- stats::model.matrix(f3, X)
    }

    all.tags <- colnames(X)

    tags <- incl.excl(include.only, exclude, tags=all.tags, fixed=fixed)

    effClust.default(X, cluster, tags, rho, nominal, XpXinv=NULL)

}
