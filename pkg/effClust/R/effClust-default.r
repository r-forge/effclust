#' @rdname effClust
#' @export
effClust.default <- function(object, cluster, tags=colnames(object), rho=0.999,
                             nominal=FALSE,  XpXinv=NULL, ...) {

    if ( !is.matrix(object) ) object <- as.matrix(object)
    # object is really the X matrix!
    if (is.null(XpXinv)) {
        qrX <- qr(object)
        if (qrX$rank < ncol(object)) stop("object matrix includes linear dependencies.")
        XpXinv <- chol2inv(qr.R(qrX))
        dimnames(XpXinv) <- list(colnames(object), colnames(object))

    }

    if (!is.integer(cluster)) { # It is real, character, or factor
        
        # Deal with the unlikely possibility that cluster is real numbers with decimals.
        if ( is.numeric(cluster) && !all(cluster==as.integer(cluster)) )
            cluster <- as.integer(factor(cluster))

        if (is.character(cluster)) cluster <- factor(cluster)

        cluster <- as.integer(cluster)
    }
    clusters <- unique(cluster)
    G <- length(clusters)
    gamma <- matrix(NA_real_,  nrow=G, ncol=length(tags))
        # coef's gammas are in a column.
    #########################################################################
    # LOOP OVER CLUSTERS
    # Each loop iteration calculates gamma_g for all coefs in tags.
    for (g in 1:G) {
        Xg <- object[cluster == clusters[g], , drop=FALSE] # drop=FALSE in case ng==1
        ng <- NROW(Xg)

        # Now get to making the gammas.
        # Calculation follows equations (38) in MacKinnon, Nielsen, and Webb,
        # "Leverage, Influence, etc."
        # rather than direct calculation based on Carter et al., which prevents
        # problem of huge assumed variance matrix when ng is large.
        # Note: Using sapply() below calculates only the diagonal elements of the
        #       matrix multiplication in the original formula, which should be more
        #       efficient, but it makes using a "selection vector" with more than
        #       a single nonzero value impossible. The error checks and so forth
        #       for selection are in the code for the 0.5 version.
        XpXg <- crossprod(Xg)
        gam0 <- sapply(tags, function(v){XpXinv[v, ] %*% XpXg %*% XpXinv[ , v]})    # eqn (38a)
        if (rho > 0) {
            ones <- matrix(rep(1, ng), nrow=1)                               # 1 X ng
            U <- ones %*% Xg                                                 # 1 X k
            U <- t(U) %*% U                                                  # k X k
            gam1 <- sapply(tags,
                           function(v){XpXinv[v, ] %*% U %*% XpXinv[ , v]})  # eqn (38b)
        } else {
            gam1 <- 0
        }
        gamma[g, ] <- rho*gam1 + (1-rho)*gam0                                # eqn (38c)
    }

    #########################################################################
    # MAIN FORMULA
    Gamma.numer <- ((G-1)/G) * apply(gamma, MARGIN=2, FUN=var)
    Gamma.denom <- apply(gamma, MARGIN=2, FUN=mean)^2
    Gamma <- Gamma.numer/Gamma.denom
    Gstar <- G/(1 + Gamma)
    names(Gstar) <- tags
    if (nominal) Gstar <- c(nominal=G, Gstar)
    return(Gstar)
}
