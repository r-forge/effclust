#' Convenience function for checking for single NA
#' Similar in spirit to isTRUE
#' @noRd 
isNA <- function(x) { length(x)==1 && is.na(x) }

#' Process the include.only and exclude args and pass back
#' a list of coefs to compute effective number of clusters for.
#' @noRd
incl.excl <- function(include.only, exclude, tags, fixed) {
    if (isNA(include.only) & isNA(exclude)) {
        return(tags)
    }
    if (!isNA(include.only)) {
        # Expand the regex for inclusion.
        include.vars <- unlist(sapply(include.only, grep, tags,
            fixed = fixed, value = TRUE
        ))
        # What should be removed?
        exclude.vars <- setdiff(tags, include.vars)
        exclude.vars <- exclude.vars[!is.na(exclude.vars)]
    }
    if (!isNA(exclude)) {
        # Expand regex for exclusion.
        exclude.vars <- unlist(sapply(exclude, grep, tags,
            fixed = fixed, value = TRUE
        ))
        exclude.vars <- exclude.vars[!is.na(exclude.vars)]
    }
    return(setdiff(tags, exclude.vars))
}

#' Compute Approximate Effective Number of Clusters for Regression Coefficients
#'
#' Specifically, for each coefficient the function returns the quantity
#' \eqn{G^{*A}}, the feasible version of \eqn{G^*}, introduced by
#' Carter, Schnepel, and Steigerwald (2017).  \eqn{G^{*A}} does not
#' accommodate multi-way clustering.
#'
#' @param object a formula or an object of class "lm", "plm", or
#' "fixest" (with method "feols")
#' @param cluster cluster identifier:
#' a variable, a formula (see details), or, a length one
#' character vector naming a variable in \code{data}.
#' @param data a dataframe (see details)
#' @param subset an optional vector used to select cases from \code{data}
#' @param exclude a vector of regular expressions for variables to be excluded
#' from the return
#' @param include.only a vector of regular expressions for variables to be
#' included in the return (implying others excluded)
#' @param fixed logical indicates how regular expressions in \code{exclude}
#' or \code{include.only} should be evaluated (as in \code{grep})
#' @param nominal logical indicating whether the number of clusters should
#' be included as the first element of the return vector
#' @param rho numeric scalar that changes assumed variance matrix (see details)
#' @return vector of effective number of clusters for coefficients
#' @details
#' When \code{object} is a formula, it does not need a response variable
#' (left-hand side), but if the response variable might have different
#' missing cases than the right-hand side, it should be included. Cluster
#' fixed effects, if any, must be explicitly included in
#' the formula, or \code{data} should be within-transformed data. (The
#' \pkg{data.table} package is especially handy for the coding the
#' transformation; see the example.)
#'
#' For regression objects, when \code{cluster} is a formula or a one
#' element character vector, it is evaluated in the context of the data
#' used for the model, and the \code{data} and \code{subset} arguments
#' are ignored with a warning.
#'
#' The \code{data} argument is required when \code{object} is a formula.
#' In this case, \code{object} is evaluated in the context of \code{data}.
#' If \code{cluster} is a formula or length 1 character vector, it is
#' interpreted as the name of a column of \code{data}.
#'
#' If \code{cluster} is a variable, it must not contain \code{NA} and its
#' length must equal \code{nobs(reg)} or the number of rows in (subsetted)
#' data when \code{object} is a formula.
#'
#' By default \eqn{G^{*A}} is returned for all coefficients.
#' The output may be limited by using \code{include.only}
#' or \code{exclude}. This does not affect the computation
#' for coefficients that are included, only the vector returned.
#'
#' The regular expressions used by \code{include.only} or \code{exclude}
#' should refer to the contents of \code{names(coef(reg))}, which might differ
#' from how the regression formula refers to the same things.
#' For example, groups of variables entered by using something
#' like \code{factor(statefip)} in the regression formula show up in the
#' coefficient vector names looking like \code{factor(statefip)27}.
#' They can be excluded en mass by the regular expression
#' \code{exclude="^factor\\(statefip\\)"} (parentheses need to be
#' escaped in a regular expression). Alternatively, you can
#' specify \code{fixed=TRUE}, which has the same meaning as for \code{grep}.
#'
#' Carter et al. recommend assuming every element of the variance matrix of
#' errors in a cluster is 1, a conservative scenario. However, as explained
#' in "A Note on Computation," the default \code{rho} is 0.999, so that the
#' variance matrix has 1 on the diagonal and 0.999 off the diagonal.  It is
#' possible to set \code{rho=1}, but that produces incorrect values in
#' when there are cluster fixed effects. In other cases the
#' computation is not very sensitive to the value of \code{rho}, so there
#' is no appreciable difference between \code{rho=1} and \code{rho=0.999}.
#' MacKinnon, Nielsen, and Webb (2022) argue that cluster fixed effects
#' usually remove much of the intra-cluster correlation, justifying \code{rho}
#' close to zero. (However, this is inconsistent with the ``worst case
#' scenario'' approach that motivates Carter et al.'s recommendation.)
#'
#' @examples
#' # some data with correlated errors
#' set.seed(85914270)
#' G <- 50
#' cl.sizes <- sample(10:100, G, replace=TRUE)
#' n <- sum(cl.sizes)
#' id <- rep(1:G, cl.sizes)
#' X1 <- rchisq(n, 5)
#' X2 <- ifelse(id %% 4 == 0, 1, 0)
#' e <- rnorm(n)
#' eg <- rep(rnorm(G), cl.sizes)
#' Y <- 1 + 2*X1 + 3*X2 - 1*X1*X2 + e + eg
#' d <- data.frame(Y, X1, X2, id)
#'
#' r <- lm(Y ~ X1 + X1:X2 + factor(id), data=d)
#'
#'effClust(r, ~d$id, exclude=c("factor\\(id\\)","Intercept"))
#'
#' library(data.table)
#' setDT(d)
#' d[ , `:=`(X1dot = X1 - mean(X1),
#'           X2dot = X2 - mean(X2),
#'           X1X2dot = X1*X2 - mean(X1*X2)),
#'   by = id]
#'
#' effClust(~ -1+X1dot+X1X2dot, cluster=~id, data=d)
#' @references Andrew V. Carter, Kevin T. Schnepel, and Douglas G. Steigerwald,
#' "Asymptotic Behavior of a $t$-test Robust to Cluster Heterogeneity,"
#' \emph{The Review of Economics and Statistics,} October 2017, 99(4).
#' \doi{10.1162/REST_a_00639}.
#'
#' James G. MacKinnon, Morten Ørregaard Nielsen, and Matthew D. Webb,
#' "Leverage, Influence, and the Jackknife in Clustered Regression Models:
#' Reliable Inference Using summclust," QED Working Paper 1483, Queen's
#' University (2022). \url{https://www.econ.queensu.ca/research/working-papers/1483}.
#' @export
effClust <- function(object, cluster, data=NA, subset=NA,
                     include.only=NA, exclude=NA,
                     fixed=FALSE, nominal=FALSE, rho=0.999) {
    UseMethod("effClust")
}


#' Do the actual calculation of G*.
#' @noRd
calcGstar <- function(X, XpXinv, cluster, tags, rho, nominal) {
    cl <- unique(cluster)
    G <- length(cl)
    gamma <- matrix(NA_real_,  nrow=G, ncol=length(tags))
        # coef's gammas are in a column.
    #########################################################################
    # LOOP OVER CLUSTERS
    # Each loop iteration calculates gamma_g for all coefs in tags.

    for (g in seq_along(cl)) {
        Xg <- X[cluster == cl[g], , drop=FALSE] # drop=FALSE in case ng==1
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

#' Common error checks
#' @noRd
errorChecks <- function(object, cluster, data, subset, exclude, include.only) {

    if ( inherits(object, "formula") && isNA(data) )
            stop("Please provide data.")

    if (is.character(cluster) && length(cluster) == 1 &&
        !(cluster %in% names(data)))
        stop(paste(cluster, "not found in data."))

    if ( !isNA(exclude) && !isNA(include.only) )
       stop("Only one of exclude or include.only may be specified.")
    if ( !inherits(object, "formula") & (!isNA(data) || !isNA(subset)) )
        warning("data and subset are ignored when object is not a formula.",
                call.=FALSE)

    if ( !inherits(cluster, "formula") && anyNA(cluster) )
        stop("Can't have NA values in cluster argument
        \nIt might help to specify cluster as a formula.")

    if ( !inherits(cluster, "formula") && length(cluster != nobs(object)) )
        stop("length(cluster) is different than number of observations.
              \nIt might help to specify cluster as formula." )

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
