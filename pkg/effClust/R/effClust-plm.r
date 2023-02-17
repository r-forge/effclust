# Setup work for "plm" is identical to "lm" objects.
#' @rdname effClust
#' @export
effClust.plm <- function(object, cluster,
                     include.only=NULL, exclude=NULL,
                     fixed=FALSE, nominal=FALSE, rho=0.999, ...) {
    return(effClust.lm(object, cluster, include.only, exclude,
                       fixed, nominal, rho))

}