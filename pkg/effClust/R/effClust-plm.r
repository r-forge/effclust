# Setup work for "plm" is identical to "lm" objects.
#' @rdname effClust
#' @export
effClust.plm <- function(object, cluster, data=NA, subset=NA,
                     include.only=NA, exclude=NA,
                     fixed=FALSE, nominal=FALSE, rho=0.999) {

    return(effClust.lm(object, cluster, data))

}