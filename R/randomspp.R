#' Generate uniformly distributed point pattern on sphere
#'
#' @param n Number of points to generate.
#' @template domain
#' @param nsim Number of simulated realisations to be generated.
#' @param drop Logical. If nsim=1 and drop=TRUE (the default), the result will be a spherical point pattern, rather than a list containing a spherical point pattern.
#' @template dots_sphere
#'
#' @return Spherical point pattern of class \code{"spp"}.
#' @export
#'
#' @examples
#' runifspp(10)
#'
runifspp <- function(n, domain = NULL, nsim = 1, ..., drop = TRUE){
  if(is.null(domain)){
    domain <- do.call(sphere, list(...))
  }
  if(domain$type!="sphere")
    stop("Only implemented for full sphere at the moment.")
  # Generate polar coords
  lat <- acos(2*runif(n)-1)
  long <- 2*pi*runif(n)
  # Transform to relevant type
  latlong <- transspherecoords(lat, long, from = "polar", to = domain$coord_type)
  return(spp(latlong$lat, latlong$long, domain = domain, check = FALSE))
}

#' Generate Poisson point pattern on sphere
#'
#' @param lambda Intensity. A positive numeric.
#' @template domain
#' @param nsim Number of simulated realisations to be generated.
#' @param drop Logical. If nsim=1 and drop=TRUE (the default), the result will be a spherical point pattern, rather than a list containing a spherical point pattern.
#' @template dots_sphere
#'
#' @return List of \code{nsim} spherical point patterns (of class \code{"spp"}).
#'  Or a spherical point pattern if nsim=1 and drop=TRUE (the default).
#' @export
#'
#' @examples
#' rpoisspp(10)
#'
rpoisspp <- function(lambda, domain = NULL, nsim = 1, ..., drop = TRUE){
  if(is.null(lambda) || !is.finite(lambda) || lambda<0)
    stop("Argument lambda should be a finite non-negative number.")
  n <- rpois(1, lambda)
  return(runifspp(n, domain = domain, nsim = nsim, ..., drop = drop))
}
