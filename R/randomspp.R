#' Generate uniformly distributed point pattern on sphere
#'
#' @param n Number of points to generate. Either a single positve integer or a
#'   vector of integers of length \code{nsim} indicating the number of points in
#'   each of the \code{nsim} generated patterns.
#' @template domain
#' @param nsim Number of simulated realisations to be generated.
#' @param drop Logical. If nsim=1 and drop=TRUE (the default), the result will
#'   be a spherical point pattern, rather than a list containing a spherical
#'   point pattern.
#' @template dots_sphere
#'
#' @return Spherical point pattern of class \code{"spp"}.
#' @export
#'
#' @examples
#' runifspp(10)
#'
runifspp <- function(n, domain = NULL, nsim = 1, ..., drop = TRUE){
  # Resolve n and nsim
  if(length(n) > 1){
    if(nsim != length(n))
      stop("Mismatch between nsim and length of first argument n.")
  } else{
    n <- rep(n, nsim)
  }
  if(is.null(domain)){
    domain <- do.call(sphere, list(...))
  }
  if(domain$type!="sphere")
    stop("Only implemented for full sphere at the moment.")
  result <- vector(mode = "list", length = nsim)
  for(i in 1:nsim){
    # Generate polar coords
    lat <- acos(2*runif(n[i])-1)
    long <- 2*pi*runif(n[i])
    # Transform to relevant type
    latlong <- transspherecoords(lat, long, from = "polar", to = domain$coord_type)
    result[[i]] <- spp(latlong$lat, latlong$long, domain = domain, check = FALSE)
  }
  if(nsim == 1 && drop)
    return(result[[1]])
  names(result) <- paste("Simulation", 1:nsim)
  return(result)
}

#' Generate Poisson point pattern on sphere
#'
#' @param lambda Intensity. Either a single positve numeric or a vector of
#'   length \code{nsim} indicating the intensity in each of the \code{nsim}
#'   generated patterns.
#' @template domain
#' @param nsim Number of simulated realisations to be generated.
#' @param drop Logical. If nsim=1 and drop=TRUE (the default), the result will
#'   be a spherical point pattern, rather than a list containing a spherical
#'   point pattern.
#' @template dots_sphere
#'
#' @return List of \code{nsim} spherical point patterns (of class \code{"spp"}).
#'   Or a spherical point pattern if nsim=1 and drop=TRUE (the default).
#' @export
#'
#' @examples
#' rpoisspp(10)
#'
rpoisspp <- function(lambda, domain = NULL, nsim = 1, ..., drop = TRUE){
  if(is.null(domain)){
    domain <- do.call(sphere, list(...))
  }
  if(is.null(lambda) || !is.finite(lambda) || lambda<0)
    stop("The intensity should be finite and non-negative.")
  n <- rpois(nsim, lambda*area(domain))
  return(runifspp(n, domain = domain, nsim = nsim, drop = drop))
}
