#' Distance on sphere.
#'
#' Low-level function to calculate distance on the sphere. No checking
#' of arguments is done.
#'
#' @param lat1,lat2 Numeric of length 1. Latitude in radians between -pi/2 and
#'   pi/2.
#' @param long1,long2 Numeric of length 1. Longitude in radians between -pi and
#'   pi (or between 0 and 2*pi).
#' @param domain Sphere object of class \code{"sphericaldomain"} specifying
#'   radius of sphere and type of coordinates. Defaults to \code{sphere(1)}.
#'
#' @return Numeric of length 1.
#' @export
spheredist <- function(lat1, long1, lat2, long2, domain = sphere(1)){
  ctype <- domain$coord_type
  if(ctype!="geo_rad"){
    latlong1 <- transspherecoords(lat1, long1, from = ctype, to = "geo_rad")
    latlong2 <- transspherecoords(lat2, long2, from = ctype, to = "geo_rad")
    lat1 <- latlong1$lat
    long1 <- latlong1$long
    lat2 <- latlong2$lat
    long2 <- latlong2$long
  }
  domain$radius * unitspheredist(lat1, long1, lat2, long2)
}

#' Distance on unitsphere.
#'
#' Low-level function to calculate distance on the sphere. No checking
#' of arguments is done.
#'
#' @param lat1,lat2 Numeric of length 1. Latitude in radians between -pi/2 and
#'   pi/2.
#' @param long1,long2 Numeric of length 1. Longitude in radians between -pi and
#'   pi (or between 0 and 2*pi).
#'
#' @return Numeric of length 1.
unitspheredist <- function(lat1, long1, lat2, long2){
  2 * asin(sqrt(sin((lat1 - lat2)/2)^2 + cos(lat1)*cos(lat2)*sin((long1 - long2)/2)^2))
}

#' Pairwise distances
#'
#' Computes the matrix of great circle (geodesic) distances between all pairs
#' of points in a point pattern on the sphere
#'
#' @param X object of class \code{"spp"}.
#' @param \dots ignored.
#' @return Symmetric square matrix where the (i,j)'th entry is the great circle
#' (geodesic) distance between points i and j.
#'
#' @export
#'
#' @examples
#' X <- spp(pi*c(1/4,1/2,3/4), 2*pi*c(1/4,1/2,3/4), domain = sphere(type = "polar"))
#' pairdist(X)
pairdist.spp <- function(X, ...){
  dom <- domain(X)
  n <- npoints(X)
  X <- as.data.frame(X, coord_type = "geo_rad")
  lat <- X$lat
  long <- X$long
  distfun <- function(i,j){
    unitspheredist(lat[i], long[i], lat[j], long[j])
  }
  rslt <- outer(1:n, 1:n, distfun)
  if(!all.equal(dom$radius, 1))
    rslt <- rslt * dom$radius
  return(rslt)
}

#' Nearest neighbour distances
#'
#' Computes the vector of great circle (geodesic) distances to nearest neighbour for
#' the points in a point pattern on the sphere
#'
#' @param X object of class \code{"spp"}.
#' @param \dots ignored.
#' @return Numeric of length n where n is the number of points in \code{X}.
#'
#' @export
#'
#' @examples
#' X <- spp(pi*c(1/4,1/2,3/4), 2*pi*c(1/4,1/2,3/4), domain = sphere(type = "polar"))
#' nndist(X)
nndist.spp <- function(X, ...){
  d <- pairdist(X)
  diag(d) <- Inf
  return(apply(d, 2, min))
}

#' Distance to border for points inside a spherical cap
#'
#' @param X Spherical point pattern of class \code{"spp"}
#'
#' @return Numeric
#' @export
bdist_spherepoints <- function(X){
  stopifnot(inherits(X, "spp"))
  dom <- domain(X)
  switch(dom$type,
         "sphere" = rep(Inf, npoints(X)),
         "cap" = {
           centre <- dom$cap$centre
           X <- coords(X)
           # Border dist equals total cap dist minus dist to cap centre.
           dom$cap$dist - spheredist(X$lat, X$long, centre[1], centre[2], dom)
         },
         stop("Unrecognized domain."))
}
