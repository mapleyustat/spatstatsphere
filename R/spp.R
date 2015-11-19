#' Import spatstat
#' @import spatstat


#' @title Spherical point pattern
#'
#' @param lat Numeric. Latitudes of points.
#' @param long Numeric. Longitude of points.
#' @param domain Domain where the points occur. Object of class
#'   \code{"sphericaldomain"}.
#' @template dots_sphere
#' @param check Logical. Check that points actually are inside the provided
#'   \code{domain}.
#'
#' @return Object of class \code{"spp"}.
#' @export
#'
#' @examples
#' dom <- sphere(coord_type = "polar")
#' lat <- pi * c(1/4, 1/2, 3/4)
#' long <- 2 * pi * c(1/4, 1/2, 3/4)
#' X <- spp(lat, long, domain = dom)
#' X <- spp(lat, long, coord_type = "polar")
#'
#' lat <- c(45, 0, -45)
#' long <- c(90, 180, -90)
#' X <- spp(lat, long)
#' dom <- sphere(coord_type = "geo_deg")
#' X <- spp(lat, long, domain = dom)
#'
spp <- function(lat, long, domain = NULL, ..., check = TRUE){
  if(is.null(domain)){
    domain <- do.call(sphere, list(...))
  }
  n <- (length(lat) + length(long))/2
  if (check && n > 0) {
    ok <- inside_sphericaldomain(lat, long, domain)
    nout <- sum(!ok)
    if (nout > 0) {
      warning(paste(nout, ngettext(nout, "point was", "points were"),
                    "rejected as lying outside the specified window"))
      rejects <- data.frame(lat = lat[!ok], long = long[!ok])
      lat <- lat[ok]
      long <- long[ok]
      n <- length(lat)
    }
  }
  else nout <- 0

  rslt <- ppx(data.frame(lat=lat, long=long), domain = domain, coord.type = rep("s", 2))
  class(rslt) <- c("spp", class(rslt))

  if (check && anyDuplicated(rslt))
    warning("data contain duplicated points")
  if (nout > 0)
    attr(rslt, "rejects") <- rejects
  return(rslt)
}

#' Print object of class spp
#'
#' @param x Object of class \code{"spp"}.
#' @param ... Ignored.
#' @return NULL (invisibly)
#' @export
print.spp <- function(x, ...){
  spatstat:::splat(paste("Spherical point pattern with", npoints(x), "points."))
  spatstat:::splat("Domain:")
  print(domain(x))
  # Temporary hack: Detect attribute set by simulation algorithm for DPPs
  nmean <- attr(x, "nmean")
  if(!is.null(nmean)){
    cat(paste("It has been simulated from a model where\n",
              "the expected number of points is: ",
              signif(nmean,4), "\n", sep=""))
  }
  return(invisible(NULL))
}

#' Extract domain of spherical point pattern
#'
#' @param X Spherical point pattern of class \code{"spp"}.
#' @param ... Ignored.
#'
#' @return Object of class \code{"sphericaldomain"}.
#' @export
#'
#' @examples
#' lat <- c(45, 0, -45)
#' long <- c(90, 180, -90)
#' X <- spp(lat, long)
#' domain(X)
domain.spp <- function(X, ...){
  X$domain
}

#' Extract number of points in a spherical point pattern.
#'
#' @param x Spherical point patten of class \code{"spp"}.
#'
#' @export
#'
#' @examples
#' lat <- c(45, 0, -45)
#' long <- c(90, 180, -90)
#' X <- spp(lat, long)
#' npoints(X)
#'
npoints.spp <- function(x){
  npoints.ppx(x)
}

#' Extract Coordinates of a Point Pattern on a Sphere
#'
#' @param x Object of class \code{"spp"}
#' @param ... Ignored
#'
#' @return a \code{data.frame} with one row for each point, containing the
#'   coordinates.
#' @export
#'
coords.spp <- function(x, ...){
  coords.ppx(x, temporal = FALSE, local = FALSE)
}

#' Check for duplicated points in a spherical point pattern.
#'
#' @param x Spherical point patten of class \code{"spp"}.
#' @param ... Additional arguments passed to
#'   \code{\link{anyDuplicated.data.frame}}.
#'
#' @return An integer or real vector of length one with value the 1-based index
#'   of the first duplicate if any, otherwise 0.
#' @export
#'
anyDuplicated.spp <- function(x, ...){
  anyDuplicated(as.data.frame(x), ...)
}

#' Extract coordinates of a spherical point pattern into a data.frame
#'
#' @param x Spherical point patten of class \code{"spp"}.
#' @param row.names NULL or a character vector giving the row names for the data
#'   frame. Missing values are not allowed.
#' @param ... Ignored.
#' @param coord_type String. Type of coordinates used for return coordinates.
#'   Defaults to coordinate type of \code{x}. See \code{\link{spp}} for
#'   supported coordinate types.
#'
#' @return \code{data.frame} with coordinates.
#' @export
#'
as.data.frame.spp <- function(x, row.names = NULL, ..., coord_type = NULL){
  current <- domain(x)$coord_type
  coord_type <- ifelse(is.null(coord_type),
                       current,
                       match.arg(coord_type, c("geo_deg", "geo_rad", "polar")))
  x <- coords(x)
  if(current!=coord_type){
    x <- transspherecoords(x$lat, x$long, from = current, to = coord_type)
  }
  data.frame(lat = x$lat, long = x$long, row.names = row.names)
}

#' Convert Data To Class spp
#'
#' Tries to coerce any reasonable kind of data to a spherical point pattern (an
#' object of class \code{"spp"}) for use by the \pkg{spherespatstat} package).
#'
#' Converts the dataset \code{X} to a point pattern (an object of class
#' \code{"spp"}).
#'
#' This function is normally used to convert an existing point pattern dataset,
#' stored in another format, to the \code{"spp"} format.  To create a new point
#' pattern from raw data such as lat,long coordinates, it is normally easier to
#' use the creator function \code{\link{spp}}.
#'
#' The dataset \code{X} may be:
#'
#' \itemize{
#'   \item an object of class \code{"spp"}
#'   \item a matrix or data frame with at least two columns
#'   \item a structure with entries \code{lat}, \code{long} which are numeric
#'     vectors of equal length
#'   \item a numeric vector of length 2, interpreted as the coordinates of a
#'     single point.
#' }
#'
#' The first case is typically used to change (or ensure) a specific coordinate
#' format by specifying either the argument \code{domain} directly or indirectly
#' using the argument \code{coord_type} which is passed to \code{\link{sphere}}
#' through the additional arguments \dots. In the last three cases the default
#' behaviour is to assume the domain is the entire sphere with coordinates of
#' type \code{"geo_deg"} (see \code{\link{spp}}). Alternatively the domain can
#' be specified by the argument \code{domain} or through the additional
#' arguments \dots.
#'
#' If \code{X} is a matrix or data frame, the first and second columns will be
#' interpreted as the \eqn{lat} and \eqn{long} coordinates respectively. Any
#' additional columns will be ignored.
#'
#' The function \code{as.spp} is generic, with methods for the classes
#' \code{"spp"}, \code{"matrix"}, \code{"data.frame"} and a default method.
#'
#' Point pattern datasets can also be created by the function \code{\link{spp}}.
#'
#' @aliases as.spp as.spp.spp as.spp.matrix as.spp.data.frame as.spp.default
#' @param X Data which will be converted into a point pattern
#' @template dots_sphere
#' @return An object of class \code{"spp"} (see \code{\link{spp}}) describing
#'   the spherical point pattern and its domain.
#' @keywords spatial manip
#' @export
#' @examples
#'
#' dom <- sphere(coord_type = "polar")
#' lat <- pi * c(1/4, 1/2, 3/4)
#' long <- 2 * pi * c(1/4, 1/2, 3/4)
#' X <- spp(lat, long, dom)
#' Y <- as.spp(X, coord_type = "geo_deg")
#' as.data.frame(Y)
as.spp <- function(X, ...){
  UseMethod("as.spp")
}

#' @describeIn as.spp Method for spp objects.
#' @template domain
#' @export
as.spp.spp <- function(X, ..., domain = NULL){
  if(is.null(domain)){
    sphere_args <- resolve.defaults(list(...), unclass(domain(X)))
    domain <- do.call(sphere, sphere_args)
  }
  co <- as.data.frame(X, coord_type = domain$coord_type)
  spp(co$lat, co$long, domain = domain, check = FALSE)
}

#' @describeIn as.spp Method for data frames.
#' @export
as.spp.data.frame <- function(X, ..., domain = NULL){
  if(is.null(domain))
    domain <- do.call(sphere, list(...))
  spp(X[[1]], X[[2]], domain = domain, check = TRUE)
}

#' @describeIn as.spp Method for matrices.
#' @export
as.spp.matrix <- function(X, ..., domain = NULL){
  if(is.null(domain))
    domain <- do.call(sphere, list(...))
  spp(X[,1], X[,2], domain = domain, check = TRUE)
}

#' @describeIn as.spp Default method.
#' @export
as.spp.default <- function(X, ..., domain = NULL){
  if(is.null(domain))
    domain <- do.call(sphere, list(...))
  lat <- X$lat
  long <- X$long
  if(is.null(lat) || is.null(long)){
    if(is.vector(X) && length(X) == 2){
      lat <- X[1]
      long <- X[2]
    }
  }
  spp(lat, long, domain = domain, check = TRUE)
}

transspherecoords <- function(lat, long, from, to){
  stopifnot(all(c(from, to) %in% c("geo_deg", "geo_rad", "polar")))
  # Convert input to geo_rad if necessary
  if(from=="geo_deg"){
    lat <- lat/180*pi
    long <- long/180*pi
  }
  if(from=="polar"){
    lat <- pi/2-lat
    long[long>pi] <- long[long>pi] - 2*pi
  }
  # Convert from deg_rad to relevant type
  if(to=="geo_deg"){
    lat <- lat*180/pi
    long <- long*180/pi
  }
  if(to=="polar"){
    lat <- 90-lat
    long[long<0] <- 2*pi + long[long<0]
  }
  return(list(lat = lat, long = long))
}

#' Extract a subset of a point pattern on a sphere. Extraction of a subset has
#' the effect of thinning the points.
#'
#' @param x object of class \code{"spp"}.
#' @param i Subset index. A valid subset index in the usual R sense, indicating
#'   which points should be retained.
#' @param j Ignored. (Required for compatibility with the generic function.)
#' @param drop Ignored. (Required for compatibility with the generic function.)
#' @param \dots Ignored. (Required for compatibility with the generic function.)
#' @return object of class \code{"spp"}.
#' @export
#' @examples
#'
#' X <- spp(c(-45,0,45), c(-10,0,160))
#' X[1:2]
#'
"[.spp" <- function(x, i, j, drop, ...) {
  verifyclass(x, "spp")
  if(missing(i) || npoints(x) == 0)
    return(x)
  dom <- domain(x)
  x <- coords(x)
  return(spp(x$lat[i], x$long[i], domain = dom))
}
