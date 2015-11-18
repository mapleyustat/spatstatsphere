#' Title Generate a sphere
#'
#' @param radius Numeric. Defaults to 1.
#' @param coord_type String. Type of coordinates used. See details.
#' @param unitname Optional. Name of unit of length. Either a single character
#'   string, or a vector of two character strings giving the singular and plural
#'   forms, respectively.
#' @param ... Ignored.
#'
#' @details Currently three coordinate types/conventions are supported:
#' \describe{
#'   \item{\code{"geo_deg"}}{
#'     Geographical coordinates in degrees. Latitude from -90 (south pole) to 90
#'     (north pole). Longitude from -180 to 180 where 0 is the prime meridian
#'     (Greenwich).
#'   }
#'   \item{\code{"geo_rad"}}{
#'     Geographical coordinates in radians. Latitude from -\eqn{\pi/2}{pi/2}
#'     (south pole) to \eqn{\pi}{pi}/2 (north pole). Longitude from
#'     -\eqn{\pi/2}{pi/2} to \eqn{\pi/2}{pi/2} where 0 is the prime meridian
#'     (Greenwich).
#'   }
#'   \item{\code{"polar"}}{
#'     Polar coordinates in radians. Latitude from 0 (north pole) to
#'     \eqn{\pi}{pi} (south pole). Longitude from 0 to 2*\eqn{\pi}{pi} where 0
#'     is the prime meridian (Greenwich).
#'   }
#' }
#'
#' @return Object of class \code{"sphericaldomain"}
#' @export
sphere <- function(radius = 1, coord_type = c("geo_deg", "geo_rad", "polar"),
                   unitname = NULL, ...){
  coord_type <- match.arg(coord_type)
  unitname <- as.units(unitname)
  lat_range <- switch(coord_type,
                      geo_deg = c(-90, 90),
                      geo_rad = c(-pi/2, pi/2),
                      polar = c(0, pi))
  long_range <- switch(coord_type,
                       geo_deg = c(-180, 180),
                       geo_rad = c(-pi, pi),
                       polar = c(0, 2*pi))
  s <- list(type = "sphere", radius = radius,
            coord_type = coord_type, unitname = unitname,
            lat_range = lat_range, long_range = long_range)
  class(s) <- "sphericaldomain"
  return(s)
}

#' Print object of class sphericaldomain
#'
#' @param x Object of class \code{"sphericaldomain"}.
#' @param ... Ignored.
#'
#' @return NULL (invisibly)
#' @export
print.sphericaldomain <- function(x, ...){
  type <- x$type
  if(type=="cap"){
    msg <- paste0("Spherical cap containing points within distance ",
                  signif(x$cap$dist, 3), " of the centre (",
                  paste(signif(x$cap$centre,5), collapse = ","),
                  ") defined on:")
    spatstat:::splat(msg)
  }
  spatstat:::splat(paste("Sphere of radius", x$radius,
                         "using lat,long coordinates of type:",
                         sQuote(x$coord_type)))
  return(invisible(NULL))
}

#' Test whether points are inside spherical domain
#'
#' @param lat Latitude.
#' @param long Longitude.
#' @template domain
#'
#' @return Logical
#' @export
inside_sphericaldomain <- function(lat, long, domain){
  # Check args
  if(!(is.numeric(lat) & is.numeric(long)))
    stop("lat and long must be numeric vectors.")
  if(length(lat) != length(long))
    stop("lat and long differ in length!")

  # Check inside sphere
  latok <- lat >= domain$lat_range[1] & lat <= domain$lat_range[2]
  longok <- long >= domain$long_range[1] & long <= domain$long_range[2]
  ok <- latok & longok

  # Check inside cap if relevant
  if(domain$type=="cap"){
    centre <- domain$cap$centre
    ok[ok] <- spheredist(lat[ok], long[ok], centre[1], centre[2], domain = domain) < domain$cap$dist
  }

  return(ok)
}

#' Define spherical cap
#'
#' @param centre Numeric of length two with latitude and longitude of point on
#'   sphere to be used as centre for the cap.
#' @param dist Numeric of length one defining the extend of the cap. Should be
#'   positive and at most half the circumfrence of the sphere.
#' @param domain Sphere containing the cap (object of class
#'   \code{"sphericaldomain"}).
#' @param ... Arguments passed to \code{\link{sphere}} to specify \code{domain}
#'   if this is NULL (the default).
#' @param simplify Logical. Whether to simplify to full sphere if the cap
#'   extends over the full sphere.
#'
#' @return Object of class \code{"sphericaldomain"}.
#' @export
sphericalcap <- function(centre, dist, domain, ..., simplify = TRUE){
  # Resolve the containing sphere
  if(missing(domain) || is.null(domain)){
    domain <- do.call(sphere, list(...))
  }
  stopifnot(inherits(domain, "sphericaldomain"))
  # Resolve dist missing or NULL
  if(missing(dist) || is.null(dist)){
    # Return sphere if simplify is TRUE
    if(simplify)
      return(domain)
    # Otherwise set dist to maximal distance on sphere
    dist <- pi*domain$radius
  }
  stopifnot(is.numeric(dist) && length(dist)==1 && dist>0)
  # Resolve centre missing or NULL
  if(missing(centre) || is.null(centre)){
    centre <- c(0, 0)
  }
  centre <- spatstat::ensure2vector(centre)
  if(!inside_sphericaldomain(centre[1], centre[2], domain))
    stop("Centre not inside sphere.")
  # Check dist is at most circumfrence/2
  if(dist > pi*domain$radius)
    stop("dist is bigger than the maximal distance on the given sphere.")

  # Change domain to cap and add centre and dist
  domain$type <- "cap"
  domain$cap <- list(centre = centre, dist = dist)
  return(domain)
}

#' Convert object to a spherical domain with desired coordinate type
#'
#' @param x Object of class \code{"sphericaldomain"}.
#' @param coord_type String. Type of coordinates used for return coordinates.
#'   Defaults to coordinate type of \code{x}. See \code{\link{spp}} for
#'   supported coordinate types.
#'
#' @return Spherical domain of class \code{"sphericaldomain"}
#' @export
#'
as.sphericaldomain <- function(x, coord_type = NULL){
  stopifnot(inherits(x, "sphericaldomain"))
  coord_type <- ifelse(is.null(coord_type),
                       x$coord_type,
                       match.arg(coord_type, c("geo_deg", "geo_rad", "polar")))
  dom <- sphere(x$radius, x$coord_type, unitname = x$unitname)
  if(x$type=="cap"){
    dom$type <- "cap"
    centre <- transspherecoords(x$centre[1], x$centre[2],
                                from = x$coord_type, to = coord_type)
    dom$centre <- unlist(centre)
    dom$dist <- x$dist
  }
  return(dom)
}

#' Area of spherical domain
#'
#' @param w Objcet of class \code{"sphericaldomain"}.
#' @export
#'
#' @examples
#' dom <- sphere(2)
#' area(dom)
#' dom <- sphericalcap(c(0,0), pi/2, radius = 2)
#' area(dom)
area.sphericaldomain <- function(w){
  R <- w$radius
  a <- switch(w$type,
              cap = 2*pi*R^2*(1-cos(w$cap$dist/R)),
              sphere = 4*pi*R^2,
              stop("Unknown domain type."))
  return(a)
}
