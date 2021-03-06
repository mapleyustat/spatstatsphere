% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sphericaldomain.R
\name{sphere}
\alias{sphere}
\title{Title Generate a sphere}
\usage{
sphere(radius = 1, coord_type = c("geo_deg", "geo_rad", "polar"),
  unitname = NULL, ...)
}
\arguments{
\item{radius}{Numeric. Defaults to 1.}

\item{coord_type}{String. Type of coordinates used. See details.}

\item{unitname}{Optional. Name of unit of length. Either a single character
string, or a vector of two character strings giving the singular and plural
forms, respectively.}

\item{...}{Ignored.}
}
\value{
Object of class \code{"sphericaldomain"}
}
\description{
Title Generate a sphere
}
\details{
Currently three coordinate types/conventions are supported:
\describe{
  \item{\code{"geo_deg"}}{
    Geographical coordinates in degrees. Latitude from -90 (south pole) to 90
    (north pole). Longitude from -180 to 180 where 0 is the prime meridian
    (Greenwich).
  }
  \item{\code{"geo_rad"}}{
    Geographical coordinates in radians. Latitude from -\eqn{\pi/2}{pi/2}
    (south pole) to \eqn{\pi}{pi}/2 (north pole). Longitude from
    -\eqn{\pi/2}{pi/2} to \eqn{\pi/2}{pi/2} where 0 is the prime meridian
    (Greenwich).
  }
  \item{\code{"polar"}}{
    Polar coordinates in radians. Latitude from 0 (north pole) to
    \eqn{\pi}{pi} (south pole). Longitude from 0 to 2*\eqn{\pi}{pi} where 0
    is the prime meridian (Greenwich).
  }
}
}

