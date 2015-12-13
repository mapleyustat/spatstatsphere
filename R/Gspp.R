#' G-function for spherical point pattern
#'
#' @template spp_X
#' @param \dots Ignored.
#' @param r Optional. Vector of values for the argument \eqn{r} at which
#'   \eqn{G(r)} should be evaluated. Users are advised \emph{not} to specify
#'   this argument; there is a sensible default. If necessary, specify
#'   \code{rmax}.
#' @param rmax Optional. Maximum desired value of the argument \eqn{r}.
#' @param breaks This argument is for internal use only.
#' @param correction Optional. A character vector containing any selection of
#'   the options \code{"none"} and \code{"border"} (and possibly others in the
#'   future).
#' @param lambda Intensity (usually not supplied). Estimated from data by default.
#' @param unnormalized Logical. Only for internal use at the moment.
#' @param angular_arg Logical. Interpret r values as angles (possibly in degrees
#'   rather than radians depending on coordinate type of \code{X}).
#' @return An object of class \code{"fv"}, see \code{\link{fv.object}}, which
#'   can be plotted directly using \code{\link{plot.fv}}.
#'
#'   Essentially a data frame containing columns \item{r}{the vector of values
#'   of the argument \eqn{r} at which the function \eqn{G} has been estimated }
#'   \item{theo}{the theoretical value \eqn{G(r) = 2 \pi (1-cos(r))}{G(r) = 2 *
#'   pi * (1 - cos(r))} for a stationary Poisson process (where r is in
#'   radians)} together with columns named \code{"none"}, and/or
#'   \code{"border"}, according to the selected edge corrections. These columns
#'   contain estimates of the function \eqn{G(r)} obtained by the edge
#'   corrections named.
#' @export
#'
#' @examples
#' X <- runifspp(100)
#' G <- Gspp(X)
#'
Gspp <- function(X, r = NULL, rmax = NULL, breaks = NULL, correction = NULL,
                 lambda = NULL, unnormalized = FALSE, angular_arg = FALSE){
  ratio <- FALSE # This could potentially be an argument of the function in the long run.

  # Check whether r values should be degrees (radians will be used and the chaged to degrees at the end)
  dom <- domain(X)
  degrees <- angular_arg && grepl("deg", dom$coord_type)
  if(!is.null(r) && degrees) r <- r*2*pi/360
  if(!is.null(rmax) && degrees) rmax <- rmax*2*pi/360

  # Figure out r-values
  rmaxdefault <- rmax %orifnull% pi
  breaks <- handle.r.b.args(r, breaks, dom, pixeps = rmaxdefault/128, rmaxdefault=rmaxdefault)
  r <- breaks$r
  rmax <- breaks$max
  zeroes <- numeric(length(r))

  # Choose correction based on domain if unset by user
  best_correction <- ifelse(dom$type=="sphere", "none", "border")
  if(missing(correction) || is.null(correction)){
    correction <- best_correction
  }
  correction <- pickoption("correction", correction,
                           c(none="none", border="border", best = best_correction),
                           multi=TRUE)

  # Intensity
  areaW <- area(dom)
  n <- npoints(X)
  if(is.null(lambda))
    lambda <- n/areaW

  ## initialise output object
  theo <- 1 - exp( -lambda * 2*pi * (1-cos(r)) )
  if(unnormalized) theo <- lambda * theo
  Odf <- data.frame(r = r, theo = theo)
  desc <- c(paste("distance argument", ifelse(!angular_arg, "r", "theta")),
            "theoretical uncorrected %s")
  ylab <- if(!angular_arg) quote(G(r)) else quote(G(theta))
  OO <- fv(Odf,
           argu = "r",
           ylab = ylab,
           valu = "theo",
           fmla = . ~ r,
           alim = c(0, ifelse(!degrees, rmax, rmax*360/(2*pi))),
           labl = c(ifelse(!angular_arg, "r", "theta"),
                    ifelse(!angular_arg, "%s[pois](r)", "%s[pois](theta)")),
           desc = desc,
           fname = "G",
           yexp = ylab
  )

  nnd <- nndist(X)

  if ("none" %in% correction) {
    ## Uncorrected for data on entire sphere.
    if (n <= 1)
      edf <- zeroes
    else {
      hh <- hist(nnd[nnd <= rmax], breaks = breaks$val,
                 plot = FALSE)$counts
      edf <- cumsum(hh)/length(nnd)
    }
    labl <- ifelse(!angular_arg, "hat(%s)(r)", "hat(%s)(theta)")
    OO <- bind.fv(OO, data.frame(raw = edf),
                  labl,
                  "uncorrected estimate of %s",
                  "raw")
  }

  if(degrees) OO$r <- OO$r*360/(2*pi)

  return(OO)
}
