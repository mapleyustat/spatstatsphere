#' K-function for spherical point pattern
#'
#' @template spp_X
#' @param \dots Ignored.
#' @param r Optional. Vector of values for the argument \eqn{r} at which
#'   \eqn{K(r)} should be evaluated. Users are advised \emph{not} to specify
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
#'   of the argument \eqn{r} at which the function \eqn{K} has been estimated }
#'   \item{theo}{the theoretical value \eqn{K(r) = 2 \pi (1-cos(r))}{K(r) = 2 *
#'   pi * (1 - cos(r))} for a stationary Poisson process (where r is in
#'   radians)} together with columns named \code{"none"}, and/or
#'   \code{"border"}, according to the selected edge corrections. These columns
#'   contain estimates of the function \eqn{K(r)} obtained by the edge
#'   corrections named.
#' @export
#'
#' @examples
#' X <- runifspp(100)
#' K <- Kspp(X)
#'
Kspp <- function(X, r = NULL, rmax = NULL, breaks = NULL, correction = NULL,
                 lambda = NULL, unnormalized = FALSE, angular_arg = FALSE){
  ratio <- FALSE # This could potentially be an argument of the function in the long run.

  # Check whether r values should be degrees (radians will be used and the chaged to degrees at the end)
  dom <- domain(X)
  degrees <- angular_arg && grepl("deg", dom$coord_type)
  if(!is.null(rmax) && degrees) rmax <- rmax*2*pi/360

  # Figure out r-values
  rmaxdefault <- rmax %orifnull% pi
  breaks <- handle.r.b.args(r, breaks, dom, pixeps = rmaxdefault/128, rmaxdefault=rmaxdefault)
  r <- breaks$r
  rmax <- breaks$max

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
  if(is.null(lambda)){
    lambda <- n/areaW
    lambda2 <- lambda * (n-1)/areaW
  } else{
    lambda2 <- lambda^2
  }

  ## initialise output object
  theo <- 2*pi*(1-cos(r))
  if(unnormalized) theo <- lambda2 * theo
  Odf <- data.frame(r = r, theo = theo)
  desc <- c(paste("distance argument", ifelse(!angular_arg, "r", "theta")),
            "theoretical uncorrected %s")
  ylab <- if(!angular_arg) quote(K(r)) else quote(K(theta))
  OO <- fv(Odf,
           argu = "r",
           ylab = ylab,
           valu = "theo",
           fmla = . ~ r,
           alim = c(0, ifelse(!degrees, rmax, rmax*360/(2*pi))),
           labl = c(ifelse(!angular_arg, "r", "theta"),
                    ifelse(!angular_arg, "%s[pois](r)", "%s[pois](theta)")),
           desc = desc,
           fname = "K",
           yexp = ylab
  )

  delta <- pairdist(X)
  diag(delta) <- Inf

  if(any(correction == "none")) {
    ## Uncorrected for data on entire sphere.
    hh <- hist(delta, breaks = c(r, Inf), plot = FALSE)$counts
    khat <- c(0, cumsum(hh[-length(hh)]))
    if(!unnormalized) khat <- khat/(lambda2*areaW)

    ## uncorrected estimate
    labl <- ifelse(!angular_arg, "hat(%s)(r)", "hat(%s)(theta)")
    OO <- bind.fv(OO, data.frame(un=khat),
                  labl,
                  "uncorrected estimate of %s",
                  "un")
  }

  if(any(correction == "border")) {
    ## Create list of r-close information
    id <- delta<=rmax
    close <- list(d = delta[id], i = row(delta)[id], j = col(delta)[id])
    ## border method
    ## Compute distances to boundary
    b <- bdist_spherepoints(X)
    I <- close$i
    bI <- b[I]
    ## apply reduced sample algorithm
    RS <- Kount(close$d, bI, b, breaks)
    numKb <- RS$numerator
    denKb <- lambda * RS$denom.count
    if(unnormalized){
      denKb / lambda
      numKb <- numKb * (n-1)
    }
    labl <- ifelse(!angular_arg, "hat(%s)[bord](r)", "hat(%s)[bord](theta)")
    OO <- bind.ratfv(OO,
                     data.frame(border=numKb),
                     data.frame(border=denKb),
                     labl,
                     "border-corrected estimate of %s",
                     "border",
                     ratio=ratio)
  }

  if(degrees) OO$r <- OO$r*360/(2*pi)

  return(OO)
}
