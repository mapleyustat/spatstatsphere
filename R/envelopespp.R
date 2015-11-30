#  Envelopes for 'spp' objects (based on envelope.lpp in spatstat)
#
#' Envelope for Point Patterns on Sphere
#'
#' Enables envelopes to be computed for point patterns on a sphere.
#'
#' This is a method for the generic function \code{\link{envelope}} applicable
#' to point patterns on a sphere.
#'
#' The argument \code{Y} should be a point pattern on a sphere. The function
#' \code{fun} will be evaluated for the data and also for \code{nsim} simulated
#' point patterns on the same spherical domain. The upper and lower envelopes of
#' these evaluated functions will be computed as described in
#' \code{\link{envelope}}.
#'
#' The type of simulation is determined as follows.
#' \itemize{
#'   \item
#'     \code{simulate} is missing or \code{NULL}, then random point patterns
#'     will be generated according to a Poisson point process on the spherical
#'     domain on which \code{Y} is defined, with intensity estimated from
#'     \code{Y}.
#'   \item
#'     If \code{simulate} is present, it should be an expression that can be
#'     evaluated to yield random point patterns on the same spherical domain as
#'     \code{Y}.
#' }
#'
#' The function \code{fun} should accept as its first argument a point pattern
#' on a spherical domain (object of class \code{"spp"}) and should have another
#' argument called \code{r} or a \code{\dots{}} argument.
#'
#' @param Y A point pattern on a spherical domain (object of class \code{"spp"})
#' @param fun Function that is to be computed for each simulated pattern.
#' @param nsim Number of simulations to perform.
#' @param nrank Integer. Rank of the envelope value amongst the \code{nsim}
#'   simulated values. A rank of 1 means that the minimum and maximum simulated
#'   values will be used.
#' @param \dots Extra arguments passed to \code{fun}.
#' @param funargs A list, containing extra arguments to be passed to \code{fun}.
#' @param simulate Optional. Specifies how to generate the simulated point
#'   patterns.  If \code{simulate} is an expression in the R language, then this
#'   expression will be evaluated \code{nsim} times, to obtain \code{nsim} point
#'   patterns which are taken as the simulated patterns from which the envelopes
#'   are computed.  If \code{simulate} is a list of point patterns, then the
#'   entries in this list will be treated as the simulated patterns from which
#'   the envelopes are computed.  Alternatively \code{simulate} may be an object
#'   produced by the \code{envelope} command: see Details.
#' @param verbose Logical flag indicating whether to print progress reports
#'   during the simulations.
#' @param transform Optional. A transformation to be applied to the function
#'   values, before the envelopes are computed.  An expression object (see
#'   Details).
#' @param global Logical flag indicating whether envelopes should be pointwise
#'   (\code{global=FALSE}) or simultaneous (\code{global=TRUE}).
#' @param ginterval Optional.  A vector of length 2 specifying the interval of
#'   \eqn{r} values for the simultaneous critical envelopes. Only relevant if
#'   \code{global=TRUE}.
#' @param use.theory Logical value indicating whether to use the theoretical
#'   value, computed by \code{fun}, as the reference value for simultaneous
#'   envelopes. Applicable only when \code{global=TRUE}.
#' @param alternative Character string determining whether the envelope
#'   corresponds to a two-sided test (\code{side="two.sided"}, the default) or a
#'   one-sided test with a lower critical boundary (\code{side="less"}) or a
#'   one-sided test with an upper critical boundary (\code{side="greater"}).
#' @param scale Optional. Scaling function for global envelopes.  A function in
#'   the language which determines the relative scale of deviations, as a
#'   function of distance \eqn{r}, when computing the global envelopes.
#'   Applicable only when \code{global=TRUE}.  Summary function values for
#'   distance \code{r} will be \emph{divided} by \code{scale(r)} before the
#'   maximum deviation is computed. The resulting global envelopes will have
#'   width proportional to \code{scale(r)}.
#' @param clamp Logical value indicating how to compute envelopes when
#'   \code{alternative="less"} or \code{alternative="greater"}.  Deviations of
#'   the observed summary function from the theoretical summary function are
#'   initially evaluated as signed real numbers, with large positive values
#'   indicating consistency with the alternative hypothesis.  If
#'   \code{clamp=FALSE} (the default), these values are not changed.  If
#'   \code{clamp=TRUE}, any negative values are replaced by zero.
#' @param savefuns Logical flag indicating whether to save all the simulated
#'   function values.
#' @param savepatterns Logical flag indicating whether to save all the simulated
#'   point patterns.
#' @param nsim2 Number of extra simulated point patterns to be generated if it
#'   is necessary to use simulation to estimate the theoretical mean of the
#'   summary function. Only relevant when \code{global=TRUE} and the simulations
#'   are not based on CSR.
#' @param VARIANCE Logical. If \code{TRUE}, critical envelopes will be
#'   calculated as sample mean plus or minus \code{nSD} times sample standard
#'   deviation.
#' @param nSD Number of estimated standard deviations used to determine the
#'   critical envelopes, if \code{VARIANCE=TRUE}.
#' @param Yname Character string that should be used as the name of the data
#'   point pattern \code{Y} when printing or plotting the results.
#' @param do.pwrong Logical. If \code{TRUE}, the algorithm will also estimate
#'   the true significance level of the \dQuote{wrong} test (the test that
#'   declares the summary function for the data to be significant if it lies
#'   outside the \emph{pointwise} critical boundary at any point). This estimate
#'   is printed when the result is printed.
#' @param envir.simul Environment in which to evaluate the expression
#'   \code{simulate}, if not the current environment.
#' @return Function value table (object of class \code{"fv"}) with additional
#'   information, as described in \code{\link{envelope}}.
#' @seealso \code{\link{envelope}}, \code{\link{Kspp}}
#' @export
#' @keywords spatial
#' @examples
#'
#'    if(interactive()) {
#'      ns <- 39
#'      np <- 40
#'    } else { ns <- np <- 3 }
#'    X <- runifspp(np)
#'    # uniform Poisson
#'    envelope(X, nsim=ns)
#'
envelope.spp <-
  function(Y, fun=Kspp, nsim=99, nrank=1, ...,
           funargs=list(),
           simulate=NULL, verbose=TRUE,
           transform=NULL, global=FALSE, ginterval=NULL, use.theory=NULL,
           alternative=c("two.sided", "less", "greater"),
           scale=NULL, clamp=FALSE,
           savefuns=FALSE, savepatterns=FALSE, nsim2=nsim,
           VARIANCE=FALSE, nSD=2,
           Yname=NULL, do.pwrong=FALSE, envir.simul=NULL) {
  cl <- short.deparse(sys.call())
  if(is.null(Yname)) Yname <- short.deparse(substitute(Y))
  if(is.null(fun)) fun <- Kspp

  if("clipdata" %in% names(list(...)))
    stop(paste("The argument", sQuote("clipdata"),
               "is not available for envelope.spp"))

  envir.user <- if(!is.null(envir.simul)) envir.simul else parent.frame()
  envir.here <- sys.frame(sys.nframe())

  if(is.null(simulate)) {
    # ...................................................
    # Realisations of complete spatial randomness
    # will be generated by rpoisspp
    # Data pattern X is argument Y
    # Data pattern determines intensity of Poisson process
    X <- Y
    nY <- npoints(Y)
    dom <- domain(Y)
    Yintens <- nY/area(dom)
    # expression that will be evaluated
    simexpr <- expression(rpoisspp(Yintens, dom))
    # evaluate in THIS environment
    simrecipe <- simulrecipe(type = "csr",
                             expr = simexpr,
                             envir = envir.here,
                             csr   = TRUE)
  } else {
    # ...................................................
    # Simulations are determined by 'simulate' argument
    # Processing is deferred to envelopeEngine
    simrecipe <- simulate
    # Data pattern is argument Y
    X <- Y
  }
  envelopeEngine(X=X, fun=fun, simul=simrecipe,
                 nsim=nsim, nrank=nrank, ..., funargs=funargs,
                 verbose=verbose, clipdata=FALSE,
                 transform=transform,
                 global=global, ginterval=ginterval, use.theory=use.theory,
                 alternative=alternative, scale=scale, clamp=clamp,
                 savefuns=savefuns, savepatterns=savepatterns, nsim2=nsim2,
                 VARIANCE=VARIANCE, nSD=nSD,
                 Yname=Yname, cl=cl,
                 envir.user=envir.user, do.pwrong=do.pwrong,
                 foreignclass = "spp")
}
