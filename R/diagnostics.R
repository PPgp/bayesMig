###############
#MIGRATION
###############

#' @title MCMC convergence diagnostics
#'
#' @description Runs convergence diagnostics of existing migration Markov chains 
#'     using the \code{raftery.diag} function from the \code{coda} package.
#' 
#' @param sim.dir Directory with MCMC simulation results.
#' @param thin Thinning interval.
#' @param burnin Number of iterations to discard from the beginning of the parameter traces.
#' @param express Logical. If \code{TRUE}, the convergence diagnostic is run only on the country-independent
#' parameters. If \code{FALSE}, the country-specific parameters are included in the diagnostics. The number of
#' countries can be controlled by \code{country.sampling.prop}.
#' @param country.sampling.prop Proportion of countries to include in the diagnostics. If it is \code{NULL} and
#' \code{express=FALSE}, all countries are included. Setting a number between 0 and 1 will determine the proportion of countries
#' to be randomly sampled. For long Markov chains, this argument may significantly influence the runtime of this function.
#' @param keep.thin.mcmc Logical. If \code{TRUE}, the thinned traces used for computing the diagnostics are stored on disk.
#' @param verbose Logical value. Switches log messages on and off.
#' 
#' @details The \code{mig.diagnose} function invokes the \code{\link{mig.raftery.diag}} 
#'     function separately for country-independent parameters and for country-specific 
#'     parameters. It results in two possible states: red, i.e. it did not converge, and green, 
#'     i.e. it converged. The resulting object is stored in 
#'     \file{\{sim.dir\}/diagnostics/bayesMig.convergence_\{thin\}_\{burnin\}.rda} 
#'     and can be accessed using the function \code{\link{get.mig.convergence}}.
#'     
#'     Function \code{\link[bayesTFR]{has.mcmc.converged}} from the \pkg{bayesTFR} package 
#'     can be used to check if the existing diagnostics converged.
#' @return \code{mig.diagnose} returns an object of class \code{bayesMig.convergence} 
#'     containing summaries of the convergence check inputs and outputs. It has the 
#'     same structure as \code{\link[bayesTFR]{bayesTFR.convergence}}. 
#'     In addition it has an element \code{a.hw.est} which is the estimated value for 
#'     the \code{a.half.width} argument in \code{\link{run.mig.mcmc}}.
#' @seealso \code{\link[bayesTFR]{tfr.raftery.diag}}, \code{\link[coda]{raftery.diag}}
#' @examples
#' \dontrun{
#' mig.diagnose(sim.dir='./bayesMig.output', burnin=100, thin=1)
#' }
#' @aliases bayesMig.convergence
#' @rdname diagnose
#' @export

mig.diagnose <- function(sim.dir, thin=80, burnin=2000, express=FALSE, 
                         country.sampling.prop=NULL, keep.thin.mcmc=FALSE, verbose=TRUE) {
  diag <- bayesTFR:::.do.diagnose(type='mig', class.name='bayesMig.convergence', 
                         sim.dir=sim.dir, thin=thin, burnin=burnin, express=express,
                         country.sampling.prop=country.sampling.prop, keep.thin.mcmc=keep.thin.mcmc,	verbose=verbose)
  m <- get.mig.mcmc(sim.dir)
  diag$a.hw.est <- estimate.a.hw(m, burnin = burnin, thin = thin)
  invisible(diag)
}

#' @param mcmc A \code{\link{bayesMig.mcmc}} or \code{\link{bayesMig.mcmc.set}} object. 
#'     If not given, the object is loaded from the simulation directory given by 
#'     \code{sim.dir}.
#' @param country Name or code of a country. If it is given, only country-specific 
#'     parameters parameters of that country are considered.
#' @param par.names Names of country-independent parameters for which the Raftery 
#'     diagnostics should be computed. By default all parameters are used.
#' @param par.names.cs Names of country-specific parameters for which the Raftery 
#'     diagnostics should be computed. By default all parameters are used.
#' @param \dots Additional arguments passed to the \code{\link{mig.coda.list.mcmc}} function.
#'     
#' @details For details on the \code{mig.raftery.diag} function, see \code{\link[bayesTFR]{tfr.raftery.diag}}. 
#' @rdname diagnose
#' @export
mig.raftery.diag <- function(mcmc=NULL, 
                             sim.dir=file.path(getwd(), 'bayesMig.output'),
                             burnin=0, country=NULL,
                             par.names = NULL,
                             par.names.cs = NULL,
                             country.sampling.prop=1,
                             verbose=TRUE, ...
                            ) {
  mcmc.set <- if (is.null(mcmc)) get.mig.mcmc(sim.dir = sim.dir) else mcmc
  if(bayesTFR:::is.missing(par.names)) 
    par.names <- mig.parameter.names()
  if(bayesTFR:::is.missing(par.names.cs)) 
    par.names.cs <- mig.parameter.names.cs()
  return(bayesTFR::tfr.raftery.diag(mcmc = mcmc.set, sim.dir = sim.dir, burnin = burnin,
                                    country = country, par.names = par.names, par.names.cs = par.names.cs,
                                    country.sampling.prop = country.sampling.prop, verbose = verbose, ...))
}

#' @details The \code{estimate.a.hw} function estimates an optimal value for the \code{a.half.width}
#'     argument in \code{\link{run.mig.mcmc}}. 
#' @rdname diagnose
#' @export
estimate.a.hw <- function(mcmc, burnin = 0, thin = NULL) {
  tr <- get.mig.parameter.traces(mcmc, burnin = burnin, par.names = c("a", "b"))
  return(sd(tr[, "a"])*sqrt(1-abs(cor(tr[, "a"], tr[, "b"]))^2)*2.3 * 2)
}
