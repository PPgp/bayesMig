###############
#MIGRATION
###############

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
  

#' @title MCMC convergence diagnostics
#'
#' @description Runs convergence diagnostics of existing migration Markov chains using the \code{raftery.diag} function from the \code{coda} package.
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
#' @return An object of class \code{bayesMig.convergence} containing summaries of the convergence check inputs and outputs
#' @examples
#' \dontrun{
#' mig.diagnose(sim.dir='./bayesMig.output', burnin=100, thin=1)
#' }
#' @export

mig.diagnose <- function(sim.dir, thin=80, burnin=2000, express=FALSE, 
                         country.sampling.prop=NULL, keep.thin.mcmc=FALSE, verbose=TRUE) {
  invisible(bayesTFR:::.do.diagnose(type='mig', class.name='bayesMig.convergence', 
                         sim.dir=sim.dir, thin=thin, burnin=burnin, express=express,
                         country.sampling.prop=country.sampling.prop, keep.thin.mcmc=keep.thin.mcmc,	verbose=verbose))
}

