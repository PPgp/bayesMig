#' @md
#' @title Bayesian Projection of Migration
#'
#' @description Collection of functions for making probabilistic projections of net migration rate for all countries of the world,
#' using a Bayesian hierarchical model (BHM) and the United Nations demographic time series. The model can be also applied
#' to user-defined data for other locations, such as subnational units.
#' Methodological details are provided in Azose & Raftery (2015).
#'
#' @author Jon Azose, Hana Sevcikova and Adrian Raftery
#'
#' Maintainer: Hana Sevcikova <hanas@uw.edu>
#' 
#' @docType package
#' @import coda truncnorm wpp2019
#' @import grDevices graphics stats utils bayesTFR
#' @importFrom utils data
#' 
#' @details The package is implemented in a similar way as the \pkg{bayesTFR} 
#'     package and thus, many functions have their equivalents in \pkg{bayesTFR}. 
#'     The main functions of the package are:
#' * \code{\link{run.mig.mcmc}}: Runs a Markov Chain Monte Carlo (MCMC) simulation.
#' It results in posterior samples of the model parameters.
#' * \code{\link{mig.predict}}: Using the posterior parameter samples, trajectories of future 
#' net migration rates are generated for all countries or given locations.
#' 
#' The following functions can be used to analyze results:
#' * \code{\link{mig.trajectories.plot}}: Shows the posterior trajectories for a given location, including the median and given probability intervals.
#' * \code{\link{mig.trajectories.table}}: Shows a tabular form of the posterior trajectories for a given location.
#' * \link{mig.map} and \link{mig.map.gvis}: Show a world map of migration rates 
#'     for a given projection or observed period, or for country-specific parameter estimates.
#' * \code{\link{mig.partraces.plot}} and \code{\link{mig.partraces.cs.plot}}: Plot the MCMC traces
#'  of country-independent parameters and country-specific parameters, respectively.
#' * \code{\link{mig.pardensity.plot}} and \code{\link{mig.pardensity.cs.plot}}: Plot the posterior density of the 
#'  country-independent parameters and country-specific parameters, respectively. 
#' * \code{\link{summary.bayesMig.mcmc.set}}: Summary method for the MCMC results.
#' * \code{\link{summary.bayesMig.prediction}}: Summary method for the prediction results.
#' 
#' For MCMC diagnostics, function \code{\link{mig.coda.list.mcmc}} creates an object of type 
#' \dQuote{mcmc.list} that can be used with the \pkg{coda} package. 
#' Furthermore, function \code{\link{mig.diagnose}} analyzes the MCMCs using the 
#' Raftery diagnostics implemented in the \pkg{coda} package and gives information 
#' about parameters that did not converge.
#' 
#' Existing results can be accessed using the \code{\link{get.mig.mcmc}} (estimation) and 
#' \code{\link{get.mig.prediction}} (prediction) functions. 
#' Existing convergence diagnostics can be accessed using the \code{\link{get.mig.convergence}} and
#' \code{\link{get.mig.convergence.all}} functions.
#' 
#' Historical data on migration rates are taken from the \pkg{wpp2019} (default), \pkg{wpp2022} or \pkg{wpp2017} package, 
#' depending on users settings. Alternatively, users can input their own data.
#' 
#' @examples
#' \dontrun{
#' # Run a real simulation (can take long time)
#' sim.dir <- tempfile()
#' m <- run.mig.mcmc(nr.chains = 4, iter = 10000, thin = 10, output.dir = sim.dir,
#'         verbose.iter = 1000)
#' 
#' # Prediction for all countries
#' pred <- mig.predict(sim.dir = sim.dir, nr.traj = 1000, burnin = 1000)
#' 
#' # Explore results
#' summary(pred, country = "Germany")
#' mig.trajectories.plot(pred, country = "Germany", nr.traj = 50)
#' 
#' # Check convergence diagnostics
#' mig.diagnose(sim.dir, burnin = 4000, thin = 1)
#' 
#' unlink(sim.dir, recursive = TRUE)
#' }
#' 
#' @references Azose, J. J., & Raftery, A. E. (2015). 
#' Bayesian probabilistic projection of international migration. Demography, 52(5), 1627-1650.
#' \doi{10.1007/s13524-015-0415-0}.
#' @name bayesMig-package
#' @aliases bayesMig
NULL
