#' _bayesMig
#'
#' @title Bayesian Projection of Migration
#'
#' @description Collection of functions for making probabilistic projections of net migration rate for all countries of the world,
#' using a Bayesian hierarchical model (BHM) and the United Nations demographic time series.
#' Methodological details are provided in Azose & Raftery (2015) 
#'
#' @author Jon Azose, Hana Sevcikova and Adrian Raftery
#'
#' @maintainer: Hana Sevcikova <hanas@uw.edu>
#' 
#' @docType package
#' @import coda truncnorm wpp2019
#' @import grDevices graphics stats utils bayesTFR
#' @importFrom utils data
#' 
#' @examples
#' \dontrun{
#' #Run simulations and write to file
#' run.mig.mcmc(nr.chains=4, iter=10000, thin=10, replace.output = TRUE)
#' mig.diagnose(sim.dir='./bayesMig.output', burnin=100, thin=1)
#' 
#' #Predict and plot
#' mig.pred=mig.predict(burnin=100, replace.output = TRUE)
#' mig.trajectories.plot.all(mig.pred)
#' }
#' 
#' @references Azose, J. J., & Raftery, A. E. (2015). 
#' Bayesian probabilistic projection of international migration. Demography, 52(5), 1627-1650.
#' @name bayesMig-package
NULL
