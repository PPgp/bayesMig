#' @title MCMC Meta Object
#' @description Meta object \code{bayesMig.mcmc.meta} used by all chains of the same MCMC simulation. 
#' It contains information that is common to all chains. It is a part of a {\code{\link{bayesMig.mcmc.set}}} object.
#' @details The object is in standard cases not to be manipulated by itself.
#' @return Stores values of the various input arguments of the {\code{\link{run.mig.mcmc}}} function.
#' These are \code{nr.chains}, \code{start.year}, \code{present.year}, \code{wpp.year}, \code{my.e0.file}
#' @aliases bayesMig.mcmc.meta
#' @name mcmc-meta-object
NULL

#' @title MCMC Object
#' @description MCMC object \code{bayesMig.mcmc} containing information about one MCMC chain. 
#' A set of such objects belonging to the same simulation, together with a {\code{\link{bayesMig.mcmc.meta}}} object, 
#' constitute a {\code{\link{bayesMig.mcmc.set}}} object.
#' @aliases bayesMig.mcmc
#' @name mcmc-object
NULL
