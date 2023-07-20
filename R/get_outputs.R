
##################
#MIGRATION
##################


#' @title Access MCMC results
#'
#' @description Function \code{get.mig.mcmc} retrieves results of an MCMC simulation and creates an object of class
#' \code{\link{bayesMig.mcmc.set}}. Function \code{has.mig.mcmc} checks the existence of such results.
#' 
#' @param sim.dir Directory where simulation results are stored.
#' @param chain.ids Chain identifiers in case only specific chains should be included
#'     in the resulting object. By default, all available chains are included.
#' @param low.memory Logical. If \code{FALSE} full MCMC traces are loaded into memory.
#' @param burnin Burn-in used for loading traces. Only relevant, if \code{low.memory=FALSE}.
#' @param verbose Logical value. Switches log messages on and off.
#' 
#' @seealso \code{\link{run.mig.mcmc}}
#' 
#' @examples
#' \dontrun{
#' # Toy simulation
#' sim.dir <- tempfile()
#' m <- run.mig.mcmc(nr.chains = 1, iter = 10, output.dir = sim.dir)
#' 
#' # can be later accessed via
#' m <- get.mig.mcmc(sim.dir)
#' summary(m)
#' 
#' has.mig.mcmc(sim.dir) # should be TRUE
#' 
#' unlink(sim.dir, recursive = TRUE)
#' }
#' @export
#' @rdname get-mcmc

get.mig.mcmc <- function(sim.dir=file.path(getwd(), 'bayesMig.output'), chain.ids=NULL,
                         low.memory = TRUE, burnin = 0, verbose = FALSE) {
  ############
  # Returns an object of class bayesMig.mcmc.set
  ############
  mcmc.file.path <- file.path(sim.dir, 'bayesMig.mcmc.meta.rda')
  if(!file.exists(mcmc.file.path)) {
    warning('File ', mcmc.file.path, ' does not exist.')
    return(NULL)
  }
  load(file=mcmc.file.path)
  bayesMig.mcmc.meta$output.dir <- normalizePath(sim.dir)
  if(is.null(bayesMig.mcmc.meta$user.data)) bayesMig.mcmc.meta$user.data <- FALSE
  
  if (is.null(chain.ids)) {
    mc.dirs.short <- list.files(sim.dir, pattern='^mc[0-9]+', full.names=FALSE)
    chain.ids <- as.integer(substring(mc.dirs.short, 3))
  } else {
    mc.dirs.short <- paste('mc', chain.ids, sep='')
  }
  ord.idx <- order(chain.ids)
  mc.dirs.short <- mc.dirs.short[ord.idx]
  chain.ids <- chain.ids[ord.idx]
  mcmc.chains <- list()
  counter<-1
  for (imc.d in chain.ids) {
    if (verbose)
      cat('Loading chain', imc.d, 'from disk. ')
    bayesMig.mcmc <- local({
      load(file=file.path(sim.dir, mc.dirs.short[counter], 'bayesMig.mcmc.rda'))
      bayesMig.mcmc})
    mc <- c(bayesMig.mcmc, list(meta=bayesMig.mcmc.meta))
    class(mc) <- class(bayesMig.mcmc)
    if (!low.memory) { # load full mcmc traces
      #load full mcmc traces
      th.burnin <- bayesTFR:::get.thinned.burnin(mc, burnin)
      mc$traces <- load.mig.parameter.traces.all(mc, burnin=th.burnin)
    } else {
      th.burnin <- 0 # low memory (traces will be loaded as they are needed)
      mc$traces <- 0
    }
    mc$traces.burnin <- th.burnin
    
    mc$output.dir <- mc.dirs.short[counter]
    if (verbose)
      cat('(mcmc.list[[', counter, ']]).\n')
    mcmc.chains[[counter]] <- mc
    counter <- counter+1
  }
  names(mcmc.chains) <- chain.ids
  return(structure(list(meta=bayesMig.mcmc.meta, 
                        mcmc.list=mcmc.chains), class='bayesMig.mcmc.set'))
}

#' @title Access Prediction Object
#'
#' @description Function \code{get.mig.prediction} retrieves results of a prediction and creates an object of class
#'     \code{\link{bayesMig.prediction}}. Function \code{has.mig.prediction} checks an existence of such results.
#' 
#' @param mcmc Object of class \code{\link{bayesMig.mcmc.set}} used to make the prediction. 
#'     If it is \code{NULL}, the prediction is loaded from directory given by \code{sim.dir}.
#' @param sim.dir Directory where the prediction is stored.
#' @param mcmc.dir Optional argument to be used only in a special case when the mcmc object 
#'     contained in the prediction object was estimated in different directory than in the one 
#'     to which it points to (for example due to moving or renaming the original directory). 
#'     The argument causes that the mcmc is redirected to the given directory. 
#'     It can be set to \code{NA} if no loading of the mcmc object is desired.
#'     
#' @details If \code{mcmc} is not \code{NULL}, the search directory is set to 
#'     \code{mcmc$meta$output.dir}. This approach assumes that the prediction was stored in the 
#'     same directory as the MCMC simulation, i.e. the \code{output.dir} argument of the 
#'     \code{\link{mig.predict}} function was set to \code{NULL}. If it is not the case, 
#'     the argument \code{mcmc.dir} should be used.
#'     
#' @return Function \code{get.mig.prediction} returns an object of class \code{\link{bayesMig.prediction}}.
#' @export
#' @rdname get-prediction
#' 
get.mig.prediction <- function(mcmc=NULL, sim.dir=NULL, mcmc.dir=NULL) {
  ############
  # Returns an object of class bayesMig.prediction
  # Set mcmc.dir to NA, if the prediction object should not have a pointer 
  # to the corresponding mcmc traces
  ############
  if (!is.null(mcmc)) 
    sim.dir <- if(is.character(mcmc)) mcmc else mcmc$meta$output.dir
  if (is.null(sim.dir)) stop('Either mcmc or directory must be given.')
  output.dir <- file.path(sim.dir, 'predictions')
  pred.file <- file.path(output.dir, 'prediction.rda')
  if(!file.exists(pred.file)) {
    warning('File ', pred.file, ' does not exist.')
    return(NULL)
  }
  load(file=pred.file)
  bayesMig.prediction$output.directory <- output.dir

  pred <- bayesMig.prediction
  # re-route mcmcs if necessary
  if(!is.null(mcmc.dir) || !has.mig.mcmc(pred$mcmc.set$meta$output.dir)) {
    if((!is.null(mcmc.dir) && !is.na(mcmc.dir)) || is.null(mcmc.dir)) {
      new.path <- file.path(sim.dir, basename(pred$mcmc.set$meta$output.dir))
      if (has.mig.mcmc(new.path)) pred$mcmc.set <- get.mig.mcmc(new.path)
      else {
        est.dir <- if(is.null(mcmc.dir)) sim.dir else mcmc.dir
        pred$mcmc.set <- get.mig.mcmc(est.dir)
      }
    }
  }
  if(is.null(pred$nr.imputed)) pred$nr.imputed <- rep(0, pred$mcmc.set$meta$nr.countries)
  return(pred)
}

#' @rdname get-mcmc
#' @export
#' 
has.mig.mcmc <- function(sim.dir) {
  return(file.exists(file.path(sim.dir, 'bayesMig.mcmc.meta.rda')))
}


#' @title Accessing Convergence Diagnostics Object
#'
#' @description The function retrieves results of convergence diagnostics 
#'     (created by \code{\link{mig.diagnose}}) from disk.
#' 
#' @param sim.dir Simulation directory used for computing the diagnostics.
#' @param thin Thinning interval used with this diagnostics.
#' @param burnin Burn-in used for computing the diagnostics.
#' 
#' @details Function \code{get.mig.convergence} loads an object of class 
#'     \code{\link{bayesMig.convergence}} for the specific \code{thin} and \code{burnin} 
#'     used in \code{\link{mig.diagnose}} to generate this object. 
#'     Function \code{get.mig.convergence.all} loads all \code{\link{bayesMig.convergence}} objects 
#'     available in \code{sim.dir}.
#'     
#' @return \code{get.mig.convergence} returns an object of class \code{\link{bayesMig.convergence}}; \cr
#'     \code{get.mig.convergence.all} returns a list of objects of class \code{\link{bayesMig.convergence}}. 
#'     
#' @seealso \code{\link{mig.diagnose}}, \code{\link{summary.bayesMig.convergence}}
#' 
#' @examples 
#' \dontrun{
#' # Run a real simulation (can take long time)
#' sim.dir <- tempfile()
#' m <- run.mig.mcmc(nr.chains = 2, iter = 10000, thin = 10, output.dir = sim.dir)
#' 
#' # Run convergence diagnostics with different burning and thin
#' mig.diagnose(sim.dir, burnin = 1000, thin = 2)
#' mig.diagnose(sim.dir, burnin = 500, thin = 1)
#' 
#' diags <- get.mig.convergence.all(sim.dir)
#' for(i in 1:length(diags))
#'     print(summary(diags[[i]]))
#' 
#' unlink(sim.dir, recursive = TRUE)
#' }
#' @rdname get-convergence
#' @export
#' 
get.mig.convergence <- function(sim.dir=file.path(getwd(), 'bayesMig.output'), 
                               thin=225, burnin=10000) {
  file.name <- file.path(sim.dir, 'diagnostics', paste('bayesMig.convergence_', 
                                                       thin, '_', burnin, '.rda', sep=''))
  if(!file.exists(file.name)){
    warning('Convergence diagnostics in ', sim.dir, ' for burnin=', burnin, 
            ' and thin=', thin, ' does not exist.')
    return(NULL)
  }
  bayesMig.convergence <- local({load(file.name)
                                  bayesMig.convergence})
  return(bayesMig.convergence)
}

#' @rdname get-convergence
#' @export
#' 
get.mig.convergence.all <- function(sim.dir=file.path(getwd(), 'bayesMig.output')) {
    return(bayesTFR:::.do.get.convergence.all('mig', 'bayesMig', sim.dir=sim.dir))
}

get.thinned.burnin <- function(mcmc, burnin) {
  if (burnin==0) return(0)
  if (mcmc$thin == 1) return(burnin)
  return(1 + if(burnin >= mcmc$thin) length(seq(mcmc$thin, burnin, by=mcmc$thin)) else 0)
}


load.mig.parameter.traces.all <- function(mcmc, par.names=mig.parameter.names(), 
                                          par.names.cs=mig.parameter.names.cs(),
                                          burnin=0, thinning.index=NULL) {
  result <- load.mig.parameter.traces(mcmc, par.names, burnin=burnin, thinning.index=thinning.index)
  if(!is.null(par.names.cs))
    for (country in get.countries.index(mcmc$meta)) {
      result <- cbind(result, 
                      load.mig.parameter.traces.cs(mcmc, 
                                                   get.country.object(country, 
                                                                      mcmc$meta, index=TRUE)$code, 
                                                   par.names.cs, burnin=burnin,
                                                   thinning.index=thinning.index))
    }
  return (result)
}

#' @title Accessing Parameter Names
#'
#' @description Functions for accessing names of the MCMC parameters, 
#'     either country-independent or country-specific.
#' 
#' @return \code{mig.parameter.names} returns names of the world parameters.
#' 
#' @examples 
#' mig.parameter.names()
#' 
#' @export
#' @rdname parnames
#' 
mig.parameter.names <- function() {
  # Return all country-independent parameter names. 
  return(c("a","b","mu_global","sigma2_mu"))
}

#' @param country.code Country code. If it is given, the country-specific parameter names contain 
#'     the suffix \sQuote{_country\eqn{X}} where \eqn{X} is the \code{country.code}.
#'     
#' @return \code{mig.parameter.names.cs} returns names of the country-specific parameters.
#' 
#' @examples 
#' mig.parameter.names.cs()
#' 
#' @export
#' @rdname parnames
#' 
mig.parameter.names.cs <- function(country.code=NULL) {
  #Return all country-specific parameter names. 
  #If country is not NULL, it must be a country code.
  #It is attached to the parameter name.
  par.names <- c("mu_c","phi_c","sigma2_c")
  if (is.null(country.code))
    return(par.names)
  return(paste0(par.names, '_country', country.code))
}

load.mig.parameter.traces <- function(mcmc, par.names=NULL, burnin=0, thinning.index=NULL) {
  if(missing(par.names)) par.names <- mig.parameter.names()
  return(bdem.parameter.traces(mcmc, par.names, burnin=burnin, thinning.index=thinning.index))
}

load.mig.parameter.traces.cs <- function(mcmc, country, par.names=NULL, burnin=0, 
                                         thinning.index=NULL) {
  #The par.names input should look something like c("mu_c","phi_c","sigma2_c")
  if(missing(par.names)) par.names <- mig.parameter.names.cs()
  return(bdem.parameter.traces(mcmc, par.names, paste0("_country",country),
                               burnin=burnin, thinning.index=thinning.index))
}

#' @export
bdem.parameter.traces.bayesMig.mcmc <- function(mcmc, par.names, ...) {
  # Load traces from the disk
  all.standard.names <- c(mig.parameter.names(), mig.parameter.names.cs())
  return(bayesTFR:::.do.get.traces(mcmc, par.names=par.names, ..., all.standard.names = all.standard.names))
}

#' @export
coda.mcmc.bayesMig.mcmc <- function(mcmc, country=NULL, par.names=NULL, 
                                    par.names.cs=NULL, burnin=0, thin=1, ...
                                    ) {
  if(missing(par.names)) par.names <- mig.parameter.names()
  if(missing(par.names.cs)) par.names.cs <- mig.parameter.names.cs()
  return(bayesTFR:::coda.mcmc.bayesTFR.mcmc(mcmc, country = country, par.names = par.names, 
                                            par.names.cs = par.names.cs, ...))
}

#' @title Conversion to coda-formatted objects
#' 
#' @description The functions convert MCMC traces (simulated using \code{\link{run.mig.mcmc}}) into 
#' objects that can be used with the \pkg{coda} package.
#' @param mcmc.list A list of objects of class \code{bayesMig.mcmc}, or an object of class \code{\link{bayesMig.mcmc.set}} or \code{\link{bayesMig.prediction}}.
#' If \code{NULL}, the MCMCs are
#' loaded from \code{sim.dir}. Either \code{mcmc} or \code{sim.dir} must be given.
#' @param country Country name or code. Used in connection with the \code{par.names.cs} argument
#' (see below).
#' @param chain.ids Vector of chain identifiers. By default, all chains available in the \code{mcmc.list}
#' object are included.
#' @param sim.dir Directory with the MCMC simulation results. Only used if \code{mcmc.list} is \code{NULL}.
#' @param par.names Names of country-independent parameters to be included. Default names are
#' those returned by the \code{mig.parameter.names} function, which includes all country-independent
#' parameters in the BHM.
#' @param par.names.cs Names of country-specific parameters to be included. The argument \code{country}
#' is used to filter out traces that correspond to a specific country. If \code{country} is not given, 
#' traces of each parameter are given for all countries. Default names are those returned by 
#' \code{mig.parameter.names.cs()}, which includes all country-specific parameters in the BHM.
#' @param low.memory Logical indicating if the function should run in a memory-efficient mode.
#' @param \dots Additional arguments passed to the \pkg{coda}'s \code{\link[coda]{mcmc}} function, such as \code{burnin} and \code{thin}.
#' 
#' @return Returns an object of class \dQuote{mcmc.list} defined in the \pkg{coda} package.
#' @export

mig.coda.list.mcmc <- function(mcmc.list = NULL, country = NULL, chain.ids = NULL,
                              sim.dir = file.path(getwd(), 'bayesMig.output'), 
                              par.names = NULL, par.names.cs = NULL, 
                              low.memory = FALSE, ...) {
  # return a list of mcmc objects that can be analyzed using the coda package
  if (is.null(mcmc.list)) {
    mcmc.list <- get.mig.mcmc(sim.dir, chain.ids=chain.ids, low.memory = low.memory)$mcmc.list
  } else {
    mcmc.list <- get.mcmc.list(mcmc.list)
    if (!is.null(chain.ids)) {
      mcmc.list <- mcmc.list[chain.ids]
    }
  }
  if(missing(par.names)) par.names <- mig.parameter.names()
  if(missing(par.names.cs)) par.names.cs <- mig.parameter.names.cs()
  result <- list()
  i <- 1
  for(mcmc in mcmc.list) {
    result[[i]] <- coda.mcmc(mcmc, country = country, par.names = par.names, 
                             par.names.cs = par.names.cs, ...)
    i <- i+1
  }
  return(mcmc.list(result))
}

get.full.par.names.cs <- function(par.names, full.par.names, country=NULL, index=FALSE) {
  # Return full name of par.names that are included in full.par.names
  # which are suppose to be all country-specific parameters.
  # E.g. for 'phi_c', it would return 'phi_c_country1', ..., 'phi_c_country233'
  # If index is TRUE, return the index of the matches.
  result <- c()	
  for (name in par.names) {
    pattern <- paste('^',name,'_.*',sep='')
    if (!is.null(country)) {
      pattern <- paste(pattern, 'country', country, '$', sep='')
    } else {
      pattern <- paste(pattern, '[[:digit:]]$', sep='')
    }
    result <- c(result, grep(pattern, full.par.names, value=!index))
  }
  return(result)
}

get.full.par.names <- function(par.names, full.par.names, index=FALSE) {
  # Return full name of par.names that are included in full.par.names
  # which are suppose to be all country-independent parameters.
  # E.g. for 'a', it would return just 'a' (and not 'alpha')
  # If index is TRUE, return the index of the matches.
  result <- c()	
  for (name in par.names) {
    pattern <- paste('^', name, '$', sep='') # matches exactly
    result <- c(result, grep(pattern, full.par.names, value=!index))
  }
  return(result)
}

get.burned.mig.traces <- function(mcmc, par.names, burnin=0, thinning.index=NULL) {
  # get traces that are already loaded in the mcmc object
  traces <- matrix(mcmc$traces[, par.names],ncol=length(par.names), dimnames=list(NULL, par.names))
  discard <- burnin - mcmc$traces.burnin
  if (discard > 0)
    traces <- traces[-seq(1, discard),, drop = FALSE]
  if(!is.null(thinning.index))
    traces <- traces[thinning.index, , drop = FALSE]
  return(traces)
}

#' @export
get.mcmc.list.bayesMig.mcmc.set <- function(mcmc.list, ...) return(mcmc.list$mcmc.list)
#' @export
get.mcmc.list.bayesMig.mcmc <- function(mcmc.list, ...) return(list(mcmc.list))
#' @export
get.mcmc.list.bayesMig.prediction <- function(mcmc.list, ...) return(mcmc.list$mcmc.set$mcmc.list)


get.thinned.mig.mcmc <- function(mcmc.set, thin=1, burnin=0) {
  dir.name <- file.path(mcmc.set$meta$output.dir, paste('thinned_mcmc', thin, burnin, sep='_'))
  if(file.exists(dir.name)) return(get.mig.mcmc(dir.name))
  return(NULL)
}

#' @return Function \code{has.mig.prediction} returns  a logical indicating if a prediction exists.
#' @rdname get-prediction
#' @export
#' 
has.mig.prediction <- function(mcmc=NULL, sim.dir=NULL) {
  if (!is.null(mcmc)) sim.dir <- if(is.character(mcmc)) mcmc else mcmc$meta$output.dir
  if (is.null(sim.dir)) stop('Either mcmc or directory must be given.')
  if(file.exists(file.path(sim.dir, 'predictions', 'prediction.rda'))) return(TRUE)
  return(FALSE)
}

#' @title Accessing MCMC Parameter Traces
#' @description Functions for accessing traces of the MCMC parameters, either country-independent 
#'     or country-specific.
#' @param mcmc.list List of \code{\link{bayesMig.mcmc}} objects.
#' @param par.names Names of country-independent parameters (in case of 
#'     \code{get.mig.parameter.traces}) or country-specific parameters 
#'     (in case of \code{get.mig.parameter.traces.cs}) to be included. 
#'     By default all parameters are included, given either by \code{\link{mig.parameter.names}()} 
#'     (for global parameters) or \code{\link{mig.parameter.names.cs}()} (for location-specific parameters).
#' @param burnin Burn-in indicating how many iterations should be removed 
#'     from the beginning of each chain.
#' @param thinning.index Index of the traces for thinning. If it is \code{NULL}, 
#'     \code{thin} is used. \code{thinning.index} does not include \code{burnin} 
#'     and should be flattened over all chains. For example, if there are two MCMC 
#'     chains of length 1000, \code{burnin=200} and we want an equidistantly spaced 
#'     sample of length 400, then the value should be \cr
#'     \code{thinning.index = seq(1, 1600, length = 400)}.
#' @param thin An integer value for thinning. It is an alternative to 
#'     \code{thinning.index}. The above example is equivalent to \code{thin=4}.
#'     
#' @return Both functions return a matrix with columns being the parameters and 
#'     rows being the MCMC values, attached to one another in case of multiple chains. 
#'     \code{get.mig.parameter.traces} returns country-independent parameters, 
#'     \code{get.mig.parameter.traces.cs} returns country-specific parameters.
#'     
#' @seealso \code{\link{mig.coda.list.mcmc}} for another way of retrieving parameter traces;
#'     \code{\link{mig.parameter.names}} and \code{\link{mig.parameter.names.cs}} for parameter names.
#' 
#' @examples 
#' # Toy simulation for US states
#' us.mig.file <- file.path(find.package("bayesMig"), "extdata", "USmigrates.txt")
#' sim.dir <- tempfile()
#' m <- run.mig.mcmc(nr.chains = 2, iter = 30, thin = 1, my.mig.file = us.mig.file, 
#'         output.dir = sim.dir, present.year = 2017, annual = TRUE)
#' # obtain traces of hierarchical parameters     
#' par.values <- get.mig.parameter.traces(m$mcmc.list, burnin = 5)
#' dim(par.values) # matrix 50 x 4
#' hist(par.values[, "mu_global"], main = "mu")
#' 
#' # obtain traces of location-specific traces for California
#' mig.parameter.names.cs() # allowed parameter names 
#' par.values.cs <- get.mig.parameter.traces.cs(m$mcmc.list, 
#'         country.obj = get.country.object("California", meta = m$meta),
#'         burnin = 5, par.names = "phi_c")
#' dim(par.values.cs) # matrix 50 x 1
#' hist(par.values.cs, main = colnames(par.values.cs))
#' unlink(sim.dir, recursive = TRUE)
#'     
#' @export
#' @rdname get-traces
#' 
get.mig.parameter.traces <- function(mcmc.list, par.names=NULL, 
                                     burnin=0, thinning.index=NULL, thin=NULL) {
  # get parameter traces either from disk or from memory, if they were already loaded
  mcmc.list <- get.mcmc.list(mcmc.list)
  if(missing(par.names)) par.names <- mig.parameter.names()
  return(bayesTFR:::do.get.tfr.parameter.traces(is.cs=FALSE, mcmc.list=mcmc.list, par.names=par.names, 
                                     burnin=burnin, thinning.index=thinning.index, thin=thin))
}

#' @param country.obj Country object (see \code{\link[bayesTFR]{get.country.object}}).
#' @export
#' @rdname get-traces
get.mig.parameter.traces.cs <- function(mcmc.list, country.obj, par.names=NULL, 
                                        burnin=0, thinning.index=NULL, thin=NULL) {
  # country.obj is result of get.country.object()
  # get traces for country-specific parameters either from disk or from memory, if they were already loaded
  mcmc.list <- get.mcmc.list(mcmc.list)
  if(missing(par.names)) par.names <- mig.parameter.names.cs()
  return(bayesTFR:::do.get.tfr.parameter.traces(is.cs=TRUE, mcmc.list=mcmc.list, par.names=par.names, country.obj=country.obj,
                                     burnin=burnin, thinning.index=thinning.index, thin=thin))
}


create.thinned.mig.mcmc <- function(mcmc.set, thin=1, burnin=0, output.dir=NULL, verbose=TRUE) {
  #Return a thinned mcmc.set object with burnin removed and all chanins collapsed into one
  mcthin <- max(sapply(mcmc.set$mcmc.list, function(x) x$thin))
  thin <- max(c(thin, mcthin))
  meta <- mcmc.set$meta
  total.iter <- get.stored.mcmc.length(mcmc.set$mcmc.list, burnin=burnin)
  meta$is.thinned <- TRUE
  meta$parent.iter <- get.total.iterations(mcmc.set$mcmc.list, burnin)
  meta$parent.meta <- mcmc.set$meta
  meta$nr.chains <- 1
  
  if(verbose) cat('\nStoring thinned mcmc:')
  # store the meta object
  meta$output.dir <- file.path(
    if(is.null(output.dir)) meta$output.dir else output.dir, 
    paste('thinned_mcmc', thin, burnin, sep='_'))
  if(!file.exists(meta$output.dir)) 
    dir.create(meta$output.dir, recursive=TRUE)
  store.bayesMig.meta.object(meta, meta$output.dir)
  
  thin.index <- if(thin > mcthin) unique(round(seq(1, total.iter, by=thin/mcthin))) else 1:total.iter
  nr.points <- length(thin.index)
  
  #create one collapsed mcmc
  thinned.mcmc <- mcmc.set$mcmc.list[[1]]
  thinned.mcmc$meta <- meta
  thinned.mcmc$thin <- 1
  thinned.mcmc$id <- 1
  thinned.mcmc$traces <- 0
  thinned.mcmc$length <- nr.points
  thinned.mcmc$finished.iter <- nr.points
  thinned.mcmc$output.dir <- 'mc1'	
  outdir.thin.mcmc <- file.path(meta$output.dir, 'mc1')
  if(!file.exists(outdir.thin.mcmc)) dir.create(outdir.thin.mcmc)
  
  store.bayesMig.object(thinned.mcmc, outdir.thin.mcmc)
  
  if(verbose) cat('\nStoring country-independent parameters ...')
  for (par in mig.parameter.names()) {
    values <- get.mig.parameter.traces(mcmc.set$mcmc.list, par, burnin,
                                       thinning.index=thin.index)
    #values = as.vector(t(values))#JA: This seems to be necessary to switch from multiple rows to a single vector.
    bayesTFR:::write.values.into.file.cindep(par, values, outdir.thin.mcmc, compression.type='None')
  }
  if(verbose) cat('done.\nStoring country-specific parameters ...')
  par.names.cs <- mig.parameter.names.cs()
  for (country in 1:mcmc.set$meta$nr.countries){
    country.obj <- get.country.object(country, mcmc.set$meta, index=TRUE)
    for (par in par.names.cs) {
      values <- get.mig.parameter.traces.cs(mcmc.set$mcmc.list, country.obj, par, 
                                            burnin=burnin, thinning.index=thin.index)
      #values = as.vector(t(values))#JA: This seems to be necessary to switch from multiple rows to a single vector.
      bayesTFR:::write.values.into.file.cdep(par, values, outdir.thin.mcmc, country.code=country.obj$code,
                                  compression.type='None')
    }
  }
  if(verbose) cat('done.\n')
  #JA: We'll assume all countries are used for estimation.
  # if (mcmc.set$meta$nr_countries > mcmc.set$meta$nr_countries_estimation) {
  #   .update.thinned.extras(mcmc.set, (mcmc.set$meta$nr_countries_estimation+1):mcmc.set$meta$nr_countries,
  #                          burnin=burnin, nr.points=nr.points, dir=outdir.thin.mcmc, verbose=verbose)
  # }
  invisible(structure(list(meta=meta, mcmc.list=list(thinned.mcmc)), class='bayesMig.mcmc.set'))
}


#' @export
#' @rdname summary-mcmc
#' 
summary.bayesMig.mcmc <- function(object, country = NULL, 
                                   par.names = NULL, par.names.cs = NULL, thin = 1, burnin = 0, ...) {
  res <- list()
  class(res) <- "summary.bayesTFR.mcmc"
  if(missing(par.names)) par.names <- mig.parameter.names()
  if(missing(par.names.cs)) par.names.cs <- mig.parameter.names.cs()
  if (!is.null(country)) {
    country.obj <- get.country.object(country, object$meta)
    if(is.null(country.obj$name)) stop("Country ", country, " not found.")
    res$country.name <- country.obj$name
    country <- country.obj$code
  } 
  res$results <- summary(coda.mcmc(object, country=country, par.names=par.names,
                                   par.names.cs=par.names.cs, thin=thin, burnin=burnin), ...)
  return(res)
}

#' @title Summary Statistics for Migration Markov Chain Monte Carlo 
#'
#' @description Summary of an object \code{\link{bayesMig.mcmc.set}} or \code{\link{bayesMig.mcmc}},
#'     computed via \code{\link{run.mig.mcmc}}.  It can be obtained either for all countries or 
#'     for a specific country, and either for all parameters or for specific parameters.
#'     The function uses the \code{\link[coda]{summary.mcmc}} function of the \pkg{coda} package.
#' 
#' @param object Object of class \code{\link{bayesMig.mcmc.set}} or \code{\link{bayesMig.mcmc}}.
#' @param country Country name or code if a country-specific summary is desired. 
#'     The code can be either numeric or ISO-2 or ISO-3 characters. By default, summary 
#'     for all countries is generated.
#' @param chain.id Identifiers of MCMC chains. By default, all chains are considered.
#' @param par.names Country independent parameters to be included in the summary. 
#'     The default names are given by \code{\link{mig.parameter.names}()}.
#' @param par.names.cs Country-specific parameters to be included in the summary.
#'     The default names are given by \code{\link{mig.parameter.names.cs}()}.
#' @param meta.only Logical. If it is \code{TRUE}, only meta information of the simulation is included.
#' @param thin Thinning interval. Only used if larger than the \code{thin} argument used in \code{\link{run.mig.mcmc}}.
#' @param burnin Number of iterations to be discarded from the beginning of each chain before computing the summary.
#' @param \dots Additional arguments passed to the \code{\link[coda]{summary.mcmc}} function of the \pkg{coda} package.
#' @export
#' @rdname summary-mcmc
#' 
summary.bayesMig.mcmc.set <- function(object, country=NULL, chain.id=NULL, 
                                       par.names = NULL, 
                                       par.names.cs = NULL, 
                                       meta.only=FALSE, thin=1, burnin=0, ...) {
  if(missing(par.names)) par.names <- mig.parameter.names()
  if(missing(par.names.cs) && !is.null(country)) par.names.cs <- mig.parameter.names.cs()
  
  res <- list(meta = summary(object$meta))
  class(res) <- "summary.bayesMig.mcmc.set"
  if(meta.only) {
    res$chain.info <- bayesTFR:::chain.info(object)
    return(res)
  }
  if (!is.null(chain.id)) {
    res$mcmc <- summary(object$mcmc.list[[chain.id]], country = country, 
                        par.names = par.names, par.names.cs = par.names.cs, 
                        thin = thin, burnin = burnin, ...)
    return(res)
  }
  if (!is.null(country)) {
    country.obj <- get.country.object(country, object$meta)
    if(is.null(country.obj$name)) stop("Country ", country, " not found.")
    res$country.name <- country.obj$name
    country <- country.obj$code
  }
  res$results <- summary(coda.list.mcmc(object, country = country, par.names = par.names,
                                        par.names.cs = par.names.cs, thin = thin, 
                                        burnin = burnin), ...)
  return(res)
}


#' 
#' @export
summary.bayesMig.mcmc.meta <- function(object, ...) {
  res <- list(est.period = paste(object$start.year, object$present.year, sep = '-'),
              nr.countries = object$nr.countries,
              nr.countries.est = object$nr.countries.est,
              data.source = if(object$user.data) "user-defined" else "WPP ",
              wpp.year = if(object$user.data) NULL else object$wpp.year
  )
  class(res) <- "summary.bayesMig.mcmc.meta"
  return(res)
}

#' @export
print.summary.bayesMig.mcmc.meta <- function(x, ...) {
  cat('\nNumber of locations:', x$nr.countries)
    if(x$nr.countries != x$nr.countries.est)
        cat('\nNumber of locations influencing world posterior:', x$nr.countries.est)
  cat('\nData source:', x$data.source, x$wpp.year)
  cat('\nInput data: migration for period', x$est.period)
  cat('\n')
}

#' @export
#' 
print.summary.bayesMig.mcmc.set <- function(x, ...) {
  print(x$meta)
  if(!is.null(x$chain.info)) bayesTFR:::print.summary.chain.info(x$chain.info)
  if(!is.null(x$mcmc)) print(x$mcmc)
  if(!is.null(x$country.name)){
    cat('\nLocation:', x$country.name, '\n')
    if (is.null(x$results))
      cat('\tnot used for estimation.\n')
  }
  if(!is.null(x$results)) print(x$results)
}

#' @export
#' 
print.bayesMig.mcmc <- function(x, ...) {
  print(summary(x, ...))
}

#' @export
print.summary.bayesMig.mcmc <- function(x, ...) {
  if(!is.null(x$country.name)){
    cat('\nCountry:', x$country.name, '\n')
    if (is.null(x$results))
      cat('\tnot used for estimation.\n')
  }
  if(!is.null(x$results))
    print(x$results)
}

#' @export
print.bayesMig.mcmc.set <- function(x, ...) {
  print(summary(x, ...))
}

#' @export
print.bayesMig.mcmc.meta <- function(x, ...) {
  print(summary(x, ...))
}

#' @export
print.bayesMig.prediction <- function(x, ...) {
  print(summary(x, ...))
}


#' @title Summary of Prediction of Net Migration Rate
#'
#' @description Country-specific summary of an object of class \code{\link{bayesMig.prediction}}, 
#'     created using the function \code{\link{mig.predict}}. The summary contains the mean, 
#'     standard deviation and several commonly used quantiles of the simulated trajectories.
#' 
#' @param object Object of class \code{\link{bayesMig.prediction}}.
#' @param country Country name or code if a country-specific summary is desired. 
#'     The code can be either numeric or ISO-2 or ISO-3 characters. If it is \code{NULL}, only prediction parameters are included.
#' @param compact Logical switching between a smaller and larger number of displayed quantiles.
#' @param \dots A list of further arguments.
#' @export
#' @rdname summary-prediction
#' 
summary.bayesMig.prediction <- function(object, country = NULL, compact = TRUE, ...) {
  res <- bayesTFR:::get.prediction.summary.data(object, 
                                                unchanged.pars=c('burnin', 'nr.traj'), 
                                                country=country, compact=compact)
  class(res) <- 'summary.bayesMig.prediction'
  return(bayesTFR:::.update.summary.data.by.shift(res, object, country))
}

#' @param x A result of the \code{summary} function.
#' @param digits Minimal number of significant digits.
#' @rdname summary-prediction
#' @export
#' 
print.summary.bayesMig.prediction <- function(x, digits = 3, ...) {
  cat('\nProjections:', length(x$projection.years), '(', x$projection.years[1], '-', 
      x$projection.years[length(x$projection.years)], ')')
  cat('\nTrajectories:', x$nr.traj)
  cat('\nBurnin:', x$burnin)
  
  if(!is.null(x$country.name)) {
    cat('\nCountry:', x$country.name, '\n')
    cat('\nProjected Migration Rate:')
    cat('\n')
    print(x$projections, digits=digits, ...)
  }
}

#' @title Summary of Convergence Diagnostics
#'
#' @description Summary of an object of class \code{\link{bayesMig.convergence}} created 
#'     using the \code{\link{mig.diagnose}} function. It gives an overview about parameters 
#'     that did not converge.
#' 
#' @param object Object of class \code{\link{bayesMig.prediction}}.
#' @param expand By default, the function does not show country-specific parameters 
#'     for which there was no convergence (only country-independent parameters), 
#'     if the status is \sQuote{red}. This argument can switch that option on.
#' @param \dots Not used.
#' @export
#' 

summary.bayesMig.convergence <- function(object, expand=FALSE, ...) {
  return(bayesTFR:::summary.bayesTFR.convergence(object, expand=expand, ...))
}

#' @export
get.nr.countries.bayesMig.mcmc.meta <- function(meta, ...) 
  return(meta$nr.countries)

#' @export
get.nrest.countries.bayesMig.mcmc.meta <- function(meta, ...) 
  return(meta$nr.countries)

country.names.bayesMig.mcmc.meta <- function(meta) {
  return(meta$regions$country_name)
}

#' @export
get.countries.index.bayesMig.mcmc.meta  <- function(meta, ...) 
  return (1:meta$nr.countries)

#' @export
get.countries.table.bayesMig.mcmc.set <- function(object, ...) 
  return(bayesTFR:::get.countries.table.bayesTFR.mcmc.set(object,...))

#' @export
get.countries.table.bayesMig.prediction <- function(object, ...) 
  return(bayesTFR:::get.countries.table.bayesTFR.prediction(object,...))

#' @export
get.data.matrix.bayesMig.mcmc.meta <- function(meta, ...) return (t(meta$mig.rates))

#' @export
get.mcmc.meta.bayesMig.mcmc.set <- function(meta, ...) return(meta$meta)

#' @export
get.mcmc.meta.bayesMig.mcmc.meta <- function(meta, ...) return(meta)

#' @export
get.mcmc.meta.bayesMig.mcmc <- function(meta, ...) return(meta$meta)

#' @export
get.mcmc.meta.list <- function(meta, ...) return(meta[[1]]$meta)

#' @export
get.mcmc.meta.bayesMig.prediction <- function(meta, ...) return(meta$mcmc.set$meta)

