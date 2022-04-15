###############
#MIGRATION
###############

#' @title Run Markov chain Monte Carlo for parameters of net migration rate model
#'
#' @description Runs MCMCs for simulating the net migration rate of all countries of the
#' world (or locations specified by users), using the Bayesian hierarchical model of Azose & Raftery (2015).
#' 
#' @param nr.chains An integer number of independent Markov chains to run.
#' @param iter The number of iterations to run per Markov chain.
#' @param thin Thinning interval -- A chain with 1000 iterations thinned by 20 will return a 
#' final count of 50 iterations.
#' @param start.year Start year for using historical data.
#' @param present.year End year for using historical data.
#' @param wpp.year Year for which WPP data is used if no user data is provided via \code{my.mig.file}. 
#' In such a case, the function loads a package called \pkg{wpp}\eqn{x} where \eqn{x} is the \code{wpp.year} and generates 
#' historical migration rates using the 
#' \code{\link[wpp2019]{migration}} and \code{\link[wpp2019]{pop}} datasets.
#' @param my.mig.file File name containing user-specified historical time series of migration rates 
#' for all locations that should be included in the simulation. 
#' @param sigma.c.min,a.up,a.ini,a.half.width,mu.range,sigma.mu.range,mu.ini Settings for the parameters
#' of the model (see Azose & Raftery 2015), such as minimum value, truncation ranges, slicing half width and initial values.
#' Initial values (*.ini) can be given as a vector of length \code{nr.chains}, giving one initial value per chain.
#' By default the initial values are equidistantly spread between their respective ranges.
#' @param exclude.from.world Vector of country codes that should not influence the hyperparameters. 
#' However, country-specific parameters will be generated for these countries.
#' @param seed Seed of the random number generator. If \code{NULL} no seed is set. It can be used to generate reproducible results.
#' @param verbose Whether or not to print status updates to console window while code is running.
#' @param verbose.iter If verbose is TRUE, the number of iterations to wait between printing updates.
#' @param output.dir A file path pointing to the directory in which to store results.
#' @param replace.output If the specified output directory already exists, should it be overwritten?
#' @param annual If \code{TRUE}, the model assumes the underlying data is on annual time scale. 
#'     In such a case, argument \code{my.mig.file} must be used to provide the annual observed data.
#' @param parallel Whether to run code in parallel.
#' @param nr.nodes Relevant only if \code{parallel} is \code{TRUE}. It gives the number of nodes for running the simulation in parallel. 
#' By default it equals to the number of chains.
#' @param buffer.size Buffer size (in number of iterations) for keeping data in the memory before flushing to disk.
#' @param \dots Additional parameters to be passed to the function \code{\link[snowFT]{performParallel}}, if \code{parallel} is \code{TRUE}.
#' 
#' @return An object of class \code{bayesMig.mcmc.set} which is a list with two components:
#' \item{meta}{An object of class \code{\link{bayesMig.mcmc.meta}}.}
#' \item{mcmc.list}{A list of objects of class \code{\link{bayesMig.mcmc}}, one for each MCMC.}
#' 
#' @details The function creates an object of class \code{\link{bayesMig.mcmc.meta}} and 
#' stores it in \code{output.dir}. It launches \code{nr.chains} MCMCs, either sequentially or 
#' in parallel. Parameter traces of each chain are stored as ASCII files in a subdirectory 
#' of \code{output.dir}, called \code{mc}\emph{x} where \emph{x} is the identifier of that chain. 
#' There is one file per parameter, named after the parameter with the suffix \dQuote{.txt}.
#' Country-specific parameters have the suffix \code{_country}\emph{c} where \emph{c} is the country code.
#' In addition to the trace files, each \code{mc}\emph{x} directory contains the object 
#' \code{\link{bayesMig.mcmc}} in binary format.  
#' All chain-specific files  are written onto disk after the first, last and each 
#' \eqn{i}-th (thinned) iteration, where \eqn{i} is given by the argument \code{buffer.size}.
#' 
#' By default (if no data is passed via the \code{my.mig.file} argument), the function 
#' loads observed data (further denoted as WPP dataset), from the \code{\link[wpp2019]{migration}} 
#' and \code{\link[wpp2019]{pop}} datasets in the \pkg{wpp}\eqn{x} package where \eqn{x} is 
#' the \code{wpp.year}. Net migration rates are computed as migration(\eqn{t}) / (population(\eqn{t_e}) - migration(\eqn{t})) 
#' where \eqn{t_e} means the end of time period \eqn{t}. 
#' 
#' The argument \code{my.mig.file} can be used to overwrite the default data. 
#' If it is used, it should contain the net migration rate for all locations to be used in the simulation, as no WPP data is used 
#' in such a case. The structure of the file has the same format as the \code{\link[wpp2019]{migration}} dataset,
#' but the values should be rates (instead of counts). Each row corresponds to a location. It does not have 
#' to be necessarily a country - it can be for example a subnational unit. It must contain columns 
#' \dQuote{country_code} or \dQuote{code} (unique identifier of the location), \dQuote{name}, and columns representing 
#' 5-year time intervals (if \code{annual} is \code{FALSE}), e.g., \dQuote{1995-2000}, \dQuote{2000-2005} etc. 
#' 
#' If \code{annual} is \code{TRUE} the default WPP dataset is not used and the \code{my.mig.file} argument 
#' must provide the dataset to be used for estimation. Its time-related columns should be single years.
#' 
#' If there are countries or locations that should be excluded from influencing the hyperparameters,
#' for example small countries or locations with unique migration patterns, their codes 
#' should be included in the argument \code{exclude.from.world}. These locations will still get 
#' their parameters simulated and thus, can be included in a projection.
#' 
#' @aliases bayesMig.mcmc.set
#' 
#' @references Azose, J. J., & Raftery, A. E. (2015). 
#' Bayesian probabilistic projection of international migration. Demography, 52(5), 1627-1650.
#' 
#' @seealso \code{\link{get.mig.mcmc}}, \code{\link{summary.bayesMig.mcmc.set}}, \code{\link{mig.partraces.plot}},
#' \code{\link{mig.pardensity.plot}}, \code{\link{mig.predict}}
#' 
#' @examples
#' m <- run.mig.mcmc(nr.chains = 2, iter = 30, thin = 1)
#' summary(m)
#' 
#' @export
run.mig.mcmc <- function(nr.chains=3, iter=50000, output.dir=file.path(getwd(), 'bayesMig.output'), 
                         thin=1, replace.output=FALSE, annual = FALSE,
                         start.year = 1950, present.year=2020, wpp.year=2019, my.mig.file = NULL,
                         # starting values and ranges for truncations
                         sigma.c.min = 0.0001, a.up = 10, a.ini = NULL, a.half.width = 0.3,
                         mu.range = c(-0.5, 0.5), sigma.mu.range = c(0, 0.5), mu.ini = NULL,
                         # other settings
                         exclude.from.world = NULL,
                         seed = NULL, parallel = FALSE, nr.nodes = nr.chains, 
                         buffer.size = 1000, verbose=TRUE, verbose.iter=10, ...){
  
  if(file.exists(output.dir)) {
    if(length(list.files(output.dir)) > 0 & !replace.output)
      stop('Non-empty directory ', output.dir, 
           ' already exists.\nSet replace.output=TRUE if you want to overwrite existing results.')
    unlink(output.dir, recursive=TRUE)
  }
  dir.create(output.dir)
  
  #Auto run stuff goes here.
  
  init.values.between.low.and.up <- function(low, up)
    ifelse(rep(nr.chains==1, nr.chains), (low + up)/2, seq(low, to=up, length=nr.chains))
  
  
  if (verbose) {
    cat('\nStarting Bayesian Hierarchical Model Migration.\n')
    cat('========================================================\n')
    cat('Initialize simulation -', nr.chains, 'chain(s) in total.\n')
  }
  
  if(!is.null(seed)) set.seed(seed)
  
  #A bunch of initializations, which should be fed into some initialization later.
  # starting values (length of 1 or nr.chains)
  if(missing(mu.ini) || is.null(mu.ini)){
      mu.ini <- init.values.between.low.and.up(mu.range[1], mu.range[2])
  }
  if(missing(a.ini) || is.null(a.ini)){
    a.ini <- init.values.between.low.and.up(1.001, 5)
  }
  
  bayesMig.mcmc.meta <- mcmc.meta.ini(output.dir=output.dir, wpp.year = wpp.year,
                                      start.year=start.year, present.year = present.year, 
                                      annual.simulation = annual,
                                      my.mig.file = my.mig.file, 
                                     sigma.c.min = sigma.c.min, a.up = a.up,
                                     mu.range = mu.range, sigma.mu.range = sigma.mu.range,
                                     mu.ini = mu.ini, a.ini = a.ini, a.half.width = a.half.width,
                                     exclude.from.world = exclude.from.world, buffer.size = buffer.size)
  #cat(bayesMig.mcmc.meta$mig.rates)
  
  #Storage  
  store.bayesMig.meta.object(bayesMig.mcmc.meta, output.dir)
  
  # propagate initial values for all chains if needed

  if (parallel) { # run chains in parallel
    chain.set <- bayesTFR:::bDem.performParallel(nr.nodes, 1:nr.chains, mcmc.run.chain.mig, 
                                      initfun = mig.init.nodes, seed = seed,
                                      meta = bayesMig.mcmc.meta, 
                                      thin = thin, iter = iter, verbose = verbose, 
                                      verbose.iter = verbose.iter, ...)
  } else { # run chains sequentially
    chain.set <- list()
    for (chain in 1:nr.chains) {
      chain.set[[chain]] <- mcmc.run.chain.mig(chain, bayesMig.mcmc.meta, thin = thin, 
                                           iter = iter, verbose = verbose, verbose.iter = verbose.iter)
    }
  }
  names(chain.set) <- 1:nr.chains
  mcmc.set <- structure(list(meta=bayesMig.mcmc.meta, mcmc.list=chain.set), class='bayesMig.mcmc.set')
  cat('\nResults stored in', output.dir,'\n')
  
  if (verbose)
    cat('\nSimulation successfully finished!!!\n')
  invisible(mcmc.set)
}

mig.init.nodes <- function(){library(bayesMig)}

mcmc.run.chain.mig <- function(chain.id, meta, thin=1, iter=100,
                        #In the final version, probably also take initial values
                        verbose=FALSE, verbose.iter=10){
  cat('\n\nChain nr.', chain.id, '\n')
  if (verbose) {
    cat('************\n')
    cat('Starting values:\n')
    sv <- c(meta$mu.ini[chain.id], meta$a.ini[chain.id])
    names(sv) <- c('mu', 'a')
    print(sv)
  }
  
  #This will eventually get updated with a lot more parameters.
  mcmc=mcmc.ini(chain.id=chain.id, mcmc.meta=meta, iter=iter)  
  
  if (verbose){
    cat('Store initial values into ', mcmc$output.dir, '\n')
  } 
  
  #Storage goes here.
  store.mcmc(mcmc, append=FALSE, flush.buffer=TRUE, verbose=verbose)
  
  if (verbose) 
    cat('Start sampling -', mcmc$iter, 'iterations in total.\n')
  mcmc <- mig.mcmc.sampling(mcmc, thin=thin, verbose=verbose, verbose.iter=verbose.iter)
  return(mcmc)
}

