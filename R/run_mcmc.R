###############
#MIGRATION
###############

#' @title Run Markov chain Monte Carlo for parameters of net migration rate model
#'
#' @description Runs MCMCs for simulating the net migration rate of all countries of the
#' world or for locations specified by users, using the Bayesian hierarchical model of Azose & Raftery (2015).
#' 
#' @param nr.chains An integer number of independent Markov chains to run.
#' @param iter The number of iterations to run per Markov chain.
#' @param thin Thinning interval. A chain with 1000 iterations thinned by 20 will return a 
#' final count of 50 iterations.
#' @param start.year Start year for using historical data.
#' @param present.year End year for using historical data.
#' @param wpp.year Year for which WPP data is used if no user data is provided via \code{my.mig.file}. 
#' In such a case, the function loads a package called \pkg{wpp}\eqn{x} where \eqn{x} is the \code{wpp.year} and generates 
#' historical migration rates using the 
#' \code{\link[wpp2019]{migration}} and \code{\link[wpp2019]{pop}} datasets.
#' @param my.mig.file File name containing user-specified historical time series of migration rates 
#'      for all locations that should be included in the simulation. It should be a tab-separated file.
#'      For structure, see Details below.
#' @param sigma.c.min,a.ini,mu.ini Settings for the parameters
#' of the model (see Azose & Raftery 2015), such as minimum value and initial values.
#' Initial values (*.ini) can be given as a vector of length \code{nr.chains}, giving one initial value per chain.
#' By default the initial values are equidistantly spread between their respective ranges.
#' @param a.half.width Half width for Metropolis proposals of the a parameter. This argument can greatly influence 
#'      the convergence and it is dependent on the scale of the data. By default it is set to 0.01 for 5-year data 
#'      defined as rate per population; to 0.03 for 5-year data defined as per 1000; to 0.3 for 
#'      annual data per population; to 0.5 for annual data per 1000. If the default does not 
#'      yield satisfactory results, use the function \code{\link{estimate.a.hw}} to estimate 
#'      an appropriate value, based on an existing simulation. Also it is important to set the \code{pop.denom}
#'      argument correctly. 
#' @param exclude.from.world Vector of location codes that should be excluded from estimating the hyperparameters. 
#'      These would be for example small locations or locations with unusual patters. 
#'      Note that location-specific parameters are generated for all locations, regardless of this setting.
#' @param pop.denom Denominator used to generate the input migration rates. It is used to derive an appropriate scaler 
#'      for the priors and conditional distributions. Typically, this will be either 1 (default) if the rates are 
#'      defined as per population, or 1000, if the rates are per 1000 population. 
#'      Use this argument only if user-specified rates are supplied via the \code{my.mig.file} argument. 
#' @param seed Seed of the random number generator. If \code{NULL} no seed is set. It can be used to generate reproducible results.
#' @param verbose Whether or not to print status updates to console window while the code is running.
#' @param verbose.iter If verbose is TRUE, the number of iterations to wait between printing updates.
#' @param output.dir A file path pointing to the directory in which to store results.
#' @param replace.output If the specified output directory already exists, should it be overwritten?
#' @param annual If \code{TRUE}, the model assumes the underlying data is on annual time scale. 
#' @param parallel Whether to run code in parallel.
#' @param nr.nodes Relevant only if \code{parallel} is \code{TRUE}. It gives the number of nodes for running the simulation in parallel. 
#' By default it equals to the number of chains.
#' @param buffer.size Buffer size (in number of iterations) for keeping data in the memory before flushing to disk.
#' @param \dots Additional parameters to be passed to the function \code{\link[snowFT]{performParallel}}, if \code{parallel} is \code{TRUE}.
#' 
#' @return An object of class \code{bayesMig.mcmc.set} which is a list with two components:
#' \item{meta}{An object of class \code{bayesMig.mcmc.meta}. It contains information that is common to all chains.
#'     Most items are the same as in \code{\link[bayesTFR]{bayesTFR.mcmc.meta}}. In addition, \code{mig.rates}
#'     is a matrix of the observed migration rates with \code{NA}s in spots that were not used 
#'     for estimation. \code{mig.rates.all} is a similar matrix but contains all data, regardless
#'     if used for estimation or not. Item \code{user.data} is a logical indicating 
#'     if the migration rates are given by the user (\code{TRUE}) or are taken from a \pkg{wpp} package
#'     (\code{FALSE}).}
#' \item{mcmc.list}{A list of objects of class \code{bayesMig.mcmc}, one for each MCMC. 
#'     Information stored here is specific to each MCMC chain, similarly to \code{\link[bayesTFR]{bayesTFR.mcmc}}.}
#' 
#' @details The function creates an object of class \code{\link{bayesMig.mcmc.meta}} and 
#' stores it in \code{output.dir}. It launches \code{nr.chains} MCMCs, either sequentially or 
#' in parallel. Parameter traces of each chain are stored as ASCII files in a subdirectory 
#' of \code{output.dir}, called \code{mc}\emph{x} where \emph{x} is the identifier of that chain. 
#' There is one file per parameter, named after the parameter with the suffix \dQuote{.txt}.
#' Location-specific parameters have the suffix \code{_country}\emph{c} where \emph{c} is the location code.
#' In addition to the trace files, each \code{mc}\emph{x} directory contains the object 
#' \code{\link{bayesMig.mcmc}} in binary format.  
#' All chain-specific files  are written onto disk after the first, last and each 
#' \eqn{i}-th (thinned) iteration, where \eqn{i} is given by the argument \code{buffer.size}.
#' 
#' By default (if no data is passed via the \code{my.mig.file} argument), the function 
#' loads observed data (further denoted as WPP dataset), from the \code{\link[wpp2019]{migration}} 
#' and \code{\link[wpp2019]{pop}} datasets in the \pkg{wpp}\eqn{x} package where \eqn{x} is 
#' the \code{wpp.year}. Net migration rates are computed as migration(\eqn{t}) / (population(\eqn{t_e}) - migration(\eqn{t})) 
#' where \eqn{t_e} means the end of time period \eqn{t}. For an annual simulation and 
#' \code{wpp.year} set to 2022, \eqn{t = t_e} because the population in year \eqn{t} 
#' is considered at the end of the year. If \code{wpp.year} is smaller than 2022 and \code{annual} is \code{TRUE}
#' the default dataset is interpolated from 5-year data.
#' 
#' The argument \code{my.mig.file} can be used to overwrite the default data. It should be a tab-separated file.
#' If it is used, it should contain net migration rates for all locations to be used in the simulation, as no WPP data is used 
#' in such a case. The structure of the file has the same format as the \code{\link[wpp2019]{migration}} dataset,
#' but the values should be rates (instead of counts). Use the argument \code{pop.denom} to define the scale of the 
#' denominator in these rates, i.e. if the rates are to be interpreted as per population (default) or some other scale. 
#' Each row in the \code{my.mig.file} file corresponds to a location. It does not have 
#' to be necessarily a country - it can be for example a subnational unit. It must contain columns 
#' \dQuote{country_code} or \dQuote{code} (unique identifier of the location), \dQuote{name}, and columns representing 
#' 5-year time intervals (if \code{annual} is \code{FALSE}), e.g., \dQuote{1995-2000}, \dQuote{2000-2005} etc., or single years
#' (if \code{annual} is \code{TRUE}). An example dataset of annual net migration rates for US states is included in the package, 
#' see example below. 
#' 
#' Optionally, the \code{my.mig.file} can contain columns called \dQuote{first.observed} and/or \dQuote{last.observed}, containing 
#' for each location the year of the first and last observation, respectively. In such a case, any data 
#' before and after those time points will be ignored. Furthermore, the function \code{\link{mig.predict}} fills in the missing values 
#' after the last observation, using the median of the BHM procedure.
#' 
#' If there are countries or locations that should be excluded from influencing the hyperparameters,
#' for example small countries or locations with unique migration patterns, their codes 
#' should be included in the argument \code{exclude.from.world}. These locations will still get 
#' their parameters simulated and thus, will be included in a projection. Alternatively 
#' if \code{my.mig.file} is used, these locations can be determined using an additional column, called 
#' \dQuote{include_code}. Value 2 means the location is included in the BHM; value 1 means it's 
#' excluded but location-specific parameters are generated; value 0 means the location is ignored. 
#' 
#' @aliases bayesMig.mcmc.set bayesMig.mcmc bayesMig.mcmc.meta
#' 
#' @references Azose, J. J., & Raftery, A. E. (2015). 
#' Bayesian probabilistic projection of international migration. Demography, 52(5), 1627-1650.
#' 
#' @seealso \code{\link{get.mig.mcmc}}, \code{\link{summary.bayesMig.mcmc.set}}, \code{\link{mig.partraces.plot}},
#' \code{\link{mig.pardensity.plot}}, \code{\link{mig.predict}}
#' 
#' @examples
#' \dontrun{
#' # Toy simulation for US states
#' us.mig.file <- file.path(find.package("bayesMig"), "extdata", "USmigrates.txt")
#' sim.dir <- tempfile()
#' m <- run.mig.mcmc(nr.chains = 3, iter = 100, thin = 1, my.mig.file = us.mig.file, 
#'         annual = TRUE, output.dir = sim.dir)
#' summary(m)
#' summary(m, "Washington")
#' 
#' mig.partraces.plot(m)
#' mig.partraces.cs.plot("California", m)
#' 
#' # later one can access the object from disk
#' m <- get.mig.mcmc(sim.dir)
#'  
#' unlink(sim.dir, recursive = TRUE)
#' # For a country-level simulation, see example in ?bayesMig. 
#' }
#' @export
#' 
run.mig.mcmc <- function(nr.chains=3, iter=50000, output.dir=file.path(getwd(), 'bayesMig.output'), 
                         thin=1, replace.output=FALSE, annual = FALSE,
                         start.year = 1950, present.year=2020, wpp.year=2019, my.mig.file = NULL,
                         # starting values and ranges for truncations
                         sigma.c.min = 0.0001, #a.up = 10, 
                         a.ini = NULL, a.half.width = NULL,
                         #mu.range = c(-0.5, 0.5), sigma.mu.range = c(0, 0.5), 
                         mu.ini = NULL,
                         # other settings
                         exclude.from.world = NULL, pop.denom = 1,
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
  prior.scaler <- (if(annual) 1 else 5) * pop.denom
  mu.range <- c(-1/10 * prior.scaler, 1/10 * prior.scaler)
  sigma.mu.range <- c(0, 1/10 * prior.scaler)
  a.up <- 10
  if(is.null(a.half.width)) {
    if(pop.denom == 1) {
      if (!annual) a.half.width <- 0.01 # 5-year per population
      else a.half.width <- 0.3 # annual per population
    } else {
      if (!annual) a.half.width <- 0.03 # 5-year per 1000
      else a.half.width <- 0.5 # annual per thousand
    }
  }
  
  #A bunch of initializations, which should be fed into some initialization later.
  # starting values (length of 1 or nr.chains)
  if(missing(mu.ini) || is.null(mu.ini)){
      mu.ini <- init.values.between.low.and.up(mu.range[1], mu.range[2])
  }
  if(missing(a.ini) || is.null(a.ini)){
    a.ini <- init.values.between.low.and.up(1.1, a.up/2)
  }
  
  bayesMig.mcmc.meta <- mcmc.meta.ini(nr.chains=nr.chains,
                                    output.dir=output.dir, wpp.year = wpp.year,
                                      start.year=start.year, present.year = present.year, 
                                      annual.simulation = annual,
                                      my.mig.file = my.mig.file, 
                                     sigma.c.min = sigma.c.min, a.up = a.up,
                                     mu.range = mu.range, sigma.mu.range = sigma.mu.range,
                                     mu.ini = mu.ini, a.ini = a.ini, a.half.width = a.half.width,
                                     prior.scaler = prior.scaler, exclude.from.world = exclude.from.world, 
                                    buffer.size = buffer.size, verbose=verbose)
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

