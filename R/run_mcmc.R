###############
#MIGRATION
###############

#' @title Run Markov chain Monte Carlo for parameters of net migration rate model
#'
#' @description Runs (or continues running) MCMCs for simulating the net migration rate of all countries of the
#' world, using a Bayesian hierarchical model.
#' 
#' @param nr.chains An integer number of independent Markov chains to run.
#' @param iter The number of iterations to run per Markov chain.
#' @param thin Thinning interval -- A chain with 1000 iterations thinned by 20 will return a 
#' final count of 50 iterations.
#' @param verbose Whether or not to print status updates to console window while code is running.
#' @param verbose.iter If verbose is TRUE, the number of iterations to wait between printing updates.
#' @param output.dir A file path pointing to the directory in which to store results.
#' @param replace.output If the specified output directory already exists, should it be overwritten?
#' @param parallel (NOT CURRENTLY IMPLEMENTED) Whether to run code in parallel, if available.
#' @param ... Other arguments passed to \code{\link{run.mig.mcmc}}
#' @return An object of class \code{bayesMig.mcmc.set} containing the sampled posterior parameter values
#' @examples
#' run.mig.mcmc(nr.chains=2, iter=30, thin=1)
#' @export
run.mig.mcmc <- function(nr.chains=3, iter=50000, output.dir=file.path(getwd(), 'bayesMig.output'), 
                         thin=1, replace.output=FALSE, 
                         start.year = 1950, present.year=2020, wpp.year=2019, my.mig.file = NULL,
                         # starting values and ranges for truncations
                         sigma.c.min = 0.0001, a.up = 10, a.ini = NULL, a.half.width = 0.3,
                         mu.range = c(-0.5, 0.5), sigma.mu.range = c(0, 0.5), mu.ini = NULL,
                         # other settings
                         seed = NULL, parallel=FALSE, nr.nodes=nr.chains, 
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
                                      my.mig.file = my.mig.file, 
                                     sigma.c.min = sigma.c.min, a.up = a.up,
                                     mu.range = mu.range, sigma.mu.range = sigma.mu.range,
                                     mu.ini = mu.ini, a.ini = a.ini, a.half.width = a.half.width,
                                     buffer.size = buffer.size)
  #cat(bayesMig.mcmc.meta$mig.rates)
  
  #Storage  
  store.bayesMig.meta.object(bayesMig.mcmc.meta, output.dir)
  
  # propagate initial values for all chains if needed

  if (parallel) { # run chains in parallel
#    chain.set <- bDem.performParallel(nr.nodes, 1:nr.chains, mcmc.run.chain, 
#                                      initfun=init.nodes, meta=bayesTFR.mcmc.meta, 
#                                      thin=thin, iter=iter, S.ini=S.ini, a.ini=a.ini,
#                                      b.ini=b.ini, sigma0.ini=sigma0.ini, Triangle_c4.ini=Triangle_c4.ini, const.ini=const.ini,
#                                      gamma.ini=gamma.ini, save.all.parameters=save.all.parameters, verbose=verbose, 
#                                      verbose.iter=verbose.iter, ...)
  } else { # run chains sequentially
    chain.set <- list()
    for (chain in 1:nr.chains) {
      chain.set[[chain]] <- mcmc.run.chain(chain, bayesMig.mcmc.meta, thin=thin, 
                                           iter=iter,
                                           verbose=verbose, verbose.iter=verbose.iter)#FIX THIS!!
    }
  }
  names(chain.set) <- 1:nr.chains
  mcmc.set <- structure(list(meta=bayesMig.mcmc.meta, mcmc.list=chain.set), class='bayesMig.mcmc.set')
  cat('\nResults stored in', output.dir,'\n')
  
#if(auto.run) features go here.

  if (verbose)
    cat('\nSimulation successfully finished!!!\n')
  invisible(mcmc.set)
}

mcmc.run.chain=function(chain.id, meta, thin=1, iter=100,
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

