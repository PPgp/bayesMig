###############
#MIGRATION
###############

run.mig.mcmc <- function(nr.chains=3, iter=1000, 
                         thin=1, verbose=TRUE, verbose.iter=10,
                         # meta parameters
                         #                         start.year=1950, present.year=2010, wpp.year=2012,
                         # starting values (length of 1 or nr.chains)
                         #                         S.ini=NULL, a.ini=NULL, b.ini=NULL, 
                         #                         verbose=FALSE, verbose.iter = 10, 
                         output.dir=file.path(getwd(), 'bayesMig.output'), replace.output=FALSE,
                         parallel=FALSE, ...){
  
  if(file.exists(output.dir)) {
    if(length(list.files(output.dir)) > 0 & !replace.output)
      stop('Non-empty directory ', output.dir, 
           ' already exists.\nSet replace.output=TRUE if you want to overwrite existing results.')
    unlink(output.dir, recursive=TRUE)
  }
  dir.create(output.dir)
  
  #Auto run stuff goes here.
  
  if (verbose) {
    cat('\nStarting Bayesian Hierarchical Model Migration.\n')
    cat('========================================================\n')
    cat('Initialize simulation -', nr.chains, 'chain(s) in total.\n')
  }
  
  #  if(!is.null(seed)) set.seed(seed)
  
  #A bunch of initializations, which should be fed into some initialization later.
  # starting values (length of 1 or nr.chains)
  #  if(missing(mu_c.ini) || is.null(mu_c.ini)){
  #    mu_c.ini=something
  #  }
  
  bayesMig.mcmc.meta = mcmc.meta.ini(output.dir=output.dir)
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
    cat('(Not implemented yet.)\n')
    #Print starting values here.
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



##############
#Test code here
##############

#source("mcmc_ini.R")
#source("mcmc_update.R")
#source("mcmc_sampling.R")
#source("mcmc_storage.R")
#source("get_outputs.R")

#Run some chains
#run.mig.mcmc(replace.output=TRUE)
#run.mig.mcmc(replace.output=TRUE, thin=10, iter=1000, verbose.iter=100,
#             output.dir=file.path(getwd(), 'bayesMig.output'))
