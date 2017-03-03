###############
#MIGRATION
###############

mig.mcmc.sampling <- function(mcmc, thin=1, start.iter=2, verbose=FALSE, verbose.iter=10) {
  #Draw a sample of length mcmc$nSimulations from all parameter values.
  
  #Look at mcmc_sampling in bayesTFR. A bunch of initializations happen here.
  nSimulations=mcmc$iter
  bigC=length(mcmc$mu_c)
  
  mcenv <- as.environment(mcmc) # Create an environment for the mcmc stuff in order to avoid 
  # copying of the mcmc list 
  
  mcenv$thin <- thin
  
  ################################################################### 
  # Start MCMC
  ############
  for (simu in start.iter:nSimulations) {
    if(verbose.iter > 0 && (simu %% verbose.iter == 0)){
      cat('\nIteration:', simu, '--', date(),"\n")      
    }

    unblock.gtk('bDem.Migmcmc')
    #Something about unblock.gtk here. (Necessary GUI stuff for bayesDem?)
    #Update all the mu_c values with Gibbs sampling.
    for(c in 1:bigC){
      mcmc.update.mu.c(c,mcenv)
    }
    #Update sigma^2_mu with Gibbs sampling
    mcmc.update.sigma2.mu(mcenv)
    #Update all the phi_c values with Gibbs sampling.
    for(c in 1:bigC){
      mcmc.update.phi.c(c,mcenv)
    }

    #Update all the sigma^2_c values with Gibbs sampling.
    for(c in 1:bigC){
      mcmc.update.sigma2.c(c,mcenv)
    }
    
    #Update b with Gibbs sampling.
    mcmc.update.b(mcenv)
    
    #Update a with Metropolis.
    mcmc.update.a(mcenv)
    
    #Update mu_global with Gibbs sampling
    mcmc.update.mu.global(mcenv)
    
    ################################################################### 
    # write samples simu/thin to disk
    ##################################################################
    mcenv$finished.iter <- simu
    mcenv$rng.state <- .Random.seed
    if (simu %% thin == 0){
      mcenv$length <- mcenv$length + 1
      flush.buffer <- FALSE
      if (simu + thin > nSimulations) flush.buffer <- TRUE
      store.mcmc(mcenv, append=TRUE, flush.buffer=flush.buffer, verbose=verbose)
    }
  }       # end simu loop MCMC

  #.cleanup.mcmc(mcenv) #Cleanup step goes here.
  resmc <- as.list(mcenv)
  class(resmc) <- class(mcmc)
  return(resmc)
}

unblock.gtk <- function(option, options.list=NULL) {
  if(!getOption(option, default=FALSE)) return()
  if(!is.null(options.list)) options(options.list)
  # This is to unblock the GUI, if the run is invoked from bayesDem
  # In such a case the gtk libraries are already loaded
  while(do.call('gtkEventsPending', list()))
    do.call('gtkMainIteration', list())
  
}


