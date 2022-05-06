###############
#MIGRATION
###############

mig.mcmc.sampling <- function(mcmc, thin=1, start.iter=2, verbose=FALSE, verbose.iter=10) {
  #Draw a sample of length mcmc$iter from all parameter values.
  
  #Look at mcmc_sampling in bayesTFR. A bunch of initializations happen here.
  nSimulations=mcmc$iter
  nC=length(mcmc$mu_c)
  
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
    #Update all the mu_c values with Gibbs sampling.
    for(c in 1:nC){
      mcmc.update.mu.c(c,mcenv)
    }
    #Update sigma^2_mu with Gibbs sampling
    mcmc.update.sigma2.mu(mcenv)
    #Update all the phi_c values with Gibbs sampling.
    for(c in 1:nC){
      mcmc.update.phi.c(c,mcenv)
    }

    #Update all the sigma^2_c values with Gibbs sampling.
    for(c in 1:nC){
      mcmc.update.sigma2.c(c,mcenv)
    }
    
    #Update b with Gibbs sampling.
    #mcmc.update.b(mcenv)
    
    #Update a with Metropolis.
    #mcmc.update.a(mcenv)
    mcmc.update.r(mcenv)
    mcmc.update.s(mcenv)
    
    #Update mu_global with Gibbs sampling
    mcmc.update.mu.global(mcenv)
    
    mcenv$a <- sqrt(mcenv$r*mcenv$s)
    mcenv$b <- sqrt(mcenv$s/mcenv$r)
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



