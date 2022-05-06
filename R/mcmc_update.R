###############
#MIGRATION
###############

mcmc.update.mu.c <- function(country, mcmc){
  #Update mu_c conditional on everything else.
  #The posterior is a normal distribution with known mean and variance given by the formulas below.
  
  #Check for constraints first that fix some values.
  if(mcmc$meta$constraints.logical$mu[country]){
    #If there is a constraint, assume we've already got the correct value, and return.
    return()
  }

  bigT=mcmc$meta$bigT
  mig.rates.c=mcmc$meta$mig.rates[country,]
  
  inverseVariance <- (bigT-1)*(1-mcmc$phi_c[country])^2/mcmc$sigma2_c[country] + 1/mcmc$sigma2_mu
  variance <- 1/inverseVariance
  
  mean <- variance*((1-mcmc$phi_c[country])/mcmc$sigma2_c[country]*sum(mig.rates.c[2:bigT]-mcmc$phi_c[country]*mig.rates.c[1:(bigT-1)]) +
                      mcmc$mu_global/mcmc$sigma2_mu)
  #cat("Country",country,"Variance",variance,"\n")
  updated.mu.c <- rnorm(n=1,mean=mean,sd=sqrt(variance))
  mcmc$mu_c[country]=updated.mu.c
  return()
}

mcmc.update.phi.c <- function(country, mcmc){
  #Update phi_c conditional on everything else.
  #The posterior is a truncated normal distribution with known mean and variance.
  
  #Check for constraints first that fix some values.
  if(mcmc$meta$constraints.logical$phi[country]){
    #If there is a constraint, assume we've already got the correct value, and return.
    return()
  }
  
  bigT=mcmc$meta$bigT
  mig.rates.c=mcmc$meta$mig.rates[country,]
  
  #Sample from a truncated normal distribution.
  denominator=sum( (mig.rates.c[1:(bigT-1)] - mcmc$mu_c[country])^2 );
  mean=sum( (mig.rates.c[1:(bigT-1)] - mcmc$mu_c[country]) * (mig.rates.c[2:bigT] - mcmc$mu_c[country]) ) / denominator;
  variance=mcmc$sigma2_c[country] / denominator;
  updated.phi.c=rtruncnorm(n=1,a=0,b=1,mean=mean,sd=sqrt(variance));
  mcmc$phi_c[country]=updated.phi.c
  return();
}

mcmc.update.sigma2.c <- function(country, mcmc){
  #Update sigma^2_c conditional on everything else.
  #The posterior is an inverse gamma distribution.
  
  #Check for constraints first that fix some values.
  if(mcmc$meta$constraints.logical$sigma2[country]){
    #If there is a constraint, assume we've already got the correct value, and return.
    return()
  }
  
  bigT=mcmc$meta$bigT
  mig.rates.c=mcmc$meta$mig.rates[country,]
  
  #The inverse variance (tau_c) follows a gamma distribution
  #Sample a new value for tau_c and return its inverse.
  
  #tempVec is an auxiliary for one of the terms in rate parameter of the gamma distribution
  tempVec=(mig.rates.c[2:bigT] - mcmc$mu_c[country]) - mcmc$phi_c[country] * (mig.rates.c[1:(bigT-1)] - mcmc$mu_c[country]);
  #rate = mcmc$b + 0.5*sum(tempVec^2);
  rate <- sqrt(mcmc$s/mcmc$r) + 0.5*sum(tempVec^2)
  #shape = mcmc$a + (bigT-1)/2;
  shape <- sqrt(mcmc$r * mcmc$s) + (bigT-1)/2;
  tau_c=rgamma(n=1,shape=shape,rate=rate);
  
  #If we end up with sigma^2_c too small, just default to a reasonable small value.
  #Rather than a truncated inverse gamma, this just condenses a bit of the lower tail into a point mass.
  #This is a problem we run into because some countries have only zeroes for migration rate data.
  #With all zeroes, we tend to push towards sigma^2_c=0 (exactly zero), which causes variance problems.
  if(sqrt(1/tau_c)<mcmc$meta$sigma.c.min){
    mcmc$sigma2_c[country]=mcmc$meta$sigma.c.min^2
  }else{
    mcmc$sigma2_c[country]=1/tau_c;    
  }
  
  return()
}

mcmc.update.mu.global <- function(mcmc){
  #Update mu_global conditional on everything else.
  #The posterior is a truncated normal distribution.

  bigC=mcmc$meta$nr.countries.est
  mean=sum(mcmc$mu_c[mcmc$meta$country.indices.est])/bigC
  variance=mcmc$sigma2_mu/bigC
  updated.mu.global=rtruncnorm(n=1,a=mcmc$meta$mu.global.lower,b=mcmc$meta$mu.global.upper,mean=mean,sd=sqrt(variance))
  mcmc$mu_global=updated.mu.global
  return()
}

mcmc.update.sigma2.mu <- function(mcmc){
  #The conditional posterior for sigma^2_mu is a truncated inverse gamma distribution.
  #Try drawing from the appropriate inverse gamma until we get a value within bounds.
  bigC=mcmc$meta$nr.countries.est
  
  #Calculate shape and rate parameters for an inverse gamma distribution
  shape=(bigC-1)/2;rate=sum(mcmc$mu_c[mcmc$meta$country.indices.est]^2)/2;
  
  #Try simulating tau values from a gamma distribution, reject if they fall outside the acceptable region.
  tau.cutoff.lower=mcmc$meta$sigma.mu.upper^(-2)
  tau.cutoff.upper=mcmc$meta$sigma.mu.lower^(-2)
  tau=rgamma(n=1,shape=shape,rate=rate);
  i <- 1
  while(tau<tau.cutoff.lower || tau > tau.cutoff.upper){
    tau=rgamma(n=1,shape=shape,rate=rate)
    if(i > 100) {
      tau <- min(max(tau, tau.cutoff.lower), tau.cutoff.upper)
      break
    }
    i <- i + 1
  }
  mcmc$sigma2_mu=1/tau
  return()
}

mcmc.update.b <- function(mcmc){
  bigC=mcmc$meta$nr.countries.est
  
  #The posterior distribution for b is a truncated gamma.
  shape=mcmc$a*bigC+1;
  rate=sum(1/mcmc$sigma2_c[mcmc$meta$country.indices.est]);
  cutoff=(mcmc$a-1)/2;
  
  #Simple rejection sampling to draw from the truncated gamma distribution.
  b=rgamma(n=1,shape=shape,rate=rate);
  while(b>cutoff){
    b=rgamma(n=1,shape=shape,rate=rate);    
  }
  mcmc$b=b
  return();
}

#This is the old code for updating a. It uses too wide a range of options and gets stuck a lot.
#It also doesn't correctly ignore the small countries.
# mcmc.update.a <- function(mcmc){
#   bigC=length(mcmc$sigma2_c)
#   
#   #Metropolis step with a flat proposal distribution on the possible range of (2b+1, a.upper)
#   proposal=runif(n=1,min=2*mcmc$b+1,max=mcmc$meta$a.upper);
#   
#   sumLogTau=-sum(log(mcmc$sigma2_c));
#   
#   logPosterior_old=bigC*mcmc$a*log(mcmc$b) + (mcmc$a-1)*sumLogTau - bigC*lgamma(mcmc$a) - log(mcmc$a-1)
#   logPosterior_new=bigC*proposal*log(mcmc$b) + (proposal-1)*sumLogTau - bigC*lgamma(proposal) - log(proposal-1)
#   A=logPosterior_new-logPosterior_old
#   
#   #Now the standard Metropolis check. If the proposal is an improvement over the old value, accept it immediately.
#   if(A>0){
#     mcmc$a=proposal
#     return()
#   }
#   #Otherwise, accept with log-probability A
#   U=runif(1);logU=log(U);
#   if(logU < A){
#     mcmc$a=proposal
#     return()
#   }
#   #If we didn't accept (with log-probability A), then don't change the value of a.
#   return()
# }

mcmc.update.a <- function(mcmc){
  bigC=mcmc$meta$nr.countries.est
  aCurrent=mcmc$a
  proposal.half.width=mcmc$meta$a.half.width
  #Metropolis step with a flat proposal distribution on the possible range of (a-proposal.half.width, a+proposal.half.width).
  #Automatically return the current value if we're outside the range (2b+1, a.upper)
  proposal=runif(n=1,min=aCurrent-proposal.half.width,
                 max=aCurrent+proposal.half.width);
  if(proposal<2*mcmc$b+1 || proposal>mcmc$meta$a.upper){
    return();
  }
  #Otherwise, do the usual Metropolis stuff.
  
  sumLogTau=-sum(log(mcmc$sigma2_c[mcmc$meta$country.indices.est]));
  
  logPosterior_old=bigC*mcmc$a*log(mcmc$b) + (mcmc$a-1)*sumLogTau - bigC*lgamma(mcmc$a) - log(mcmc$a-1)
  logPosterior_new=bigC*proposal*log(mcmc$b) + (proposal-1)*sumLogTau - bigC*lgamma(proposal) - log(proposal-1)
  A=logPosterior_new-logPosterior_old
  
  #Now the standard Metropolis check. If the proposal is an improvement over the old value, accept it immediately.
  if(A>0){
    mcmc$a=proposal
    return()
  }
  #Otherwise, accept with log-probability A
  U=runif(1);logU=log(U);
  if(logU < A){
    mcmc$a=proposal
    return()
  }
  #If we didn't accept (with log-probability A), then don't change the value of a.
  return()
}

mcmc.update.r <- function(mcmc){
  bigC=mcmc$meta$nr.countries.est
  current=mcmc$r
  proposal.half.width=mcmc$meta$r.half.width
  #Metropolis step with a flat proposal distribution on the possible range of (a-proposal.half.width, a+proposal.half.width).
  #Automatically return the current value if we're outside the range (2b+1, a.upper)
  min.prop <- max(1/mcmc$s, current-proposal.half.width)
  max.prop <- min(max(current, min.prop)+proposal.half.width, (mcmc$meta$a.upper^2)/mcmc$s)
  #min.prop <- max(min.prop, ((2*sqrt(mcmc$s/max.prop) + 1)^2)/mcmc$s)
  #if(min.prop > max.prop) proposal <- max.prop
  #else {
  #min.prop <- max(0, current-proposal.half.width)
  found <- FALSE
  for (i in 1:100) {
  #proposal=runif(n=1,min=min.prop, max=min(max(current, min.prop)+proposal.half.width, (mcmc$meta$a.upper^2)/mcmc$s))
    proposal=runif(n=1,min=min.prop, max=max.prop)
    if(sqrt(proposal * mcmc$s) >= 2*sqrt(mcmc$s/proposal) + 1) {
      found <- TRUE
      break
    }
  }
  if(!found) return()
  #if(proposal<2*mcmc$b+1 || proposal>mcmc$meta$a.upper){
  #  return();
  #}
  #Otherwise, do the usual Metropolis stuff.
  
  sumLogTau=-sum(log(mcmc$sigma2_c[mcmc$meta$country.indices.est]));
  
  current.a <- sqrt(mcmc$r * mcmc$s)
  current.b <- sqrt(mcmc$s / mcmc$r)
  proposal.a <- sqrt(proposal * mcmc$s)
  proposal.b <- sqrt(mcmc$s / proposal)
  
  # the first term is the Jacobian
  logPosterior_old=-log(2*mcmc$r) + bigC*current.a*log(current.b) + (current.a-1)*sumLogTau - bigC*lgamma(current.a) - log(current.a-1)
  logPosterior_new=-log(2*proposal) + bigC*proposal.a*log(proposal.b) + (proposal.a-1)*sumLogTau - bigC*lgamma(proposal.a) - log(proposal.a-1)

  A=logPosterior_new-logPosterior_old
  
  #Now the standard Metropolis check. If the proposal is an improvement over the old value, accept it immediately.
  if(A>0){
    mcmc$r=proposal
    return()
  }
  #Otherwise, accept with log-probability A
  U=runif(1);logU=log(U);
  if(logU < A){
   mcmc$r=proposal
   return()
  }
  #If we didn't accept (with log-probability A), then don't change the value of r.
  return()
}

mcmc.update.s <- function(mcmc){
  bigC=mcmc$meta$nr.countries.est
  current=mcmc$s
  proposal.half.width=mcmc$meta$s.half.width
  #Metropolis step with a flat proposal distribution on the possible range of (a-proposal.half.width, a+proposal.half.width).
  #Automatically return the current value if we're outside the range (2b+1, a.upper)
  #min.prop <- max(mcmc$r/(mcmc$r - 2)^2, current-proposal.half.width)
  min.prop <- max(1/mcmc$r, current-proposal.half.width)
  #min.prop <- max(min.prop, (2*sqrt(min.prop/mcmc$r) + 1)^2/mcmc$r)
  max.prop <- min(max(current, min.prop)+proposal.half.width, (mcmc$meta$a.upper^2)/mcmc$r)
  #if(min.prop > max.prop) proposal <- max.prop
  #else {
  found <- FALSE
  for(i in 1:100) {
    proposal=runif(n=1,min=min.prop, max=max.prop)
    if(sqrt(proposal * mcmc$r) >= 2*sqrt(proposal/mcmc$r) + 1) {
      found <- TRUE
      break
    }
  }
  if(!found) return()
  #max.prop <- min(current+proposal.half.width, 1/mcmc$r)
  #proposal=runif(n=1,min=max(0, min(current, max.prop)-proposal.half.width), max=max.prop)
  #if(proposal<2*mcmc$b+1 || proposal>mcmc$meta$a.upper){
  #  return();
  #}
  #Otherwise, do the usual Metropolis stuff.
  
  sumLogTau=-sum(log(mcmc$sigma2_c[mcmc$meta$country.indices.est]));
  
  current.a <- sqrt(mcmc$r * mcmc$s)
  current.b <- sqrt(mcmc$s / mcmc$r)
  proposal.a <- sqrt(proposal * mcmc$r)
  proposal.b <- sqrt(proposal / mcmc$r)
  
  # bigT=mcmc$meta$bigT
  # mu.logterm <- 0
  # for(country in 1:bigC){
  #   mu.logterm <- mu.logterm + (-(bigT-1)*(1-mcmc$phi_c[country])^2/mcmc$sigma2_c[country] + 1/(2*mcmc$sigma2_mu)) * mcmc$mu_c[country]^2 +
  #     ((1-mcmc$phi_c[country])/mcmc$sigma2_c[country]*sum(mig.rates.c[2:bigT]-mcmc$phi_c[country]*mig.rates.c[1:(bigT-1)]))*mcmc$mu_c[country]
  # }
  # ab.term.current <- bigC*current.a*log(current.b) + (current.a-1)*sumLogTau - bigC*lgamma(current.a) - log(current.a-1)
  # ab.term.proposal <- bigC*proposal.a*log(proposal.b) + (proposal.a-1)*sumLogTau - bigC*lgamma(proposal.a) - log(proposal.a-1)
  
  # Jacobian no needed here as it is the same for old and new
  logPosterior_old = bigC*current.a*log(current.b) + (current.a-1)*sumLogTau - bigC*lgamma(current.a) - log(current.a-1)
  logPosterior_new = bigC*proposal.a*log(proposal.b) + (proposal.a-1)*sumLogTau - bigC*lgamma(proposal.a) - log(proposal.a-1)
  A=logPosterior_new-logPosterior_old
  
  #Now the standard Metropolis check. If the proposal is an improvement over the old value, accept it immediately.
  if(A>0){
    mcmc$s=proposal
    return()
  }
  #Otherwise, accept with log-probability A
  U=runif(1);logU=log(U);
  if(logU < A){
   mcmc$s=proposal
   return()
  }
  #If we didn't accept (with log-probability A), then don't change the value of r.
  return()
}

