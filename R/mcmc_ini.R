###############
#MIGRATION
###############



mcmc.meta.ini <- function(...) {
  # Initialize meta parameters - those that are common to all chains.
  args <- list(...)
  mcmc.input <- list()
  for (arg in names(args)) mcmc.input[[arg]] <- args[[arg]]
  
  meta <- do.meta.ini(mcmc.input, verbose=FALSE)
  return(structure(meta, class='bayesMig.mcmc.meta'))
}


do.meta.ini <- function(meta, burnin=200, verbose=FALSE) {
    start.year <- meta$start.year 
    present.year <- meta$present.year
    wpp.year <- meta$wpp.year
    if(!is.null(meta$my.mig.file))
      meta$my.mig.file <- normalizePath(meta$my.mig.file)
    my.mig.file <- meta$my.mig.file
    annual <- meta$annual.simulation
    #If the user input their own migration file:
     if(is.null(my.mig.file)){
       if(! wpp.year %in% c(2017, 2019, 2022)){
         stop("Only 2017, 2019 and 2022 revisions of WPP are currently supported by bayesMig.")
       }
       if(annual && wpp.year < 2022) 
         warning("If annual is TRUE and wpp.year is not 2022, 5-year data will be interpolated. Otherwise supply annual my.mig.file.")
     }
    migdata <- get.wpp.mig.data (start.year = start.year, present.year = present.year, 
                             wpp.year = wpp.year, my.mig.file = my.mig.file, 
                             annual = annual, exclude.from.world = meta$exclude.from.world,
                             verbose = verbose)

  #Establish some parameter constraints
    nC <- ncol(migdata$mig.matrix)
  muConstraints <- rep(NA, nC);
  phiConstraints <- rep(NA, nC);
  sigma2Constraints <- rep(NA, nC);
  
  #Format for optional fixes if we're including small countries
  # #Fix mu_c=0 for the following countries
  # muConstraints[fullCountryNameVec %in% c("Montserrat")] <- 0
  # #Fix mu_c=0.1 for the following countries
  # muConstraints[fullCountryNameVec %in% c("Holy See")] <- 0.1
  # #Fix phi_c=0 for the following countries
  # phiConstraints[fullCountryNameVec %in%
  #                  c("Montserrat","Andorra","Saint Pierre and Miquelon","Niue","Tokelau","Marshall Islands")] <- 0
  # #Fix sigma_c=0 for the following countries
  # sigma2Constraints[fullCountryNameVec %in% c("Montserrat")] <- 0
  
  #Compile constraints into logical vectors (that just say whether a constraint exists), and
  #  numeric vectors (that say what the parameter value must be)
  mu.constraints.logical <- !is.na(muConstraints)
  phi.constraints.logical <- !is.na(phiConstraints)
  sigma2.constraints.logical <- !is.na(sigma2Constraints)
  constraints.logical <- list(mu=mu.constraints.logical, phi=phi.constraints.logical, sigma2=sigma2.constraints.logical)
  constraints.numeric <- list(mu=muConstraints, phi=phiConstraints, sigma2=sigma2Constraints)

  # rename a few items
  meta$mu.global.lower <- meta$mu.range[1]
  meta$mu.global.upper <- meta$mu.range[2]
  meta$sigma.mu.lower <- meta$sigma.mu.range[1]
  meta$sigma.mu.upper <- meta$sigma.mu.range[2]
  meta$a.upper <- meta$a.up
  
  meta$mu.range <- NULL
  meta$sigma.mu.range <- NULL
  meta$a.up <- NULL

  return(c(meta, list(
    country.indices.est = 1:migdata$nr.countries.estimation,
    nr.countries.est = migdata$nr.countries.estimation,
    regions= c(migdata$regions['country_code'], migdata$regions['country_name']),
    mig.rates = t(migdata$mig.matrix),
    mig.rates.all = t(migdata$mig.matrix.all),
    user.data = !is.null(my.mig.file),
    bigT=nrow(migdata$mig.matrix),
    nr.countries = nC,
    #fullCountryCodeVec = migdata$regions$country_code,
    #fullCountryNameVec = migdata$regions$country_name,
    constraints.logical = constraints.logical,
    constraints.numeric = constraints.numeric
  )))
  
}

# ini MCMC for UN estimates

mcmc.ini <- function(chain.id, mcmc.meta,iter=1000) {
  nC <- mcmc.meta$nr.countries

  # Starting points
  mu_c <- rep(mcmc.meta$mu.ini[chain.id], nC)
  phi_c <- runif(nC, 0, 1)
  mu_global <- mcmc.meta$mu.ini[chain.id]
  a <- mcmc.meta$a.ini[chain.id]
  #b <- a - 1
  b <- (a-1)/10*mcmc.meta$prior.scaler/2

  sigma2_mu=var(as.numeric(rowMeans(mcmc.meta$mig.rates, na.rm = TRUE)))

  sigma2_c=as.numeric(apply(mcmc.meta$mig.rates,1,var, na.rm = TRUE))
  #Some countries have too many zeroes in their data
  sigma2_c[sigma2_c<mcmc.meta$sigma.c.min^2]=mcmc.meta$sigma.c.min^2
  
  #Force constraints for constrained parameters
  for(c in 1:nC){
    if(mcmc.meta$constraints.logical$mu[c]){
      mu_c[c]=mcmc.meta$constraints.numeric$mu[c];
    }
    if(mcmc.meta$constraints.logical$phi[c]){
      phi_c[c]=mcmc.meta$constraints.numeric$phi[c];
    }
    if(mcmc.meta$constraints.logical$sigma2[c]){
      sigma2_c[c]=mcmc.meta$constraints.numeric$sigma2[c];
    }
  }

  mcmc <- structure(list(
    meta=mcmc.meta,
    iter=iter,
    mu_c=mu_c,
    phi_c=phi_c,
    sigma2_c=sigma2_c,
    mu_global=mu_global,
    sigma2_mu=sigma2_mu,
    a=a,
    b=b,
    iter=iter,
    length=1,
    id=chain.id,
    output.dir=paste0('mc', chain.id),
    traces = 0, traces.burnin = 0
    ),
  class='bayesMig.mcmc')
  
  return(mcmc)
}

