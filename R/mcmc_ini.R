###############
#MIGRATION
###############



mcmc.meta.ini <- function(...) {
  # Initialize meta parameters - those that are common to all chains.
  args <- list(...)
  mcmc.input <- list()
  for (arg in names(args)) mcmc.input[[arg]] <- args[[arg]]
  
  meta <- do.meta.ini(mcmc.input, verbose=FALSE)
  return(structure(c(mcmc.input, meta), class='bayesMig.mcmc.meta'))
}


do.meta.ini <- function(meta, burnin=200, verbose=FALSE) {
    start.year <- meta$start.year # TODO: take this into account
    present.year <- meta$present.year
    wpp.year <- meta$wpp.year
    my.mig.file <- meta$my.mig.file
  #If the user input their own migration file:
  if(!is.null(my.mig.file)){
      d <- read.delim(file=my.mig.file, comment.char='#', check.names=FALSE)
    
    #Extract country names and codes
    fullCountryCodeVec=d$country_code
    if("country" %in% colnames(d)) colnames(d)[colnames(d) == "country"] <- "name" # rename country column to "name"
    if(! "name" %in% colnames(d)) d$name <- d$country_code
    fullCountryNameVec=d$name
    bigC=length(fullCountryCodeVec)
    
    #Extract migration rates
    mig.rates=as.matrix(d[,setdiff(colnames(d), c("country_code", "name"))])
    rownames(mig.rates)=fullCountryCodeVec
    
    #If a user input their own rates, assume they want to use all of them.
    countryIndices=1:bigC
    bigCountryIndices = rep(TRUE, bigC)
  }else{
    #If we get here, then the user didn't input their own migration file.
    
    if(! wpp.year %in% c(2017, 2019)){
      #Only wpp2017 is supported at the moment.
      stop("Only 2017 and 2019 revisions of WPP are currently supported by bayesMig.")
    }
    
    #If we get here, then wpp.year is 2017 or 2019.
    ###########
    
    #List of all possible countries
    UNlocations=load.bdem.dataset('UNlocations', wpp.year=wpp.year)
    fullCountryCodeVec <- UNlocations$country_code[UNlocations$location_type==4]
    fullCountryNameVec <- UNlocations$name[UNlocations$location_type==4]
    
    #Pop and migration data
    pop=load.bdem.dataset('pop', wpp.year=wpp.year)
    migration=load.bdem.dataset('migration',wpp.year=wpp.year)
    countryCodeVec_bigCountries = bayesMig::countryCodeVec_bigCountries
    #data("countryCodeVec_bigCountries")# <- scan("./Data/countryCodeVec_bigCountries.txt")#This is the 201 "big" countries
    
    #Figure out the countries of overlap
    fullDataIndices=(fullCountryCodeVec %in% migration$country_code & fullCountryCodeVec %in% pop$country_code)
    fullCountryCodeVec=fullCountryCodeVec[fullDataIndices]
    fullCountryNameVec=as.character(fullCountryNameVec[fullDataIndices])
    
    countryIndices=seq(1,length(fullCountryCodeVec))
    bigCountryIndices <- fullCountryCodeVec %in% countryCodeVec_bigCountries
    
    #Construct a matrix of initial populations
    initialPopMat <- merge(data.frame(country_code=fullCountryCodeVec), pop, sort=FALSE)
    initialPopMat <- initialPopMat[,-which(colnames(pop) %in% c("country_code", "name", "1950"))] # need end-period pop
    
    #initialPopMat=matrix(0,nrow=length(fullCountryCodeVec),ncol=14)
    rownames(initialPopMat)=fullCountryCodeVec;
    #colnames(initialPopMat)=seq(1950,2015,5);
    #for(i in 1:length(fullCountryCodeVec)){
    #  initialPopMat[i,]=colSums(popM[popM$country_code==fullCountryCodeVec[i],4:17])+colSums(popF[popF$country_code==fullCountryCodeVec[i],4:17])
    #}
    
    #Convert from thousands to raw counts
    initialPopMat=initialPopMat*1000
    
    #Construct a matrix of total migration counts
    migCountMat <- merge(data.frame(country_code=fullCountryCodeVec), 
                         migration[,1:which(colnames(migration)==paste0(present.year - 5, "-", present.year))], sort=FALSE)[,-c(1,2)]
    #migCountMat=matrix(0,nrow=length(fullCountryCodeVec),ncol=13)
    rownames(migCountMat)=fullCountryCodeVec;
    #colnames(migCountMat)=seq(1950,2010,5);
    
    #for(i in 1:length(fullCountryCodeVec)){
    #  migCountMat[i,]=as.numeric(migration[which(migration$country_code==fullCountryCodeVec[i]),3:15])
    #}
    
    #Convert from thousands to raw counts
    migCountMat=migCountMat*1000
    
    #Convert migration counts and initial populations to a matrix of migration "rates" (count/initial pop)
    # mig.rates<-migCountMat[,1:13]/initialPopMat[,1:13]  
    # do count/(end pop - mig)
    mig.rates<-as.matrix(migCountMat/(initialPopMat - migCountMat))
  }
    
  colnames(mig.rates) <- as.integer(substr(colnames(mig.rates), 1, 4)) + 3 # set column names to the middle of the periods
  #Establish some parameter constraints
  muConstraints <- rep(NA,length(fullCountryNameVec));
  phiConstraints <- rep(NA,length(fullCountryNameVec));
  sigma2Constraints <- rep(NA,length(fullCountryNameVec));
  
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

  return(list(
    countryIndices=countryIndices,
    big.country.indices=bigCountryIndices,
    nr.big.countries=sum(bigCountryIndices),
    regions=list(country_code=fullCountryCodeVec,country_name=fullCountryNameVec),
    mig.rates = mig.rates,
    start.year = start.year,
    present.year = present.year,
    user.data = !is.null(my.mig.file),
    bigT=ncol(mig.rates),
    nr.countries=nrow(mig.rates),
    fullCountryCodeVec = fullCountryCodeVec,
    fullCountryNameVec = fullCountryNameVec,
    constraints.logical = constraints.logical,
    constraints.numeric = constraints.numeric,
    sigma.c.min = meta$sigma.c.min,
    mu.global.lower = meta$mu.range[1],
    mu.global.upper = meta$mu.range[2],
    sigma.mu.lower=meta$sigma.mu.range[1],
    sigma.mu.upper=meta$sigma.mu.range[2],
    a.upper=meta$a.up,
    mu.ini = meta$mu.ini,
    a.ini = meta$a.ini, 
    a.half.width = meta$a.half.width,
    buffer.size = meta$buffer.size
  ))
  
}

# ini MCMC for UN estimates

mcmc.ini <- function(chain.id, mcmc.meta,iter=1000) {
  bigC=nrow(mcmc.meta$mig.rates)

  # Starting points
  mu_c <- rep(mcmc.meta$mu.ini[chain.id], bigC)
  phi_c <- runif(bigC, 0, 1)
  mu_global <- mcmc.meta$mu.ini[chain.id]
  a <- mcmc.meta$a.ini[chain.id]
  b <- a - 1

  sigma2_mu=var(as.numeric(rowMeans(mcmc.meta$mig.rates)))

  sigma2_c=as.numeric(apply(mcmc.meta$mig.rates,1,var))
  #Some countries have too many zeroes in their data
  sigma2_c[sigma2_c<mcmc.meta$sigma.c.min^2]=mcmc.meta$sigma.c.min^2
  
  #Force constraints for constrained parameters
  for(c in 1:bigC){
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
    output.dir=paste('mc', chain.id, sep='')
    ),
  class='bayesMig.mcmc')
  
  return(mcmc)
}


mig.parameter.names <- function() {
  # Return all country-independent parameter names.
  return(c("a","b","mu_global","sigma2_mu"))
}

mig.parameter.names.cs <- function(){
  # Return all country-specific parameter names.
  #That is, these parameters are vectors of length bigC.
  return(c("mu_c","phi_c","sigma2_c"))
}

