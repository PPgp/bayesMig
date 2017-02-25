###############
#MIGRATION
###############



mcmc.meta.ini <- function(...) {
  # Initialize meta parameters - those that are common to all chains.
  args <- list(...)
  mcmc.input <- list()
  for (arg in names(args)) mcmc.input[[arg]] <- args[[arg]]
  
  meta <- do.meta.ini(mcmc.input, verbose=verbose)
  return(structure(c(mcmc.input, meta), class='bayesMig.mcmc.meta'))
}


do.meta.ini <- function(meta, burnin=200, verbose=FALSE,
                        wpp.year=2015, my.mig.file = NULL) {
  #If the user input their own migration file:
  if(!is.null(my.mig.file)){
    d=read.table(my.mig.file,header=TRUE)
    
    #Extract country names and codes
    fullCountryCodeVec=d$country_code
    fullCountryNameVec=d$country
    bigC=length(fullCountryCodeVec)
    
    #Extract migration rates
    mig.rates=d[,3:ncol(d)]
    rownames(mig.rates)=fullCountryCodeVec
    
    #If a user input their own rates, assume they want to use all of them.
    countryIndices=1:bigC
    big.country.indices=1:bigC
  }else{
    #If we get here, then the user didn't input their own migration file.
    
    if(! wpp.year == 2015){
      #Only wpp2015 is supported at the moment.
      error(stop("Only 2015 revision of WPP is supported by bayesMig."))
    }
    
    #If we get here, then wpp.year is 2015.
    ###########
    #I suspect this solution won't play nicely if the user has already loaded a different version of wpp.
    ###########
    if(wpp.year==2015){
      library(wpp2015)
    }
    
    #List of all possible countries
    data(UNlocations)
    fullCountryCodeVec <- UNlocations$country_code[UNlocations$location_type==4]
    fullCountryNameVec <- UNlocations$name[UNlocations$location_type==4]
    
    #Pop and migration data
    data(pop)
    data(migration);
    data(countryCodeVec_bigCountries)# <- scan("./Data/countryCodeVec_bigCountries.txt")#This is the 201 "big" countries
    
    #Figure out the countries of overlap
    fullDataIndices=(fullCountryCodeVec %in% migration$country_code & fullCountryCodeVec %in% pop$country_code)
    fullCountryCodeVec=fullCountryCodeVec[fullDataIndices]
    fullCountryNameVec=as.character(fullCountryNameVec[fullDataIndices])
    
    countryIndices=seq(1,length(fullCountryCodeVec))
    bigCountryIndices <- fullCountryCodeVec %in% countryCodeVec_bigCountries[,1]
    
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
    migCountMat <- merge(data.frame(country_code=fullCountryCodeVec), migration[,1:which(colnames(migration)=="2010-2015")], sort=FALSE)[,-c(1,2)]
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
    mig.rates<-migCountMat/(initialPopMat - migCountMat)
  }

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
  # #Fix sigma_c=0 for the collowing countries
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
    present.year = 2015,
    bigT=ncol(mig.rates),
    nr.countries=nrow(mig.rates),
    fullCountryCodeVec = fullCountryCodeVec,
    fullCountryNameVec = fullCountryNameVec,
    constraints.logical = constraints.logical,
    constraints.numeric = constraints.numeric,
    sigma.c.min = 0.0001,
    mu.global.lower = -0.5,
    mu.global.upper = 0.5,
    sigma.mu.lower=0,
    sigma.mu.upper=0.5,
    a.upper=10
  ))
  
}

# ini MCMC for UN estimates

mcmc.ini <- function(chain.id, mcmc.meta,iter=1000) {
  bigC=nrow(mcmc.meta$mig.rates)

  #Reasonable starting points.
  #(TO DO: Allow these to take different values)
  mu_c=rep(0,bigC)
  phi_c=rep(0.5,bigC)
  mu_global=0
  a=1.5
  b=0.5
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

###############
#Test code goes here.
###############

#mcmc.meta=mcmc.meta.ini()
#mcmc=mcmc.ini(mcmc.meta)
#names(mcmc$meta)
#names(mcmc)
#mcenv=as.environment(mcmc)


#Write example migration file based on wpp2015 data.
# d=data.frame(country=fullCountryNameVec,
#              country_code=fullCountryCodeVec,
#              mig.rates)
# colnames(d)=c("country","country_code",
#               "1950-1955","1955-1960","1960-1965","1965-1970",
#               "1970-1975","1975-1980","1980-1985","1985-1990",
#               "1990-1995","1995-2000","2000-2005","2005-2010",
#               "2010-2015")
# write.table(d,"my_mig_example.txt",sep="\t",row.names=FALSE)
