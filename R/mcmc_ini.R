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
    meta$my.mig.file <- normalizePath(meta$my.mig.file)
    my.mig.file <- meta$my.mig.file
    annual <- meta$annual.simulation
    #If the user input their own migration file:
     if(!is.null(my.mig.file)){
    #   d <- read.delim(file=my.mig.file, comment.char='#', check.names=FALSE)
    # 
    #   #Extract country names and codes
    #   if(! "code" %in% colnames(d) && ! "country_code" %in% colnames(d))
    #       stop("Columns country_code or code must be present in the data file.")
    #   if("code" %in% colnames(d)) colnames(d)[colnames(d) == "code"] <- "country_code" # rename "code" to "country_code"
    #   fullCountryCodeVec=d$country_code
    #   if("country" %in% colnames(d)) colnames(d)[colnames(d) == "country"] <- "name" # rename country column to "name"
    #   if(! "name" %in% colnames(d)) d$name <- d$country_code
    #   fullCountryNameVec=d$name
    #   nC <- length(fullCountryCodeVec)
    # 
    #   #Extract migration rates
    #   mig.rates=as.matrix(d[,setdiff(colnames(d), c("country_code", "name"))])
    #   rownames(mig.rates)=fullCountryCodeVec
    # 
    #   if(!is.null(meta$exclude.from.world)) {
    #     #Exclude locations that should not influence the world parameters
    #     bigCountryIndices <- !fullCountryCodeVec %in% meta$exclude.from.world
    #   } else bigCountryIndices = rep(TRUE, nC)
    migdata <- get.wpp.mig.data (start.year = start.year, present.year = present.year, 
                             wpp.year = wpp.year, my.mig.file = my.mig.file, 
                             annual = annual, exclude.from.world = meta$exclude.from.world,
                             verbose = verbose)

  }else{
    #If we get here, then the user didn't input their own migration file.
    
    if(! wpp.year %in% c(2017, 2019, 2022)){
      #Only wpp2017 is supported at the moment.
      stop("Only 2017, 2019 and 2022 revisions of WPP are currently supported by bayesMig.")
    }
    if(annual && wpp.year < 2022) stop("If annual is TRUE and wpp.year is not 2022, my.mig.file must be provided. No default data available.")

    #If we get here, then wpp.year is 2017, 2019 and it is a 5-year simulation or wpp.year is 2022
    ###########
    
    #List of all possible countries
    UNlocations <- bayesTFR:::load.bdem.dataset('UNlocations', wpp.year=wpp.year)
    fullCountryCodeVec <- UNlocations$country_code[UNlocations$location_type==4]
    fullCountryNameVec <- UNlocations$name[UNlocations$location_type==4]
    
    #Pop and migration data
    pop <- bayesTFR:::load.from.wpp('pop', wpp.year=wpp.year, annual = annual)
    migration <- bayesTFR:::load.from.wpp('migration',wpp.year=wpp.year, annual = annual)

    #Figure out the countries of overlap
    fullDataIndices=(fullCountryCodeVec %in% migration$country_code & fullCountryCodeVec %in% pop$country_code)
    fullCountryCodeVec=fullCountryCodeVec[fullDataIndices]
    fullCountryNameVec=as.character(fullCountryNameVec[fullDataIndices])
    
    nC <- length(fullCountryCodeVec)
    
    if(!is.null(meta$exclude.from.world)) {
      #Exclude locations that should not influence the world parameters
      bigCountryIndices <- !fullCountryCodeVec %in% meta$exclude.from.world
    } else bigCountryIndices = rep(TRUE, nC)
    
    #Construct a matrix of initial populations
    initialPopMat <- merge(data.frame(country_code=fullCountryCodeVec), pop, sort=FALSE)
    numcols <- as.numeric(setdiff(colnames(initialPopMat), c("country_code", "name")))
    initialPopMat <- initialPopMat[,-which(colnames(pop) %in% c("country_code", "name", as.character(numcols[numcols <=1950])))] # need end-period pop
    rownames(initialPopMat) <- fullCountryCodeVec;
    
    #Convert from thousands to raw counts
    initialPopMat=initialPopMat*1000
    
    #Construct a matrix of total migration counts
    migCountMat <- merge(data.frame(country_code=fullCountryCodeVec), migration, sort=FALSE)[,-c(1,2)]
    if(!annual)
        migCountMat <- migCountMat[, substr(colnames(migCountMat), 6, 9) %in% colnames(initialPopMat)]
    else migCountMat <- migCountMat[, (as.integer(colnames(migCountMat)) + 1) %in% colnames(initialPopMat)]
    rownames(migCountMat)=fullCountryCodeVec;
    
    #Convert from thousands to raw counts
    migCountMat=migCountMat*1000
    
    #Convert migration counts and initial populations to a matrix of migration "rates"
    # as count/(end pop - mig)
    mig.rates<-as.matrix(migCountMat/(initialPopMat - migCountMat))
  }
    # # restrict dataset to columns between start and present year
    # if(start.year > present.year)
    #   stop("Arguments start.year must be smaller than present.year.")
    # 
    # cols.starty <- as.integer(substr(colnames(mig.rates), 1,4))
    # cols.endy <- if(annual) cols.starty+0.5 else as.integer(substr(colnames(mig.rates), 6,9))
    # start.index <- which((cols.starty <= start.year) & (cols.endy > start.year))
    # if(length(start.index) <= 0) {
    #   if(cols.starty[1] > start.year) start.index <- 1
    #   else stop("No data for time periods >=  ", start.year, " available. Check the argument start.year.")
    # }
    # start.col <- colnames(mig.rates)[start.index[1]]
    # present.index <- which((cols.endy >= present.year) & (cols.starty <= present.year))
    # if(length(present.index) <= 0) {
    #   if(cols.endy[length(cols.endy)] < present.year) present.index <- length(cols.endy)
    #   else stop("No data for time periods >=  ", start.year, " and <= ", present.year, " available. Check the arguments start.year and present.year.")
    # }
    # present.col <- colnames(mig.rates)[present.index[1]]
    # mig.rates <- mig.rates[, which.max(colnames(mig.rates)==start.col):which.max(colnames(mig.rates)==present.col)]
    # if(!annual) {
    #   # adjust the years to the right values
    #   meta$start.year <- cols.starty[start.index[1]]
    #   meta$present.year <- cols.endy[present.index[1]]
    #   # change the column names to the middle of the periods
    #   colnames(mig.rates) <- as.integer(substr(colnames(mig.rates), 1, 4)) + 3
    # }
    #   
  #Establish some parameter constraints
  muConstraints <- rep(NA, migdata$nr.countries.estimation);
  phiConstraints <- rep(NA, migdata$nr.countries.estimation);
  sigma2Constraints <- rep(NA, migdata$nr.countries.estimation);
  
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
  #stop("")
  return(c(meta, list(
    country.indices.est = 1:migdata$nr.countries.estimation,
    nr.countries.est = migdata$nr.countries.estimation,
    regions= c(migdata$regions['country_code'], migdata$regions['country_name']),
    mig.rates = t(migdata$mig.matrix),
    mig.rates.all = t(migdata$mig.matrix.all),
    user.data = !is.null(my.mig.file),
    bigT=nrow(migdata$mig.matrix),
    nr.countries=ncol(migdata$mig.matrix),
    fullCountryCodeVec = migdata$regions$country_code,
    fullCountryNameVec = migdata$regions$country_name,
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

  sigma2_mu=var(as.numeric(rowMeans(mcmc.meta$mig.rates)))

  sigma2_c=as.numeric(apply(mcmc.meta$mig.rates,1,var))
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

