#' @title Generate posterior trajectories of net migration rates
#'
#' @description Using the posterior parameter samples simulated by \code{\link{run.mig.mcmc}},
#' generate posterior trajectories for the net migration rates for all countries of
#' the world. This code \emph{does not} adjust trajectories to ensure that net
#' migration counts sum to zero. That adjustment is handled in the \code{bayesPop} package.
#' 
#' @param mcmc.set Object of class \code{bayesMig.mcmc.set} corresponding to sampled
#' parameter values for net migration model. If it is \code{NULL}, the object
#' is loaded from the directory specified in \code{sim.dir}
#' @param end.year End year of the prediction
#' @param sim.dir Directory with MCMC simulation results. It should be the same as
#' the \code{output.dir} argument in \code{\link{run.mig.mcmc}}
#' @param replace.output Logical value. If \code{TRUE}, existing predictions in
#' \code{output.dir} will be replaced by results of this run.
#' @param start.year Start year of the prediction. By default the prediction is 
#' started at the next time period after \code{present.year} set in the estimation
#' step. If \code{start.year} is smaller than the default, projections for countries
#'and time periods that have data available after \code{start.year} are set to those data.
#' @param nr.traj Number of trajectories to be generated. 
#' If \code{NULL}, the argument \code{thin} is taken to determine the number of 
#' trajectories. If both are \code{NULL}, the number of trajectories
#' corresponds to the size of the parameter sample.
#' @param thin Thinning interval used for determining the number of trajectories. 
#' Only relevant if \code{nr.traj} is \code{NULL}.
#' @param burnin Number of iterations to be discarded from the beginning of the parameter traces.
#' @param save.as.ascii Either a number determining how many trajectories should be
#' converted into an ASCII file, or 'all' in which case all trajectories are converted.
#' It should be set to 0 if no conversion is desired.
#' @param output.dir Directory into which the resulting prediction object and the 
#' trajectories are stored. If it is \code{NULL}, it is set to either \code{sim.dir},
#' or to \code{output.dir} of \code{mcmc.set$meta} if \code{mcmc.set} is given.
#' @param seed Seed of the random number generator. If \code{NULL} no seed is set. 
#' Can be used to generate reproducible projections.
#' @param verbose Logical value. Switches log messages on and off.
#' @param ... Other arguments passed to \code{\link{mig.predict}}
#' @return A list with 9 components. Key result component is an array of quantiles with dimensions
#' (number of countries) x (number of computed quantiles) x (number of projected time points).
#' First time point in the sequence is not a projection, but an observed time period -- by default, 2010-2015.
#' 
#' Other key result components include \code{traj.mean.sd}, a summary of means and standard deviations for each country
#' at each time point.

mig.predict <- function(mcmc.set=NULL, end.year=2100,
						sim.dir=file.path(getwd(), 'bayesMig.output'),
						replace.output=FALSE,
						start.year=NULL, nr.traj = NULL, thin = NULL, burnin=2000, 
						save.as.ascii=1000, output.dir = NULL,
						seed=NULL, verbose=TRUE, ...) {
	if(!is.null(mcmc.set)) {
		if (class(mcmc.set) != 'bayesMig.mcmc.set') {
			stop('Wrong type of mcmc.set. Must be of type bayesMig.mcmc.set.')
			}
	} else {		
		mcmc.set <- get.mig.mcmc(sim.dir, verbose=verbose)
	}

	if(!is.null(seed)) set.seed(seed)
	
	invisible(make.mig.prediction(mcmc.set, end.year=end.year, replace.output=replace.output,  
					start.year=start.year, nr.traj=nr.traj, burnin=burnin, thin=thin,
					output.dir=output.dir, verbose=verbose, ...))			
}

make.mig.prediction <- function(mcmc.set, start.year=NULL, end.year=2100, replace.output=FALSE,
								nr.traj = NULL, burnin=0, thin = NULL, 
								countries = NULL,
							    save.as.ascii=1000, output.dir = NULL, write.summary.files=TRUE, 
							    is.mcmc.set.thinned=FALSE, force.creating.thinned.mcmc=FALSE,
							    write.trajectories=TRUE, 
							    verbose=verbose){
	# if 'countries' is given, it is an index
	meta <- mcmc.set$meta
	present.year <- if(is.null(start.year)) meta$present.year else start.year - 5
	nr_project <- length(seq(present.year+5, end.year, by=5))
	cat('\nPrediction from', present.year+5, 'until', end.year, '(i.e.', nr_project, 'projections)\n\n')

	burn <- if(is.mcmc.set.thinned) 0 else burnin
	total.iter <- get.total.iterations(mcmc.set$mcmc.list, burn)
	stored.iter <- get.stored.mcmc.length(mcmc.set$mcmc.list, burn)
	mcthin <- max(sapply(mcmc.set$mcmc.list, function(x) x$thin))
	if(!is.null(nr.traj) && !is.null(thin)) {
		warning('Both nr.traj and thin are given. Argument thin will be ignored.')
		thin <- NULL
	}
	if(is.null(nr.traj)) nr.traj <- min(stored.iter, 2000)
	else {
		if (nr.traj > stored.iter) 
			warning('nr.traj is larger than the available MCMC sample. Only ', stored.iter, ' trajectories will be generated.')
		nr.traj <- min(nr.traj, stored.iter)	
	}
	if(is.null(thin)) thin <- floor(stored.iter/nr.traj * mcthin)
	if(stored.iter <= 0 || thin == 0)
		stop('The number of simulations is 0. Burnin might be larger than the number of simulated values, or # trajectories is too big.')
	
	#setup output directory
	if (is.null(output.dir)) output.dir <- meta$output.dir
	outdir <- file.path(output.dir, 'predictions')
	
	if(is.null(countries)) {
		if(!replace.output && has.mig.prediction(sim.dir=output.dir))
			stop('Prediction in ', outdir,
			' already exists.\nSet replace.output=TRUE if you want to overwrite existing projections.')
		unlink(outdir, recursive=TRUE)
		write.to.disk <- TRUE
		if(!file.exists(outdir)) 
			dir.create(outdir, recursive=TRUE)
	} else write.to.disk <- FALSE
	
	if(is.mcmc.set.thinned) {
		thinned.mcmc <- mcmc.set
		has.thinned.mcmc <- TRUE
	} else {
		thinned.mcmc <- get.thinned.mig.mcmc(mcmc.set, thin=thin, burnin=burnin)
		has.thinned.mcmc <- !is.null(thinned.mcmc) #&& thinned.mcmc$meta$parent.iter == total.iter#JA: I don't understand this constraint
	}
	#unblock.gtk('bDem.migpred')#JA: ???
	if(has.thinned.mcmc && !force.creating.thinned.mcmc){
	  load.mcmc.set=thinned.mcmc
	}else{
	  load.mcmc.set=create.thinned.mig.mcmc(mcmc.set, thin=thin, burnin=burnin, output.dir=output.dir, verbose=verbose)
	}
	
	nr_simu <- load.mcmc.set$mcmc.list[[1]]$finished.iter

	prediction.countries <- if(is.null(countries)) 1:meta$nr.countries else countries
	nr_countries <- meta$nr.countries
	nr_countries_real <- length(prediction.countries)
	#JA: Assume we're only producing genuine predictions (i.e. starting at 2015-2020.)
	#present.year.index <- get.estimation.year.index(meta, present.year)

	#keep these defaults for checking the out-of-sample projections
  quantiles.to.keep <- c(0,0.025,0.05,0.1,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.9,0.95,0.975,1)
	PIs_cqp <- array(NA, c(nr_countries, length(quantiles.to.keep), nr_project+1))
	dimnames(PIs_cqp)[[2]] <- quantiles.to.keep
	proj.middleyears <- bayesTFR:::get.prediction.years(meta, nr_project+1)
	dimnames(PIs_cqp)[[3]] <- proj.middleyears
	mean_sd <- array(NA, c(nr_countries, 2, nr_project+1))
	hasNAs <- rep(FALSE, nr_simu)

  cs.par.values.list = list()
	# # country loop for preparing data for projections
	for (country in prediction.countries){
		country.obj <- get.country.object(country, meta, index=TRUE)
		cs.par.values <- get.mig.parameter.traces.cs(load.mcmc.set$mcmc.list, country.obj, 
		                                                                 mig.parameter.names.cs(), burnin=0)
		cs.par.values.list[[country]] <- cs.par.values
		
	} # end country prep loop
	
	# array for results - includes also historical data for periods with missing data
	all.mig_ps <- array(NA, dim=c(nr_countries_real, nr_project+1, nr_simu))
  #JA: This code assumes the migration rates in meta$mig.rates end at 2010-2015 and will be projected until 2100.

	# fill the result array with observed data 
  for(country in prediction.countries){
    country.obj <- get.country.object(country, meta, index=TRUE)
    all.mig_ps[country, 1,] = meta$mig.rates[which(rownames(meta$mig.rates) == country.obj$code), as.character(present.year - 2)]
  }

	mu.c <- phi.c <- sigma.c <- rep(NA, nr_countries)

	traj.counter <- 0
	country.loop.max <- 20
	if (verbose) {
		verbose.iter <- max(1, nr_simu/100)
		if(interactive()) cat('\n')
	}
	
	#########################################
	for (s in 1:nr_simu){ # Iterate over trajectories
	#########################################
		if (verbose) {
			if(interactive()) cat('\rProjecting migration trajectories ... ', round(s/nr_simu * 100), ' %')
			else {
				if (s %% verbose.iter == 0) 
					cat('Migration projection trajectory ', s, '\n')
				}
		}
	  #Pull parameter values for this trajectory
	  for (country in prediction.countries){		
	    mu.c[country] <- cs.par.values.list[[country]][s,1]
	    phi.c[country] <- cs.par.values.list[[country]][s,2]
	    sigma.c[country] <- sqrt(cs.par.values.list[[country]][s,3])
	  }
	  
	  #########################################
	  for (year in 2:(nr_project+1)) { # Iterate over time
	    #########################################
	    #ALLmig.prev <- all.mig_ps[,year-1,s]
	    #########################################
	    for (icountry in 1:nr_countries_real){ # Iterate over countries
	      #########################################
	      all.mig_ps[icountry,year,s]=mu.c[icountry] + phi.c[icountry]*(all.mig_ps[icountry,year-1,s] - mu.c[icountry]) + rnorm(n=1,mean=0,sd=sigma.c[icountry])
	    } # end countries loop
	  } # end time loop
	} # end simu loop
	if(verbose && interactive()) cat('\n')

	##############
	# Compute quantiles
	for (icountry in 1:nr_countries_real){
		country <- prediction.countries[icountry]
		country.obj <- get.country.object(country, meta, index=TRUE)
		
		# extract the future trajectories (including the present period)
		mig_ps_future <- all.mig_ps[icountry,(dim(all.mig_ps)[2]-nr_project):dim(all.mig_ps)[2],]

		if(write.trajectories) {
			trajectories <- mig_ps_future # save only trajectories simulated for the future time
  			save(trajectories, file = file.path(outdir, paste('traj_country', country.obj$code, '.rda', sep='')))
  		}
  		# compute quantiles
 		PIs_cqp[country,,] <- apply(mig_ps_future, 1, quantile, quantiles.to.keep, na.rm = TRUE)
 		mean_sd[country,1,] <- apply(mig_ps_future, 1, mean, na.rm = TRUE)
 		mean_sd[country,2,] <- apply(mig_ps_future, 1, sd, na.rm = TRUE)
 	}
	mcmc.set <- remove.mig.traces(mcmc.set)
	bayesMig.prediction <- structure(list(
				quantiles = PIs_cqp,
				traj.mean.sd = mean_sd,
				nr.traj=nr_simu,
				output.directory=outdir,
				mcmc.set=load.mcmc.set,
				nr.projections=nr_project,
				burnin=burnin, thin=thin,
				#mu=mu, rho=rho,  sigma_t = sigmas_all, sigmaAR1 = sigmaAR1,
				end.year=end.year),
				class='bayesMig.prediction')
			
	if(write.to.disk) {
		store.bayesMig.prediction(bayesMig.prediction, outdir)
	    bayesTFR:::do.convert.trajectories(pred=bayesMig.prediction, n=save.as.ascii, output.dir=outdir, verbose=verbose)
		#if(write.summary.files)
		  #JA: Used to be tfr.write.projection.summary.and.parameters
		  #    Parameter summary didn't translate nicely, so now the default only writes projection summary
		mig.write.projection.summary(pred=bayesMig.prediction, output.dir=outdir)
		cat('\nPrediction stored into', outdir, '\n')
	}
	invisible(bayesMig.prediction)
}


remove.mig.traces <- function(mcmc.set) {
	for (i in 1:length(mcmc.set$mcmc.list))
		mcmc.set$mcmc.list[[i]]$traces <- 0
	invisible(mcmc.set)
}

get.traj.ascii.header.bayesMig.mcmc.meta <- function(meta, ...) 
	return (list(country_code='LocID', period='Period', year='Year', trajectory='Trajectory', tfr='Mig'))
		

mig.write.projection.summary <- function(pred, output.dir) {
	# one summary file
	#do.write.projection.summary(pred, output.dir)
    bayesTFR:::do.write.projection.summary(pred, output.dir)
}

get.estimation.years <- function(meta)
  return(as.numeric(colnames(meta$mig.rates))+2.5)


get.projection.summary.header.bayesMig.prediction <- function(pred, ...) 
  return (list(revision='RevID', variant='VarID', country='LocID', year='TimeID', indicator='IndicatorID', sex='SexID', tfr='Value'))

get.friendly.variant.names.bayesMig.prediction <- function(pred, ...)
  return(c('median', 'lower 80', 'upper 80', 'lower 95', 'upper 95','constant'))

get.UN.variant.names.bayesMig.prediction <- function(pred, ...) 
    return(c('BHM median', 'BHM80 lower',  'BHM80 upper', 'BHM95 lower',  'BHM95 upper', 'Zero migration'))


get.mig.periods <- function(meta) {
  mid.years <- get.estimation.years(meta)
  return (paste(mid.years-2.5, mid.years+2.5, sep='-'))
}

get.data.imputed.bayesMig.prediction <- function(pred, ...)
    return(get.data.matrix(pred$mcmc.set$meta))

get.data.for.country.imputed.bayesMig.prediction <- function(pred, country.index, ...)
    return(get.data.matrix(pred$mcmc.set$meta)[, country.index])


