#' @title Generate posterior trajectories of net migration rates
#'
#' @description Using the posterior parameter samples simulated by \code{\link{run.mig.mcmc}},
#' generate posterior trajectories for the net migration rates for all countries of
#' the world, or all locations included in the estimation. This code \emph{does not} adjust trajectories to ensure that net
#' migration counts over all countries sum to zero. 
#' 
#' @param mcmc.set Object of class \code{bayesMig.mcmc.set} corresponding to sampled
#' parameter values for net migration model. If it is \code{NULL}, the object
#' is loaded from the directory specified in \code{sim.dir}
#' @param end.year End year of the prediction
#' @param sim.dir Directory with MCMC simulation results. It should be the same as
#' the \code{output.dir} argument in \code{\link{run.mig.mcmc}}
#' @param replace.output Logical value. If \code{TRUE}, existing predictions in
#' \code{output.dir} will be replaced by results of this run.
#' @param start.year Start year of the prediction, i.e. the first predicted year. By default the prediction is 
#' started at the next time period after \code{present.year} set in the estimation
#' step. If \code{start.year} is smaller than the default, the behavior is controlled by 
#' the \code{post.last.observed} argument: Either data post \code{start.year} is ignored (default)
#' or the projection is set to the available data (\code{post.last.observed = "a"}).
#' @param nr.traj Number of trajectories to be generated. 
#' If \code{NULL}, the argument \code{thin} is taken to determine the number of 
#' trajectories. If both are \code{NULL}, the number of trajectories
#' corresponds to the size of the parameter sample.
#' @param thin Thinning interval used for determining the number of trajectories. 
#' Only relevant if \code{nr.traj} is \code{NULL}.
#' @param burnin Number of iterations to be discarded from the beginning of the parameter traces.
#' @param use.cummulative.threshold If \code{TRUE} historical cummulative thresholds are applied
#'    to avoid sampling rates that are too extreme. The thresholds are
#'    derived over prior rates of all locations. As a time span for deriving the limits on projected rates, 
#'    at each projected time point, six prior time periods are used in a 5-year simulation, 
#'    corresponding to 30 years in an annual simulation.
#'    In a national simulation, prior rates of GCC countries (plus Western Sahara and Djibouti) are excluded 
#'    when deriving thresholds for non-GCC countries. If this option is used in a non-country simulation,
#'    e.g. in a sub-national settings, set the \code{ignore.gcc.in.threshold} argument to \code{TRUE}.
#' @param ignore.gcc.in.threshold If \code{use.cummulative.threshold} is \code{TRUE}, by default the GCC countries
#'    (plus Western Sahara and Djibouti) identified by numerical codes of the countries are excluded from computing 
#'    the historical cummulative thresholds for non-GCC countries. If this argument is \code{TRUE}, this distinction is not made. 
#'    It is important to set it to \code{TRUE} in a sub-national simulation to avoid any random overlaps 
#'    of UN codes and user-defined codes.
#' @param post.last.observed If a user-specific data file was used during estimation and the data 
#'     contained the \dQuote{last.observed} column, this argument determines how to treat the time periods 
#'     between the last observed point and the start year of the prediction, for locations where there is
#'     a gap between them, or if short-term predictions were included in the file. It is also relevant
#'     if \code{start.year} is set to a smaller value than \code{present.year} in the estimation.
#'     Possible values are:
#' \itemize{    
#' \item \dQuote{obsdata} or \dQuote{o} (default) uses any non-missing observed data 
#'     provided in the data file during estimation, up to the time point defined by the argument \code{start.year} 
#'     (excluding the start year itself). 
#'     
#' \item \dQuote{alldata} or \dQuote{a} would similarly use 
#'     the provided data but would use all data, even if it goes beyond the start year. This allows
#'     to use short-term deterministic projections for locations where it is available. 
#'     
#' \item \dQuote{impute} or \dQuote{i} would ignore all data beyond the last observed data point 
#'     and impute the missing time periods.
#' }
#' @param save.as.ascii Either a number determining how many trajectories should be
#' converted into an ASCII file, or 'all' in which case all trajectories are converted.
#' It should be set to 0 if no conversion is desired. If this argument 
#' is larger than zero, the resulting file can be used as input into population projection via \pkg{bayesPop},
#' see Details.
#' @param output.dir Directory into which the resulting prediction object and the 
#' trajectories are stored. If it is \code{NULL}, it is set to either \code{sim.dir},
#' or to \code{output.dir} of \code{mcmc.set$meta} if \code{mcmc.set} is given.
#' @param seed Seed of the random number generator. If \code{NULL} no seed is set. 
#' Can be used to generate reproducible projections.
#' @param verbose Logical value. Switches log messages on and off.
#' @param \dots Further arguments passed to the underlying functions.
#' 
#' @details The trajectories of net migration rates for each location are generated using the model of Azose & Raftery (2015).
#' Parameter samples  simulated via \code{\link{run.mig.mcmc}} are used from all chains, from which the given \code{burnin} 
#' was discarded. They are evenly thinned to match \code{nr.traj} or using the \code{thin} argument. 
#' Such thinned parameter traces, collapsed into one chain, if they do not already exist, are stored on disk 
#' into the sub-directory \file{thinned_mcmc_\emph{t}_\emph{b}} where \emph{t} is the value  of \code{thin} and 
#' \emph{b} the value of \code{burnin}.
#' 
#' The projection is run for all missing values before the present year, if any. 
#' Medians over the trajectories are used as imputed values and the trajectories are discarded. 
#' The process then continues by projecting the future values where all generated trajectories are kept.
#' 
#' A special case is when the argument \code{start.year} is given that is smaller than or equal to
#' the present year. In such a case, imputed missing values before present year are treated 
#' as ordinary predictions (trajectories are kept). If \code{post.last.observed} is \dQuote{a}, 
#' all historical data between start year and present year are used as projections.
#' 
#' The resulting prediction object is saved into \file{\{output.dir\}/predictions}. Trajectories 
#' for all locations are saved into the same directory in a binary format, one file per location. 
#' At the end of the projection, if \code{save.as.ascii} is larger than 0, the function converts 
#' the given number of trajectories into a CSV file, called \file{ascii_trajectories.csv} also located
#' in the \file{predictions} directory. The converted trajectories are selected by equal spacing. 
#' In addition to the converted trajectories, two summary files are created: one in a user-friendly format, the other using 
#' a UN-specific coding, as described in \code{\link{mig.write.projection.summary}}.
#' 
#' If it is desired to use these predictions as input to population projections in \pkg{bayesPop},
#' enter the full file path of the \file{ascii_trajectories.csv} file into the \code{inputs} argument 
#' of \code{bayesPop::pop.predict} as item \code{migtraj} and set the argument \code{mig.is.rate} appropriately.
#' 
#' @return Object of class \code{bayesMig.prediction} which is a list with components
#' containing details of the prediction. Key result component is an array of \code{quantiles} with dimensions
#' (number of locations) x (number of computed quantiles) x (number of projected time points).
#' First time point in the sequence is not a projection, but the last observed time period.
#' 
#' Other key result components include \code{traj.mean.sd}, a summary of means and standard deviations for each country
#' at each time point. See \code{\link[bayesTFR]{bayesTFR.prediction}} for more detail.
#' 
#' @references Azose, J. J., & Raftery, A. E. (2015). 
#' Bayesian probabilistic projection of international migration. Demography, 52(5), 1627-1650.
#' \doi{10.1007/s13524-015-0415-0}.
#' 
#' Azose, J.J., Ševčíková, H., Raftery, A.E. (2016): Probabilistic population projections with migration uncertainty. 
#' Proceedings of the National Academy of Sciences 113:6460–6465. \doi{10.1073/pnas.1606119113}.
#' 
#' @examples
#' # Toy simulation for US states
#' us.mig.file <- file.path(find.package("bayesMig"), "extdata", "USmigrates.txt")
#' sim.dir <- tempfile()
#' m <- run.mig.mcmc(nr.chains = 2, iter = 30, thin = 1, my.mig.file = us.mig.file, 
#'         output.dir = sim.dir, present.year = 2017, annual = TRUE)
#' 
#' # Prediction
#' pred <- mig.predict(sim.dir = sim.dir, burnin = 5, end.year = 2050)
#
#' # here unrealistic results since this is a toy simulation 
#' mig.trajectories.plot(pred, "Hawaii", pi = 80, ylim = c(-0.02, 0.02)) 
#' mig.trajectories.table(pred, "Hawaii")
#' summary(pred, "California")
#' 
#' # view locations included in the simulation
#' get.countries.table(pred)
#' 
#' unlink(sim.dir, recursive = TRUE)
#' # For projections on national level, see ?bayesMig.
#' @aliases bayesMig.prediction
#' @export

mig.predict <- function(mcmc.set=NULL, end.year=2100,
						sim.dir = NULL,
						replace.output=FALSE,
						start.year=NULL, nr.traj = NULL, thin = NULL, burnin=20000, 
						use.cummulative.threshold = FALSE, ignore.gcc.in.threshold = FALSE,
						fixed.thresholds = NULL,
						post.last.observed = c("obsdata", "alldata", "impute"),
						save.as.ascii=0, output.dir = NULL,
						seed=NULL, verbose=TRUE, ...) {
	if(!is.null(mcmc.set)) {
		if (! inherits(mcmc.set, 'bayesMig.mcmc.set')) {
			stop('Wrong type of mcmc.set. Must be of type bayesMig.mcmc.set.')
			}
	} else {
	    if(is.null(sim.dir))
	        stop("Either mcmc.set or sim.dir argument must be given to locate the estimation results.")
		mcmc.set <- get.mig.mcmc(sim.dir, verbose=verbose)
	}

    post.last.observed <- substr(match.arg(post.last.observed), 1, 1)
    
	if(!is.null(seed)) set.seed(seed)
	
	invisible(make.mig.prediction(mcmc.set, end.year=end.year, replace.output=replace.output,  
					start.year=start.year, nr.traj=nr.traj, burnin=burnin, thin=thin,
					use.cummulative.threshold = use.cummulative.threshold, ignore.gcc.in.threshold = ignore.gcc.in.threshold,
					fixed.thresholds = fixed.thresholds, post.last.observed = post.last.observed,
					save.as.ascii=save.as.ascii, output.dir=output.dir, verbose=verbose, ...))			
}

make.mig.prediction <- function(mcmc.set, start.year=NULL, end.year=2100, replace.output=FALSE,
								nr.traj = NULL, burnin=0, thin = NULL, 
								countries = NULL, use.cummulative.threshold = FALSE, ignore.gcc.in.threshold = FALSE,
								fixed.thresholds = NULL, post.last.observed = "o",
								save.as.ascii=0, output.dir = NULL, write.summary.files=TRUE, 
							    is.mcmc.set.thinned=FALSE, force.creating.thinned.mcmc=FALSE,
							    write.trajectories=TRUE, 
							    verbose=verbose){
	# if 'countries' is given, it is an index
	meta <- mcmc.set$meta
	year.step <- ifelse(meta$annual.simulation, 1, 5)
	present.year <- if(is.null(start.year)) meta$present.year else start.year - year.step
	nr_project <- length(seq(present.year+year.step, end.year, by = year.step))
	cat('\nPrediction from', present.year+year.step, 'until', end.year, '(i.e.', nr_project, 'projections)\n\n')

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

	present.year.index <- bayesTFR:::get.estimation.year.index(meta, present.year)
	if(is.null(present.year.index))
	    stop("present.year ", present.year, " not found in the data.")

	#keep these defaults for checking the out-of-sample projections
  quantiles.to.keep <- c(0,0.025,0.05,0.1,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.9,0.95,0.975,1)
	PIs_cqp <- array(NA, c(nr_countries, length(quantiles.to.keep), nr_project+1))
	dimnames(PIs_cqp)[[2]] <- quantiles.to.keep
	proj.middleyears <- bayesTFR:::get.prediction.years(meta, nr_project+1, 
	                                                    present.year.index = present.year.index)
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
	
    data.mtx.name <- if(post.last.observed != "i") "mig.rates.all" else "mig.rates"
    migrates <- meta[[data.mtx.name]]
    migrates.recon <- meta[["mig.rates.all"]]
    if(post.last.observed != "a" && present.year.index < ncol(migrates)) {
        migrates[, (present.year.index + 1): ncol(migrates)] <- NA
        migrates.recon[, (present.year.index + 1): ncol(migrates.recon)] <- NA
    }
    lmigrates <- ncol(migrates)
    allTend <- lmigrates
    
	# fill the first time point of the result array with observed data and find missing data
	nmissing <- all.data.list <- missing <- list()
	max.nr.project <- nr_project

    for(country in prediction.countries){
        #country.obj <- get.country.object(country, meta, index=TRUE)
        obs.rates <- migrates[country,]
        Tc <- max(which(!is.na(obs.rates)))
        Tend <- min(present.year.index, Tc)
        allTend <- min(allTend, Tend)
        if (post.last.observed != "i" && Tend < present.year.index) {
            while(Tend < present.year.index) { # shift the end as long as there are no NAs
                Tend <- Tend + 1
                Tc <-  Tc + 1
                if(is.na(obs.rates[Tend])){
                    Tend <- Tend - 1
                    Tc <-  Tc - 1
                    break
                }
            }
        }
        if(Tend < present.year.index)
            migrates.recon[country, (Tend + 1):present.year.index] <- NA
        
        nmissing[[country]] <- present.year.index - Tend
        missing[[country]] <- (Tend+1):lmigrates
        max.nr.project <- max(max.nr.project, nr_project + nmissing[[country]])
        all.data.list[[country]] <- migrates.recon[country,]

    }
	# array for results - includes also historical data for periods with missing data
	all.mig_ps <- array(NA, dim=c(nr_countries_real, max.nr.project + 1, nr_simu))
	fps.end.obs.index <- dim(migrates.recon)[2] - allTend + 1
	
	for (country in prediction.countries) {
	    for(year in 1:fps.end.obs.index) 
	        all.mig_ps[country, year,] = all.data.list[[country]][allTend + year-1]
	}
	mu.c <- phi.c <- sigma.c <- rep(NA, nr_countries)

	traj.counter <- 0
	country.loop.max <- 20
	if (verbose) {
		verbose.iter <- max(1, nr_simu/100)
		if(interactive()) cat('\n')
	}
	if(use.cummulative.threshold){
	    nperiods.for.threshold <- ifelse(meta$annual.simulation, 30, 6)
	    mig.thresholds <-  get.migration.thresholds(meta, nperiods = nperiods.for.threshold, ignore.gcc = ignore.gcc.in.threshold)
	    isGCC <- if(ignore.gcc.in.threshold) rep(FALSE, nr_countries_real) else is.gcc.plus(meta$regions$country_code)
	    fun.min <- ".min.multiplicative.pop.change"
	}
	fthresholds <- list(upper = rep(NA, nr_countries_real), lower = rep(NA, nr_countries_real))
	if(!is.null(fixed.thresholds)){
	    for(ttp in names(fthresholds)){
	        if(! ttp %in% names(fixed.thresholds)) next
	        for(cntry in names(fixed.thresholds[[ttp]])){
	           country.obj <- get.country.object(as.integer(cntry), meta)
	           fthresholds[[ttp]][country.obj$index] <- fixed.thresholds[[ttp]][[cntry]]
	        }
	    }
	}
	#########################################
	for (s in 1:nr_simu){ # Iterate over trajectories
	#########################################
		if (verbose) {
			if(interactive()) cat('\rProjecting migration trajectories ...', round(s/nr_simu * 100), '%')
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
	  for (icountry in 1:nr_countries_real){ # Iterate over countries
	  #########################################
	    if(use.cummulative.threshold) fun.max <- paste0(".max.multiplicative.pop.change", if(isGCC[icountry]) "" else ".no.gcc")
	    for (year in 2:(max.nr.project+1)) { # Iterate over time
	    #########################################
	        if(!is.na(all.mig_ps[icountry, year, s])) next
	        determ.part <- mu.c[icountry] + phi.c[icountry]*(all.mig_ps[icountry,year-1,s] - mu.c[icountry])
	        if(use.cummulative.threshold){
	            xmin <- .get.rate.mult.limit(all.mig_ps[icountry,1:(year-1),s], year-1, fun.min, max, nperiods=nperiods.for.threshold, thresholds = mig.thresholds)
	            xmax <- .get.rate.mult.limit(all.mig_ps[icountry,1:(year-1),s], year-1, fun.max, min, nperiods=nperiods.for.threshold, thresholds = mig.thresholds)
	            if(!is.na(fthresholds$lower[icountry]))
	                xmin <- min(xmin, fthresholds$lower[icountry])
	            if(!is.na(fthresholds$upper[icountry]))
	                xmax <- min(xmax, fthresholds$upper[icountry])
	            if(xmin > xmax) {
	                avg <- (xmin + xmax)/2.
	                xmin <- avg - 1e-3
	                xmax <- avg + 1e-3 
	            }
	            error <- rtruncnorm(n=1, a=xmin-determ.part, b=xmax-determ.part, mean=0,sd=sigma.c[icountry])
	        } else error <- rnorm(n=1, mean=0,sd=sigma.c[icountry])
	        all.mig_ps[icountry,year,s] <- determ.part + error
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

		# impute missing values if any
		if (nmissing[[country]] > 0) { # data imputation
		    mig_ps_future[1,] <- quantile(mig_ps_future[1,], 0.5, na.rm = TRUE) # set all trajectories in the first time period to the median
		    migrates.recon[country, (lmigrates - fps.end.obs.index + 2):lmigrates] <- apply(
		        all.mig_ps[icountry, 2:fps.end.obs.index, , drop=FALSE],
		        c(1,2), quantile, 0.5, na.rm = TRUE)
		    if (verbose) 
		        cat('\t', nmissing[[country]], 'data points reconstructed for', country.obj$name,'\n')
		}
		if(write.trajectories) {
			trajectories <- mig_ps_future # save only trajectories simulated for the future time
  			save(trajectories, file = file.path(outdir, paste('traj_country', country.obj$code, '.rda', sep='')))
  		}
  		# compute quantiles
 		PIs_cqp[country,,] <- apply(mig_ps_future, 1, quantile, quantiles.to.keep, na.rm = TRUE)
 		mean_sd[country,1,] <- apply(mig_ps_future, 1, mean, na.rm = TRUE)
 		mean_sd[country,2,] <- apply(mig_ps_future, 1, sd, na.rm = TRUE)
	}
	if (lmigrates > present.year.index)
	    migrates.recon[ , (present.year.index + 1):lmigrates] <- NA # does not need data beyond present.year as those values are now trajectories
	
	mcmc.set <- remove.mig.traces(mcmc.set)
	bayesMig.prediction <- structure(list(
				quantiles = PIs_cqp,
				traj.mean.sd = mean_sd,
				nr.traj=nr_simu,
				mig.rates.reconstructed = migrates.recon,
				nr.imputed = unlist(nmissing),
				output.directory=outdir,
				mcmc.set=load.mcmc.set,
				nr.projections=nr_project,
				burnin=burnin, thin=thin,
				end.year=end.year, present.year.index = present.year.index),
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

#' @title Converting Trajectories of Migration Rates into ACSII Files
#'
#' @description Converts trajectories of the net migration rates stored 
#'     in a binary format into two CSV files. 
#' 
#' @param sim.dir Directory containing the prediction object. It should be the same as
#'     the \code{output.dir} argument in \code{\link{mig.predict}}.
#' @param n Number of trajectories to be stored. It can be either a single number 
#'     or the word \dQuote{all} in which case all available trajectories are converted.
#'     If the number is smaller than the number of trajectories available 
#'     in the prediction object, they are selected by equal spacing.
#' @param output.dir Directory into which the resulting files will be stored. 
#'     If it is \code{NULL}, the same directory is used as for the prediction.
#' @param verbose Logical value. Switches log messages on and off.
#' 
#' @details The function creates two files. First, \dQuote{ascii_trajectories.csv}
#'     is a comma-separated table with the following columns: 
#'     \itemize{\item{\dQuote{LocID}: }{country code} 
#'         \item{\dQuote{Period}: }{prediction interval, e.g. 2015-2020} 
#'         \item{\dQuote{Year}: }{middle year of the prediction interval}
#'         \item{\dQuote{Trajectory}: }{identifier of the trajectory}
#'         \item{\dQuote{mig}: }{net migration rate}
#'     }
#'     The second file is called \dQuote{ascii_trajectories_wide.csv}, also 
#'     a comma-separated table and it contains the same information as above 
#'     but in a wide format. I.e. the data for one country are 
#'     ordered in columns, thus, there is one column per country. The country columns 
#'     are ordered alphabetically. 
#' @note This function is automatically called from the \code{\link{mig.predict}} 
#'     function, therefore in standard cases it will not be needed to call it directly. 
#'     However, it can be useful for example, if different number of trajectories are to be converted, 
#'     without having to re-run the prediction.
#' @return No return value.
#' @seealso \code{\link[bayesTFR]{convert.tfr.trajectories}}, 
#'   \code{\link{mig.write.projection.summary}}, \code{\link{get.mig.trajectories}}
#' @export
#' 
convert.mig.trajectories <- function(sim.dir = NULL, n = 1000, output.dir = NULL, verbose=FALSE) {
    # Converts all trajectory rda files into UN ascii, selecting n trajectories by equal spacing.
    if(n <= 0) return()
    pred <- get.mig.prediction(sim.dir)
    outdir <- if (is.null(output.dir)) pred$output.directory else output.dir
    bayesTFR:::do.convert.trajectories(pred=pred, n=n, output.dir=outdir, verbose=verbose)
}

remove.mig.traces <- function(mcmc.set) {
	for (i in 1:length(mcmc.set$mcmc.list))
		mcmc.set$mcmc.list[[i]]$traces <- 0
	invisible(mcmc.set)
}

#' @export
get.traj.ascii.header.bayesMig.mcmc.meta <- function(meta, ...) 
	return (list(country_code='LocID', period='Period', year='Year', trajectory='Trajectory', tfr='Mig'))

#' @title Writing Projection Summary Files
#' @description The function creates two files containing projection summaries, 
#'     such as the median, the lower and upper bound of the 80 and 90\% probability 
#'     intervals, respectively, and the constant variant. One file is in a user-friendly 
#'     format, whereas the other is in a UN-specific format with internal coding of the 
#'     time and the variants.
#' @param pred Object of class \code{bayesMig.prediction}.
#' @param output.dir Directory where output is written.
#' @param \dots Additional arguments passed to the underlying functions. 
#'     Here, argument \code{precision} can be set to determine the number 
#'     of significant digits (default is 4).
#' @return No return value.
#' @seealso \code{\link[bayesTFR]{write.projection.summary}}
#' @export
mig.write.projection.summary <- function(pred, output.dir, ...) {
	# one summary file
    bayesTFR:::do.write.projection.summary(pred, output.dir, ...)
}

get.estimation.years <- function(meta)
  return(as.numeric(colnames(meta$mig.rates))+2.5)

#' @export
get.projection.summary.header.bayesMig.prediction <- function(pred, ...) 
  return (list(revision='RevID', variant='VarID', country='LocID', year='TimeID', indicator='IndicatorID', sex='SexID', tfr='Value'))

#' @export
get.friendly.variant.names.bayesMig.prediction <- function(pred, ...)
  return(c('median', 'lower 80', 'upper 80', 'lower 95', 'upper 95', 'mean', 'constant'))

#' @export
get.UN.variant.names.bayesMig.prediction <- function(pred, ...) 
    return(c('BHM median', 'BHM80 lower',  'BHM80 upper', 'BHM95 lower',  'BHM95 upper', 'BHM mean', 'Zero migration'))


get.mig.periods <- function(meta) {
  mid.years <- get.estimation.years(meta)
  return (paste(mid.years-2.5, mid.years+2.5, sep='-'))
}

#' @export
get.data.imputed.bayesMig.prediction <- function(pred, ...)
    return(t(pred$mig.rates.reconstructed))

#' @export
get.data.for.country.imputed.bayesMig.prediction <- function(pred, country.index, ...)
    return(get.data.imputed(pred)[, country.index])

.max.multiplicative.pop.change <- function(l, thresholds) 
    thresholds$upper[l]

.max.multiplicative.pop.change.no.gcc <- function(l, thresholds) 
    thresholds$upper.nogcc[l]

.min.multiplicative.pop.change <- function(l, thresholds) 
    thresholds$lower[l]

.get.rate.mult.limit <- function(rates, n, cumfun, fun, nperiods=6, ...) {
    res <- do.call(cumfun, list(1, ...))
    for(i in 2:min(nperiods,n+1)) {
        p <- prod(1+rates[(n-i+2):n])
        res <- c(res, do.call(cumfun, c(list(i, ...)))/p)
    }
    return(do.call(fun, list(res))-1)
}

is.gcc.plus <- function(country) # GCC plus Western Sahara & Djibouti
    return(country %in% c(634, 784, 414, 48, 512, 682, 732, 262)) # Qatar, UAE, Kuwait, Bahrain, Oman, SA, Western Sahara, Djibouti

get.migration.thresholds <- function(meta, nperiods=6, ignore.gcc = FALSE) {
    # Setting cummulative thresholds
    rates <- meta$mig.rates.all
    # GCC plus Western Sahara & Djibouti
    if(!ignore.gcc)
        gcc.plus <- meta$regions$country_code[is.gcc.plus(meta$regions$country_code)]
    
    rMat <- 1 + rates
    tu <- apply(rMat, 1, max, na.rm = TRUE)
    tl <- apply(rMat, 1, min, na.rm = TRUE)
    for (i in 2:min(nperiods, ncol(rates)-1)) {
        p <- 0*rates[,1:(ncol(rates)-i+1)] + 1 # init with 1
        for(j in 1:i) 	
            p <- p * rMat[,j:(ncol(rates)-i+j)]
        tu <- cbind(tu, apply(p, 1, max, na.rm = TRUE))
        tl <- cbind(tl, apply(p, 1, min, na.rm = TRUE))
    }
    upper.bounds <- apply(tu, 2, max, na.rm = TRUE)
    upper.bounds.nogcc <- if(!ignore.gcc) apply(tu[!rownames(tu) %in% gcc.plus,], 2, max, na.rm = TRUE) else rep(NA, length(upper.bounds))
    lower.bounds <- apply(tl, 2, min, na.rm = TRUE)
    
    df <- data.frame(upper=upper.bounds, upper.nogcc=upper.bounds.nogcc, lower=lower.bounds)
    return(df)
}


#' @title Adjusting the Projection Medians
#' 
#' @description These functions are to be used by expert analysts. They allow to 
#'     change the projection medians either to specific values, or shift the medians 
#'     by a given constant or align one projection object with another.
#'     
#' @param sim.dir Directory containing the prediction object.
#' @param country Name or numerical code of a country.
#' @param values Vector of the new median values.
#' @param years Numeric vector giving years for which to change the median. 
#'     In \code{mig.median.set} it gives years which \code{values} correspond to.
#'     Ideally it should be of the same length as \code{values}. If it is \code{NULL}, 
#'     \code{values} are set starting from the first prediction time period. 
#'     If \code{values} correspond to consecutive years, only the first year might be given here.
#'     In \code{mig.align.predictions} it gives years for which the medians should be aligned.  
#' @param \dots Additional arguments passed to the underlying adjustment function. 
#'     It can be \code{verbose} to show/hide the progress of the adjustment.
#' 
#' @details The function \code{mig.median.set} can be used to set the medians of the
#'     given country to specific values.
#'    
#' @return All functions return an updated object of class \code{\link{bayesMig.prediction}}.
#' @export
#' @rdname mig.adjust
mig.median.set <- function(sim.dir, country, values, years=NULL, ...) {
    pred <- get.mig.prediction(sim.dir)
    new.pred <- bayesTFR:::.bdem.median.set(pred, type='mig', country=country, 
                                            values=values, years=years, ...)
    store.bayesMig.prediction(new.pred)
    invisible(new.pred)
}

#' @export 
#' @keywords internal
#' @rdname internal
#' 
get.mig.shift <- function(country.code, pred) return(bayesTFR::get.tfr.shift(country.code, pred))

#' @param reset Logical. If \code{TRUE} medians in a range of \code{from} and \code{to} are
#'     reset to their original values.
#' @param shift Constant by which the medians should be offset. It is not used if \code{reset} is \code{TRUE}.
#' @param from Year from which the offset/reset should start. By default, it starts at the first prediction period.
#' @param to Year until which the offset/reset should be done. By default, it is set to the last prediction period.
#' 
#' @details Function \code{mig.median.shift} can be used to offset the medians by a specific constant, or to reset
#'     the medians to their original values. 
#' @export
#' @rdname mig.adjust
mig.median.shift <- function(sim.dir, country, reset = FALSE, shift = 0, 
                             from = NULL, to = NULL) {
    pred <- get.mig.prediction(sim.dir)
    new.pred <- bayesTFR:::.bdem.median.shift(pred, type='mig', country=country, reset=reset, 
                                              shift=shift, from=from, to=to)
    store.bayesMig.prediction(new.pred)
    invisible(new.pred)
}

#' @param countries Vector of country names or codes. If this argument is \code{NULL} (default), 
#'     the reset is done for all countries.
#'     
#' @details Function \code{mig.median.reset} resets medians of the given countries
#'     to the original values. By default it deletes adjustments for all countries.
#' @export
#' @rdname mig.adjust
mig.median.reset <- function(sim.dir, countries = NULL) {
    if(is.null(countries)) {
        pred <- get.mig.prediction(sim.dir)
        pred$median.shift <- NULL
        store.bayesMig.prediction(pred)
        cat('\nMedians for all countries reset.\n')
    } else
        for(country in countries) pred <- mig.median.shift(sim.dir, country, reset=TRUE)
    invisible(pred)
}

#' @param sim.dir1 Directory with the bayesMig prediction object to be adjusted.
#' @param sim.dir2 Directory with the bayesMig prediction object used to align the medians from \code{sim.dir1} to.
#' @param country.codes Numerical codes of countries to adjust. By default all countries 
#'     found in \code{sim.dir2} are adjusted in \code{sim.dir1}.
#' @details Function \code{mig.align.predictions} shifts medians stored in \code{sim.dir1} to match 
#'     the medians in \code{sim.dir1}.
#'     
#'     In all cases, if a median is modified, the corresponding offset is stored in the prediction object 
#'     (element \code{median.shift}). All functions write the updated prediction object back to disk. All
#'     functions in the package that use trajectories and trajectory statistics use the \code{median.shift} 
#'     values to offset the results correspondingly, i.e. trajectories are shifted the same way as the
#'     medians.
#' @export
#' @rdname mig.adjust
mig.align.predictions <- function(sim.dir1, sim.dir2, country.codes = NULL, 
                                  years = NULL, ...){
    pred1 <- get.mig.prediction(sim.dir1)
    pred2 <- get.mig.prediction(sim.dir2)
    cntries1 <- if(is.null(country.codes)) get.countries.table(pred1)$code else country.codes
    cntries2 <- get.countries.table(pred2)$code
    avail.years <- dimnames(pred2$quantiles)[[3]]
    if(is.null(years)) 
        years <- avail.years
    else years <- intersect(as.character(years), avail.years)
    for(cntry in cntries1){
        cidx <- which(cntries2 == cntry)
        if(length(cidx) == 0){
            warning("Country ", cntry, " not found in sim.dir2. No adjustment for this country.")
            next
        }
        adjust.to <- get.median.from.prediction(pred2, cidx, cntry)[years]
        pred <- mig.median.set(sim.dir1, cntry, values = adjust.to, 
                               years = as.integer(names(adjust.to)), ...)
    }
    invisible(pred)
}