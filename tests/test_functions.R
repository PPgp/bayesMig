start.test <- function(name) cat('\n<=== Starting test of', name,'====\n')
test.ok <- function(name) cat('\n==== Test of', name, 'OK.===>\n')

test.run.annual.simulation <- function(parallel = FALSE) {
    # run MCMC
    test.name <- 'running annual migration MCMC for US states'
    if(parallel) test.name <- paste(test.name, "(in parallel)")
    start.test(test.name)
    us.mig.file <- file.path(find.package("bayesMig"), "extdata", "USmigrates.txt")
    sim.dir <- tempfile()
    m <- run.mig.mcmc(nr.chains = 2, iter = 30, thin = 1, my.mig.file = us.mig.file, 
             output.dir = sim.dir, present.year = 2017, annual = TRUE, parallel = parallel)
    stopifnot(m$mcmc.list[[1]]$finished.iter == 30)
    stopifnot(get.total.iterations(m$mcmc.list, 0) == 60)
    
    par.values <- get.mig.parameter.traces(m$mcmc.list, burnin = 5)
    stopifnot(all(dim(par.values) == c(50, 4)))
    par.values.cs <- get.mig.parameter.traces.cs(m$mcmc.list, 
                        country.obj = get.country.object("California", meta = m$meta),
                        burnin = 5, par.names = "phi_c")
    stopifnot(all(dim(par.values.cs) == c(50, 1)))
    test.ok(test.name)
    
    # Prediction
    test.name <- 'running annual projections'
    start.test(test.name)
    pred <- mig.predict(sim.dir = sim.dir, burnin = 10, end.year = 2050)
    spred <- summary(pred)
    stopifnot(spred$nr.traj == 40)
    stopifnot(nrow(get.countries.table(pred))== 52)
    stopifnot(dim(pred$quantiles)[3] == length(2017:2050))
    test.ok(test.name)
    
    # output
    test.name <- 'analyzing output of annual projections'
    start.test(test.name)
    tab <- mig.trajectories.table(pred, "Hawaii")
    years <- as.integer(rownames(tab))
    should.be.years <- 2001:2050
    stopifnot(length(years) == length(should.be.years))
    stopifnot(all(years == should.be.years))
    stopifnot(all(dim(tab) == c(length(should.be.years), 5)))
    test.ok(test.name)
    
    unlink(sim.dir, recursive=TRUE)
}

test.run.national.simulation <- function(parallel = FALSE) {
    # run MCMC
    test.name <- 'running national migration MCMC'
    if(parallel) test.name <- paste(test.name, "(in parallel)")
    start.test(test.name)

    sim.dir <- tempfile()
    m <- run.mig.mcmc(nr.chains = 2, iter = 60, thin = 2, output.dir = sim.dir, parallel = parallel,
                      wpp.year = 2019)
    stopifnot(m$mcmc.list[[1]]$finished.iter == 60)
    stopifnot(get.total.iterations(m$mcmc.list, 0) == 120)
    test.ok(test.name)
    
    # Prediction
    test.name <- 'running national projections'
    start.test(test.name)
    pred <- mig.predict(sim.dir = sim.dir, burnin = 10, end.year = 2050)
    spred <- summary(pred)
    stopifnot(spred$nr.traj == 50)
    stopifnot(nrow(get.countries.table(pred))== 201)
    stopifnot(dim(pred$quantiles)[3] == length(seq(2018, 2048, by = 5)))
    test.ok(test.name)
    
    # output
    test.name <- 'analyzing output of national projections'
    start.test(test.name)
    tab <- mig.trajectories.table(pred, "France")
    years <- as.integer(rownames(tab))
    should.be.years <- seq(1953, 2048, by = 5)
    stopifnot(length(years) == length(should.be.years))
    stopifnot(all(years == should.be.years))
    test.ok(test.name)
    
    unlink(sim.dir, recursive=TRUE)
}

test.run.annual.national.simulation <- function(parallel = FALSE) {
    # run MCMC using wpp2022
    test.name <- 'running annual national migration MCMC'
    if(parallel) test.name <- paste(test.name, "(in parallel)")
    start.test(test.name)
    
    sim.dir <- tempfile()
    
    # find small countries to be excluded (take it from bayesTFR include dataset)
    data(include_2022, package = "bayesTFR")
    small.countries <- subset(include_2022, include_code == 1)$country_code
    
    m <- run.mig.mcmc(nr.chains = 2, iter = 60, thin = 2, output.dir = sim.dir, 
                      parallel = parallel, annual = TRUE, wpp.year = 2022,
                      present.year = 2021, exclude.from.world = small.countries,
                      use.cummulative.threshold = TRUE)
    
    stopifnot(m$meta$nr.countries.est == 203)
    stopifnot(m$meta$nr.countries == 236)
    stopifnot(m$mcmc.list[[1]]$finished.iter == 60)
    stopifnot(get.total.iterations(m$mcmc.list, 0) == 120)
    test.ok(test.name)
    
    # Prediction
    test.name <- 'running annual national projections'
    start.test(test.name)
    pred <- mig.predict(sim.dir = sim.dir, burnin = 10, end.year = 2050)
    spred <- summary(pred)
    stopifnot(spred$nr.traj == 50)
    stopifnot(nrow(get.countries.table(pred))== 236)
    stopifnot(dim(pred$quantiles)[3] == length(seq(2021, 2050, by = 1)))
    test.ok(test.name)
    
    # output
    test.name <- 'analyzing output of annual national projections'
    start.test(test.name)
    tab <- mig.trajectories.table(pred, "Ireland")
    years <- as.integer(rownames(tab))
    should.be.years <- seq(1951, 2050, by = 1)
    stopifnot(length(years) == length(should.be.years))
    stopifnot(all(years == should.be.years))
    test.ok(test.name)
    
    unlink(sim.dir, recursive=TRUE)
}


test.run.annual.national.simulation.with.interpolation <- function(parallel = FALSE) {
    # run MCMC
    test.name <- 'running annual national migration MCMC with interpolated data'
    if(parallel) test.name <- paste(test.name, "(in parallel)")
    start.test(test.name)
    
    sim.dir <- tempfile()
    m <- run.mig.mcmc(nr.chains = 2, iter = 60, thin = 2, output.dir = sim.dir, 
                      parallel = parallel, annual = TRUE, wpp.year = 2019)
    
    stopifnot(m$mcmc.list[[1]]$finished.iter == 60)
    stopifnot(get.total.iterations(m$mcmc.list, 0) == 120)
    test.ok(test.name)
    
    # Prediction
    test.name <- 'running annual national projections with interpolated data'
    start.test(test.name)
    pred <- mig.predict(sim.dir = sim.dir, burnin = 10, end.year = 2050)
    spred <- summary(pred)
    stopifnot(spred$nr.traj == 50)
    stopifnot(nrow(get.countries.table(pred))== 201)
    stopifnot(dim(pred$quantiles)[3] == length(seq(2020, 2050, by = 1)))
    test.ok(test.name)
    
    # output
    test.name <- 'analyzing output of annual national projections with interpolated data'
    start.test(test.name)
    tab <- mig.trajectories.table(pred, "Ireland")
    years <- as.integer(rownames(tab))
    should.be.years <- seq(1950, 2050, by = 1)
    stopifnot(length(years) == length(should.be.years))
    stopifnot(all(years == should.be.years))
    test.ok(test.name)
    
    unlink(sim.dir, recursive=TRUE)
}

test.include.code.and.last.observed <- function(parallel = FALSE) {
    # run MCMC
    test.name <- 'running annual migration MCMC with states excluded and missing data'
    if(parallel) test.name <- paste(test.name, "(in parallel)")
    start.test(test.name)
    us.mig.file <- file.path(find.package("bayesMig"), "extdata", "USmigrates.txt")
    mig <- bayesTFR:::read.tfr.file(file = us.mig.file)
    mig$include_code <- 2
    mig$last.observed <- 2017
    mig[mig$name %in% c("Rhode Island", "District of Columbia"), "include_code"] <- 1 # used only for prediction
    mig[mig$name == "Hawaii", "include_code"] <- 0 # excluded
    mig[mig$name == "Washington", "last.observed"] <- 2015 # 2 data points will be imputed
    mig[mig$name == "Idaho", "2017"] <- NA # 1 data point missing without changing last.observed
    migfile <- tempfile()
    write.table(mig, file = migfile, sep='\t', row.names=FALSE)
    
    sim.dir <- tempfile()
    m <- run.mig.mcmc(nr.chains = 2, iter = 30, thin = 1, my.mig.file = migfile, 
                      output.dir = sim.dir, present.year = 2017, annual = TRUE, parallel = parallel,
                      exclude.from.world = 29) # also exclude Nevada

    stopifnot((m$meta$nr.countries - m$meta$nr.countries.est) == 3) # 3 countries excluded from estimation
    stopifnot(! "Hawaii" %in% m$meta$regions$country_name) # Hawaii is not included at all
    stopifnot(all(is.na(m$meta$mig.rates[m$meta$regions$country_name == "Washington", c("2016", "2017")]))) # missing data
    test.ok(test.name)
    
    # Prediction
    test.name <- 'running annual projections with data imputation'
    start.test(test.name)
    pred <- mig.predict(sim.dir = sim.dir, burnin = 10, end.year = 2050)
    
    imputed <- pred$mig.rates.reconstructed[m$meta$regions$country_name == "Washington", c("2016", "2017")]
    orig <- m$meta$mig.rates.all[m$meta$regions$country_name == "Washington", c("2016", "2017")]
    stopifnot(all(!is.na(imputed))) # was it imputed
    stopifnot(all.equal(imputed, orig)  > 0.1) # relative difference is big
    stopifnot(nrow(get.countries.table(pred))== 51) # 51 states included
    stopifnot(!is.na(pred$mig.rates.reconstructed[m$meta$regions$country_name == "Idaho", "2017"])) # Idaho was imputed
    stopifnot(dim(pred$quantiles)[3] == length(2017:2050))
    test.ok(test.name)
    
    unlink(migfile)
    unlink(sim.dir, recursive=TRUE)
}

test.adjustments <- function() {
    sim.dir <- tempfile()
    us.mig.file <- file.path(find.package("bayesMig"), "extdata", "USmigrates.txt")
    m <- run.mig.mcmc(nr.chains = 1, iter = 30, thin = 1, my.mig.file = us.mig.file, 
                      output.dir = sim.dir, present.year = 2017, annual = TRUE, verbose = FALSE)
    pred <- mig.predict(sim.dir = sim.dir, burnin = 10, end.year = 2030, verbose = FALSE)
    projs <- summary(pred, country = 'Hawaii')$projections
    shifted.cols <- c(1, 3:ncol(projs)) # all but SD
    years <- rownames(projs)
    
    test.name <- 'shifting the trajectories'
    start.test(test.name)
    mig.traj.shift(sim.dir, country = 'Hawaii', shift = 0.01, from = 2020, to = 2025)
    shifted.pred <- get.mig.prediction(sim.dir)
    shifted.projs <- summary(shifted.pred, country = 'Hawaii')$projections
    sidx <- years %in% as.character(2020:2025)
    stopifnot(all.equal(projs[sidx, shifted.cols] + 0.01, shifted.projs[sidx, shifted.cols]))
    stopifnot(all(projs[!sidx, shifted.cols] == shifted.projs[!sidx, shifted.cols]))
    stopifnot(all(projs[, 2] == shifted.projs[, 2])) # SD does not change
    test.ok(test.name)
    
    test.name <- 'resetting the trajectories'
    start.test(test.name)
    shifted.pred <- mig.traj.shift(sim.dir, country = 'Hawaii', reset = TRUE)
    shifted.projs <- summary(shifted.pred, country = 'Hawaii')$projections
    stopifnot(all(projs[, shifted.cols] == shifted.projs[, shifted.cols]))
    stopifnot(is.null(get.mig.shift(get.country.object('Hawaii', m$meta)$code, shifted.pred)))
    test.ok(test.name)
    
    test.name <- 'setting the median'
    start.test(test.name)
    expert.values <- c(0.01, 0.015, 0.02)
    cobj <- get.country.object('Hawaii', m$meta)
    sidx <- years %in% as.character(2020:2022)
    shift <- expert.values - pred$quantiles[cobj$index, '0.5', sidx]
    mod.pred <- mig.median.set(sim.dir, country = 'Hawaii', values = expert.values, years = 2020)
    mod.projs <- summary(mod.pred, country = 'Hawaii')$projections
    stopifnot(all.equal(mod.projs[sidx, "50%"], expert.values, check.attributes = FALSE))
    stopifnot(all.equal(mod.projs[sidx, shifted.cols], projs[sidx, shifted.cols] + shift))
    stopifnot(all(mod.projs[!sidx, shifted.cols] == projs[!sidx, shifted.cols]))
    test.ok(test.name)
    
    test.name <- 'setting the mean'
    start.test(test.name)
    shift <- expert.values - pred$traj.mean.sd[cobj$index, 1, sidx]
    mig.shift.reset(sim.dir, countries = 'Hawaii') # reset first
    mod.pred <- mig.mean.set(sim.dir, country = 'Hawaii', values = expert.values, years = 2020)
    mod.projs <- summary(mod.pred, country = 'Hawaii')$projections
    stopifnot(all.equal(mod.projs[sidx, "mean"], expert.values, check.attributes = FALSE))
    stopifnot(all.equal(mod.projs[sidx, shifted.cols], projs[sidx, shifted.cols] + shift))
    stopifnot(all(mod.projs[!sidx, shifted.cols] == projs[!sidx, shifted.cols]))
    # the mean of the adjusted trajectories matches as well
    traj <- get.mig.trajectories(mod.pred, country = 'Hawaii')
    stopifnot(all.equal(rowMeans(traj)[as.character(2020:2022)], expert.values, check.attributes = FALSE))
    test.ok(test.name)
    
    test.name <- 'aligning predictions'
    start.test(test.name)
    # a second, independent simulation to be aligned with the (mean-adjusted) first one
    sim.dir2 <- tempfile()
    run.mig.mcmc(nr.chains = 1, iter = 30, thin = 1, my.mig.file = us.mig.file, 
                 output.dir = sim.dir2, present.year = 2017, annual = TRUE, verbose = FALSE)
    pred2 <- mig.predict(sim.dir = sim.dir2, burnin = 10, end.year = 2030, verbose = FALSE)
    projs2 <- summary(pred2, country = 'Hawaii')$projections
    aligned.pred <- mig.align.predictions(sim.dir2, sim.dir, country.codes = cobj$code, verbose = FALSE)
    aligned.projs <- summary(aligned.pred, country = 'Hawaii')$projections
    stopifnot(all.equal(aligned.projs[, "50%"], mod.projs[, "50%"]))
    stopifnot(!isTRUE(all.equal(aligned.projs[sidx, "mean"], mod.projs[sidx, "mean"]))) # means are not aligned
    test.ok(test.name)
    
    test.name <- 'aligning predictions by means'
    start.test(test.name)
    mig.shift.reset(sim.dir2)
    aligned.pred <- mig.align.predictions(sim.dir2, sim.dir, country.codes = cobj$code, 
                                          stat = "mean", verbose = FALSE)
    aligned.projs <- summary(aligned.pred, country = 'Hawaii')$projections
    stopifnot(all.equal(aligned.projs[, "mean"], mod.projs[, "mean"]))
    stopifnot(!isTRUE(all.equal(aligned.projs[sidx, "50%"], mod.projs[sidx, "50%"]))) # medians are not aligned
    # align only selected years
    mig.shift.reset(sim.dir2)
    aligned.pred <- mig.align.predictions(sim.dir2, sim.dir, country.codes = cobj$code, 
                                          years = 2021:2022, stat = "mean", verbose = FALSE)
    aligned.projs <- summary(aligned.pred, country = 'Hawaii')$projections
    aidx <- years %in% as.character(2021:2022)
    stopifnot(all.equal(aligned.projs[aidx, "mean"], mod.projs[aidx, "mean"]))
    stopifnot(all(aligned.projs[!aidx, "mean"] == projs2[!aidx, "mean"]))
    unlink(sim.dir2, recursive = TRUE)
    test.ok(test.name)
    
    test.name <- 'resetting all countries'
    start.test(test.name)
    mig.traj.shift(sim.dir, country = 'Alaska', shift = 0.01)
    stopifnot(length(get.mig.prediction(sim.dir)$traj.shift) == 2)
    mig.shift.reset(sim.dir)
    new.pred <- get.mig.prediction(sim.dir)
    stopifnot(is.null(new.pred$traj.shift))
    test.ok(test.name)
    
    test.name <- 'converting old median.shift'
    start.test(test.name)
    old.shift <- rep(0.01, dim(new.pred$quantiles)[3])
    new.pred$median.shift <- list()
    new.pred$median.shift[[as.character(cobj$code)]] <- old.shift
    store.bayesMig.prediction(new.pred)
    conv.pred <- get.mig.prediction(sim.dir)
    stopifnot(is.null(conv.pred$median.shift))
    stopifnot(all(get.mig.shift(cobj$code, conv.pred) == old.shift))
    test.ok(test.name)
    
    unlink(sim.dir, recursive = TRUE)
}

test.shift.to.wpp <- function(wpp.year = 2024) {
    sim.dir <- tempfile()
    m <- run.mig.mcmc(nr.chains = 1, iter = 30, thin = 1, output.dir = sim.dir, 
                      wpp.year = 2019, verbose = FALSE)
    pred <- mig.predict(sim.dir = sim.dir, burnin = 10, end.year = 2050, verbose = FALSE)
    
    # WPP rates for Mexico
    e <- new.env()
    data("migproj5dt", package = paste0("wpp", wpp.year), envir = e)
    data("popproj5dt", package = paste0("wpp", wpp.year), envir = e)
    wppmig <- data.table::data.table(e$migproj5dt)[country_code == 484]
    wpppop <- data.table::data.table(e$popproj5dt)[country_code == 484]
    wppmig[, year := year + 2] # align migration and pop years
    wpp <- merge(wppmig[, c("year", "mig"), with = FALSE], wpppop[, c("year", "pop"), with = FALSE], by = "year")
    wpp[, rate := mig / (pop - mig)][, year := year - 2]
    
    test.name <- 'shifting medians to WPP'
    start.test(test.name)
    shifted.pred <- mig.shift.prediction.to.wpp(sim.dir, wpp.year = wpp.year, verbose = FALSE)
    shifted.projs <- summary(shifted.pred, country = 'Mexico')$projections
    dat <- merge(wpp, data.table::data.table(year = as.integer(rownames(shifted.projs)), 
                                             median = shifted.projs[, "50%"]), by = "year")
    stopifnot(nrow(dat) > 0)
    stopifnot(all.equal(dat$rate, dat$median))
    stopifnot(length(shifted.pred$traj.shift) > 0)
    stopifnot(is.null(shifted.pred$median.shift))
    test.ok(test.name)
    
    test.name <- 'shifting means to WPP'
    start.test(test.name)
    shifted.pred <- mig.shift.prediction.to.wpp(sim.dir, wpp.year = wpp.year, stat = "mean", verbose = FALSE)
    shifted.projs <- summary(shifted.pred, country = 'Mexico')$projections
    dat <- merge(wpp, data.table::data.table(year = as.integer(rownames(shifted.projs)), 
                                             mean = shifted.projs[, "mean"]), by = "year")
    stopifnot(nrow(dat) > 0)
    stopifnot(all.equal(dat$rate, dat$mean))
    test.ok(test.name)
    
    unlink(sim.dir, recursive = TRUE)
}
