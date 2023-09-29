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
