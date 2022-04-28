start.test <- function(name) cat('\n<=== Starting test of', name,'====\n')
test.ok <- function(name) cat('\n==== Test of', name, 'OK.===>\n')

test.run.annual.simulation <- function(parallel = FALSE) {
    sim.dir <- tempfile()
    
    # run MCMC
    test.name <- 'running annual migration MCMC'
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
    sim.dir <- tempfile()
    
    # run MCMC
    test.name <- 'running national migration MCMC'
    if(parallel) test.name <- paste(test.name, "(in parallel)")
    start.test(test.name)

    sim.dir <- tempfile()
    m <- run.mig.mcmc(nr.chains = 2, iter = 60, thin = 2, output.dir = sim.dir, parallel = parallel)
    stopifnot(m$mcmc.list[[1]]$finished.iter == 60)
    stopifnot(get.total.iterations(m$mcmc.list, 0) == 120)
    test.ok(test.name)
    
    # Prediction
    test.name <- 'running national projections'
    start.test(test.name)
    pred <- mig.predict(sim.dir = sim.dir, burnin = 10, end.year = 2050)
    spred <- summary(pred)
    stopifnot(spred$nr.traj == 50)
    stopifnot(nrow(get.countries.table(pred))== 200)
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