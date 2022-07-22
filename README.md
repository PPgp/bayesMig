# bayesMig

[![R build status](https://github.com/PPgp/bayesMig/workflows/R-CMD-check/badge.svg)](https://github.com/PPgp/bayesMig/actions?workflow=R-CMD-check)

R package for probabilistic projections of net migration rate for all countries of the world or for subnational units, using a Bayesian hierarchical model by [Azose and Raftery (2015)](https://doi.org/10.1007/s13524-015-0415-0).

The two main functions of the package are:

* `run.mig.mcmc`: Runs a Markov Chain Monte Carlo (MCMC) simulation. It results in posterior samples of the model parameters.
* `mig.predict`: Using the posterior parameter samples, trajectories of future net migration rates are generated for all countries or given locations.

See `?bayesMig` for more info.

Example:

```
# Run a realistic simulation 
# (can take long time; reduce iter for a toy simulation)

sim.dir <- tempfile() # directory to store results into

# Run 4 MCMCs 10,000 iterations each, thinned by 10
# (uses default national data from WPP 2019)
m <- run.mig.mcmc(nr.chains = 4, iter = 10000, 
		thin = 10, output.dir = sim.dir, verbose.iter = 1000)

# Prediction for all countries
pred <- mig.predict(sim.dir = sim.dir, nr.traj = 1000, burnin = 1000)

# Explore results
summary(pred, country = "Germany")
mig.trajectories.plot(pred, country = "Germany", nr.traj = 50)


```