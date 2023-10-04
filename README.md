# bayesMig

[![R-CMD-check](https://github.com/PPgp/bayesMig/actions/workflows/check-standard.yaml/badge.svg?branch=cran&event=push)](https://github.com/PPgp/bayesMig/actions/workflows/check-standard.yaml)

R package for probabilistic projections of net migration rate for all countries of the world or for subnational units, using a Bayesian hierarchical model by [Azose and Raftery (2015)](https://doi.org/10.1007/s13524-015-0415-0).

### Installation

```
library(devtools)
install_github("PPgp/bayesMig")

```

### Usage 

The two main functions of the package are:

* `run.mig.mcmc`: Runs a Markov Chain Monte Carlo (MCMC) simulation. It results in posterior samples of the model parameters.
* `mig.predict`: Using the posterior parameter samples, trajectories of future net migration rates are generated for all countries or given locations.

See `?bayesMig` for more info.


#### Example
Run a simulation that uses default 5-year national historical data from WPP 2019, ranging from 1950 to 2020, and projects to 2100.

```
library(bayesMig)

sim.dir <- tempfile() # directory to store results into

# Run 4 MCMCs 10,000 iterations each, thinned by 10
# (can take long time; reduce iter for a toy simulation)
m <- run.mig.mcmc(nr.chains = 4, iter = 10000, 
		thin = 10, output.dir = sim.dir, verbose.iter = 1000)

# Prediction for all countries
pred <- mig.predict(sim.dir = sim.dir, nr.traj = 1000, burnin = 5000)

# Explore results
summary(pred, country = "Germany")
mig.trajectories.plot(pred, country = "Germany", nr.traj = 50)

# Retrieve simulation objects from disk, e.g. at later time
m <- get.mig.mcmc(sim.dir)
pred <- get.mig.prediction(sim.dir)

# Remove simulation directory when not needed
unlink(sim.dir, recursive = TRUE)

```