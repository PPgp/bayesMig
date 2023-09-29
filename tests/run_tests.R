library(bayesMig)
source('test_functions.R')

#options(error=quote(dump.frames("last.dump", TRUE)))

cran <- TRUE

test.run.annual.simulation()

if(!cran) {
    test.include.code.and.last.observed()
    test.run.annual.national.simulation()
    test.run.annual.national.simulation.with.interpolation()
    test.run.national.simulation()
    test.run.annual.simulation(parallel = TRUE)
    test.run.national.simulation(parallel = TRUE)
}

#load("last.dump.rda"); debugger()