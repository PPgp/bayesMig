library(bayesMig)
source('test_functions.R')

options(error=quote(dump.frames("last.dump", TRUE)))

cran <- FALSE

test.run.annual.simulation()

if(!cran) {
    test.run.national.simulation()
    test.run.annual.simulation(parallel = TRUE)
    test.run.national.simulation(parallel = TRUE)
}

#load("last.dump.rda"); debugger()