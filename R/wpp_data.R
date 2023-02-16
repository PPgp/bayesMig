# Read in the UN estimates

load.from.wpp <- function(dataset, wpp.year, annual = FALSE, ...){
    if(wpp.year >= 2022) {
        dsname <- paste0(dataset, if(annual) 1 else 5)
        ds <- load.bdem.dataset(dsname, wpp.year, ..., check.if.exists = TRUE)
        if(!is.null(ds)) return(ds)
        # otherwise check the original name
    }
    load.bdem.dataset(dataset, wpp.year, ...)
}


load.bdem.dataset <- function(dataset, wpp.year, envir=NULL, annual = FALSE, verbose=FALSE) {
	pkg <- paste('wpp', wpp.year, sep='')
	do.call('require', list(pkg))
	if(verbose) cat('Loading ', dataset, ' from ', pkg, '.\n')
	if(is.null(envir)) envir <- new.env()
	do.call('data', list(dataset, package=pkg, envir=envir))
	return(envir[[dataset]])
}
