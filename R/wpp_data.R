# Read in the UN estimates

load.bdem.dataset <- function(dataset, wpp.year, envir=NULL, verbose=FALSE) {
	pkg <- paste('wpp', wpp.year, sep='')
	do.call('require', list(pkg))
	if(verbose) cat('Loading ', dataset, ' from ', pkg, '.\n')
	if(is.null(envir)) envir <- new.env()
	do.call('data', list(dataset, package=pkg, envir=envir))
	return(envir[[dataset]])
}
