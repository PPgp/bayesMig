if(getRversion() >= "2.15.1") utils::globalVariables(c("counter"))

###############
#MIGRATION
###############


store.mcmc <- local({
  # Writes parameter values into ascii files - one file per parameter and country (if country-specific)
  ##########################
  par.names <- mig.parameter.names()#Parameter names (not country specific)
  par.names.cs <- mig.parameter.names.cs()#Parameter names (country specific)
  
  default.buffer.size <- 100
  buffer <- buffer.cs <- NULL
  
  buffers.insert <- function(mcmc) {
    counter <<- counter + 1

    for (par in par.names) {
      #Here's how we'll eventually handle parameters that we shouldn't save.
      #        if (is.element(par, mcmc$dontsave)) next
      buffer[[par]][counter,] <<- mcmc[[par]]
    }

    for (par in par.names.cs) {
#      if (is.element(var.names[[par]], mcmc$dontsave)) next
      
      for (country in 1:mcmc$meta$nr.countries){
        result <- mcmc[[par]][country]
        buffer.cs[[par]][[country]][counter,] <<- result
      }
    }
  }
  
  buffers.ini <- function(mcmc, size) {
    #Got rid of the option for custom country lists. (It was buggy anyway.)
    buffer <<- list()
    for (par in par.names) {
      #if (is.element(par, mcmc$dontsave)) next
      buffer[[par]] <<- matrix(NA, ncol=length(mcmc[[par]]), nrow=size)
    }
    
    buffer.cs <<-list()
    for (par in par.names.cs) {
      #if (is.element(var.names[[par]], mcmc$dontsave)) next
      buffer.cs[[par]] <<- list()
      for (country in 1:mcmc$meta$nr.countries){
        v <- mcmc[[par]][country]
        buffer.cs[[par]][[country]] <<- matrix(NA, ncol=length(v), nrow=size)
      }
    }
    counter <<- 0
  }
  
  
  do.flush.buffers <- function(mcmc, append=FALSE, verbose=FALSE) {
    if (verbose)
      cat("Flushing results into disk.\n")
    output.dir <- file.path(mcmc$meta$output.dir, mcmc$output.dir)
    if(!file.exists(output.dir)) 
      dir.create(output.dir)
    open <- if(append) 'a' else 'w'

    for(par in par.names) { # write country-independent parameters
      if (is.null(buffer[[par]])) next
      if (counter == 1) {
        values <- t(buffer[[par]][1:counter,])
      } else {
        values <- buffer[[par]][1:counter,]
      }
      bayesTFR:::write.values.into.file.cindep(par, values, output.dir, mode=open, 
                                    compression.type=mcmc$compression.type)
    }

    for (par in par.names.cs) { # write country-specific parameters
      if (is.null(buffer.cs[[par]])) next
      for (country in 1:mcmc$meta$nr.countries){
        if (counter == 1) {
          values <- t(buffer.cs[[par]][[country]][1:counter,])
        } else {
          values <- buffer.cs[[par]][[country]][1:counter,]
        }
        bayesTFR:::write.values.into.file.cdep(par, values, output.dir, 
                                    get.country.object(country, meta=mcmc$meta, index=TRUE)$code, mode=open, 
                                    compression.type=mcmc$compression.type)
      }
    }
    resmc <- as.list(mcmc)
    class(resmc) <- 'bayesMig.mcmc'
    store.bayesMig.object(resmc, output.dir)
  }
  
  store <- function(mcmc, append=FALSE, flush.buffer=FALSE, verbose=FALSE) {
    buffer.size <- mcmc$meta$buffer.size
    if (is.null(buffer.size)){
      buffer.size <- default.buffer.size
    }
    if (is.null(buffer)){
      buffers.ini(mcmc, buffer.size)      
    }
    buffers.insert(mcmc)
    flushed <- FALSE
    if (flush.buffer || (counter >= buffer.size)) {
      do.flush.buffers(mcmc, append=append, verbose=verbose)
      #The pair of assignments below looks very silly. 
      #The first one is here to avoid an R CMD check note, which reads "no visible binding for '<<-' assignment to 'counter'"
      counter <- 0
      counter <<- 0 #Don't waste energy resetting the whole buffer. Just set the counter back to zero and write over it.
      flushed <- TRUE
    }
    return(flushed)
  }
  
})

store.bayesMig.object <- function(mcmc, output.dir) {
  bayesMig.mcmc <- mcmc
  bayesMig.mcmc$meta <- NULL
  save(bayesMig.mcmc, file=file.path(output.dir, 'bayesMig.mcmc.rda'))
}

store.bayesMig.meta.object <- function(meta, output.dir) {
  bayesMig.mcmc.meta <- meta
  save(bayesMig.mcmc.meta, file=file.path(output.dir, 'bayesMig.mcmc.meta.rda'))
}

#' @title Internal Functions and datasets of bayesMig
#' @description These functions and datasets are not to be used directly by the user.
#' @export
#' @keywords internal
#' @rdname internal
#' 
store.bayesMig.convergence <- function(diag, thin, burnin, output.dir){
  save.file <- file.path(output.dir, paste('bayesTFR.convergence_', thin, '_', burnin, '.rda', sep=''))
  bayesTFR.convergence <- diag
  save(bayesTFR.convergence, file=save.file)
  return(save.file)
}

#' @export
#' @keywords internal
#' @rdname internal
#' 
store.bayesMig.prediction <- function(pred, output.dir=NULL) {
  bayesMig.prediction <- pred
  if (is.null(output.dir)) output.dir <- pred$output.directory
  save(bayesMig.prediction, file=file.path(output.dir, 'prediction.rda'))
}

