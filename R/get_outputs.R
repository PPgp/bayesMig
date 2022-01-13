
##################
#MIGRATION
##################


#' @title Access MCMC results
#'
#' @description This function retrieves results of an MCMC simulation and creates an object of class
#' \code{bayesMig.mcmc.set}. 
#' 
#' @usage get.mig.mcmc(sim.dir=file.path(getwd(), 'bayesMig.output'), chain.ids=NULL,
#' burnin=0, verbose=FALSE)
#' 
#' @param sim.dir Directory where simulation results are stored
#' @param chain.ids Chain identifiers in case only specific chains should be included
#' in the resulting object. By default, all available chains are included.
#' @param burnin Burn-in used for loading traces.
#' @param verbose Logical value. Switches log messages on and off.
#' 

get.mig.mcmc <- function(sim.dir=file.path(getwd(), 'bayesMig.output'), chain.ids=NULL,
                         burnin=0, verbose=FALSE) {
  ############
  # Returns an object of class bayesMig.mcmc.set
  ############
  mcmc.file.path <- file.path(sim.dir, 'bayesMig.mcmc.meta.rda')
  if(!file.exists(mcmc.file.path)) {
    warning('File ', mcmc.file.path, ' does not exist.')
    return(NULL)
  }
  load(file=mcmc.file.path)
  bayesMig.mcmc.meta$output.dir <- normalizePath(sim.dir)
  if (is.null(chain.ids)) {
    mc.dirs.short <- list.files(sim.dir, pattern='^mc[0-9]+', full.names=FALSE)
    chain.ids <- as.integer(substring(mc.dirs.short, 3))
  } else {
    mc.dirs.short <- paste('mc', chain.ids, sep='')
  }
  ord.idx <- order(chain.ids)
  mc.dirs.short <- mc.dirs.short[ord.idx]
  chain.ids <- chain.ids[ord.idx]
  mcmc.chains <- list()
  counter<-1
  for (imc.d in chain.ids) {
    if (verbose)
      cat('Loading chain', imc.d, 'from disk. ')
    bayesMig.mcmc <- local({
      load(file=file.path(sim.dir, mc.dirs.short[counter], 'bayesMig.mcmc.rda'))
      bayesMig.mcmc})
    mc <- c(bayesMig.mcmc, list(meta=bayesMig.mcmc.meta))
    class(mc) <- class(bayesMig.mcmc)
    #load full mcmc traces
    #th.burnin <- get.thinned.burnin(mc, burnin)
    #mc$traces <- load.mig.parameter.traces.all(mc, burnin=th.burnin)
    th.burnin <- 0 # low memory (traces will be loaded as they are needed)
    mc$traces <- 0
    mc$traces.burnin <- th.burnin
    
    mc$output.dir <- mc.dirs.short[counter]
    if (verbose)
      cat('(mcmc.list[[', counter, ']]).\n')
    mcmc.chains[[counter]] <- mc
    counter <- counter+1
  }
  names(mcmc.chains) <- chain.ids
  return(structure(list(meta=bayesMig.mcmc.meta, 
                        mcmc.list=mcmc.chains), class='bayesMig.mcmc.set'))
}

get.mig.prediction <- function(mcmc=NULL, sim.dir=NULL, mcmc.dir=NULL) {
  ############
  # Returns an object of class bayesMig.prediction
  # Set mcmc.dir to NA, if the prediction object should not have a pointer 
  # to the corresponding mcmc traces
  ############
  if (!is.null(mcmc)) 
    sim.dir <- if(is.character(mcmc)) mcmc else mcmc$meta$output.dir
  if (is.null(sim.dir)) stop('Either mcmc or directory must be given.')
  output.dir <- file.path(sim.dir, 'predictions')
  pred.file <- file.path(output.dir, 'prediction.rda')
  if(!file.exists(pred.file)) {
    warning('File ', pred.file, ' does not exist.')
    return(NULL)
  }
  load(file=pred.file)
  bayesMig.prediction$output.directory <- output.dir

  pred <- bayesMig.prediction
  # re-route mcmcs if necessary
  if(!is.null(mcmc.dir) || !has.mig.mcmc(pred$mcmc.set$meta$output.dir)) {
    if((!is.null(mcmc.dir) && !is.na(mcmc.dir)) || is.null(mcmc.dir)) {
      new.path <- file.path(sim.dir, basename(pred$mcmc.set$meta$output.dir))
      if (has.mig.mcmc(new.path)) pred$mcmc.set <- get.mig.mcmc(new.path)
      else {
        est.dir <- if(is.null(mcmc.dir)) sim.dir else mcmc.dir
        pred$mcmc.set <- get.mig.mcmc(est.dir)
      }
    }
  }
  return(pred)
}

has.mig.mcmc <- function(sim.dir) {
  return(file.exists(file.path(sim.dir, 'bayesMig.mcmc.meta.rda')))
}


get.thinned.burnin <- function(mcmc, burnin) {
  if (burnin==0) return(0)
  if (mcmc$thin == 1) return(burnin)
  return(1 + if(burnin >= mcmc$thin) length(seq(mcmc$thin, burnin, by=mcmc$thin)) else 0)
}


load.mig.parameter.traces.all <- function(mcmc, par.names=mig.parameter.names(), 
                                          par.names.cs=mig.parameter.names.cs(),
                                          burnin=0, thinning.index=NULL) {
  result <- load.mig.parameter.traces(mcmc, par.names, burnin=burnin, thinning.index=thinning.index)
  if(!is.null(par.names.cs))
    for (country in get.countries.index(mcmc$meta)) {
      result <- cbind(result, 
                      load.mig.parameter.traces.cs(mcmc, 
                                                   get.country.object(country, 
                                                                      mcmc$meta, index=TRUE)$code, 
                                                   par.names.cs, burnin=burnin,
                                                   thinning.index=thinning.index))
    }
  return (result)
}

"get.countries.index" <- function(meta, ...) UseMethod("get.countries.index")

get.countries.index.bayesMig.mcmc.meta  <- function(meta, ...) 
  return (meta$countryIndices)

mig.parameter.names <- function() {
  # Return all country-independent parameter names. 
  return(c("a","b","mu_global","sigma2_mu"))
}

mig.parameter.names.cs <- function(country.code=NULL) {
  #Return all country-specific parameter names. 
  #If country is not NULL, it must be a country code.
  #It is attached to the parameter name.
  par.names <- c("mu_c","phi_c","sigma2_c")
  if (is.null(country.code))
    return(par.names)
  return(paste(par.names, '_country', country.code, sep=''))
}

load.mig.parameter.traces <- function(mcmc, par.names=NULL, burnin=0, thinning.index=NULL) {
  if(missing(par.names)) par.names <- mig.parameter.names()
  return(bdem.parameter.traces(mcmc, par.names, burnin=burnin, thinning.index=thinning.index))
}

load.mig.parameter.traces.cs <- function(mcmc, country, par.names=NULL, burnin=0, 
                                         thinning.index=NULL) {
  #The par.names input should look something like c("mu_c","phi_c","sigma2_c")
  if(missing(par.names)) par.names <- mig.parameter.names.cs()
  return(bdem.parameter.traces(mcmc, par.names, paste0("_country",country),
                               burnin=burnin, thinning.index=thinning.index))
}

bdem.parameter.traces.bayesMig.mcmc <- function(mcmc, par.names, ...) {
  # Load traces from the disk
  all.standard.names <- c(mig.parameter.names(), mig.parameter.names.cs())
  return(bayesTFR:::.do.get.traces(mcmc, par.names=par.names, ..., all.standard.names = all.standard.names))
}


coda.mcmc.bayesMig.mcmc <- function(mcmc, country=NULL, par.names=NULL, 
                                    par.names.cs=NULL, burnin=0, thin=1, ...
                                    ) {
  if(missing(par.names)) par.names <- mig.parameter.names()
  if(missing(par.names.cs)) par.names.cs <- mig.parameter.names.cs()
  return(bayesTFR:::coda.mcmc.bayesTFR.mcmc(mcmc, country = country, par.names = par.names, 
                                            par.names.cs = par.names.cs, ...))
}

mig.coda.list.mcmc <- function(mcmc.list = NULL, country = NULL, chain.ids = NULL,
                              sim.dir = file.path(getwd(), 'bayesMig.output'), 
                              par.names = NULL, par.names.cs = NULL, 
                              low.memory = FALSE, ...) {
  # return a list of mcmc objects that can be analyzed using the coda package
  if (is.null(mcmc.list)) {
    mcmc.list <- get.mig.mcmc(sim.dir, chain.ids=chain.ids)$mcmc.list
  } else {
    mcmc.list <- get.mcmc.list(mcmc.list)
    if (!is.null(chain.ids)) {
      mcmc.list <- mcmc.list[chain.ids]
    }
  }
  if(missing(par.names)) par.names <- mig.parameter.names()
  if(missing(par.names.cs)) par.names.cs <- mig.parameter.names.cs()
  result <- list()
  i <- 1
  for(mcmc in mcmc.list) {
    result[[i]] <- coda.mcmc(mcmc, country = country, par.names = par.names, 
                             par.names.cs = par.names.cs, ...)
    i <- i+1
  }
  return(mcmc.list(result))
}

get.full.par.names.cs <- function(par.names, full.par.names, country=NULL, index=FALSE) {
  # Return full name of par.names that are included in full.par.names
  # which are suppose to be all country-specific parameters.
  # E.g. for 'phi_c', it would return 'phi_c_country1', ..., 'phi_c_country233'
  # If index is TRUE, return the index of the matches.
  result <- c()	
  for (name in par.names) {
    pattern <- paste('^',name,'_.*',sep='')
    if (!is.null(country)) {
      pattern <- paste(pattern, 'country', country, '$', sep='')
    } else {
      pattern <- paste(pattern, '[[:digit:]]$', sep='')
    }
    result <- c(result, grep(pattern, full.par.names, value=!index))
  }
  return(result)
}

get.full.par.names <- function(par.names, full.par.names, index=FALSE) {
  # Return full name of par.names that are included in full.par.names
  # which are suppose to be all country-independent parameters.
  # E.g. for 'a', it would return just 'a' (and not 'alpha')
  # If index is TRUE, return the index of the matches.
  result <- c()	
  for (name in par.names) {
    pattern <- paste('^', name, '$', sep='') # matches exactly
    result <- c(result, grep(pattern, full.par.names, value=!index))
  }
  return(result)
}

get.burned.mig.traces <- function(mcmc, par.names, burnin=0, thinning.index=NULL) {
  # get traces that are already loaded in the mcmc object
  traces <- matrix(mcmc$traces[, par.names],ncol=length(par.names), dimnames=list(NULL, par.names))
  discard <- burnin - mcmc$traces.burnin
  if (discard > 0)
    traces <- traces[-seq(1, discard),, drop = FALSE]
  if(!is.null(thinning.index))
    traces <- traces[thinning.index, , drop = FALSE]
  return(traces)
}

#' Conversion to coda-formatted objects
#' 
#' @usage coda.list.mcmc(mcmc=NULL, country=NULL, chain.ids=NULL, 
#' sim.dir=file.path(getwd(), 'bayesMig.output'),
#' par.names=mig.parameter.names(), 
#' par.names.cs=mig.parameter.names.cs(), 
#' rm.const.pars=FALSE, burnin=0, ...)
#' 
#' @param mcmc A list of objects of class \code{bayesMig.mcmc}. If \code{NULL}, the MCMCs are
#' loaded from \code{sim.dir}. Either \code{mcmc} or \code{sim.dir} must be given.
#' @param country Country name or code. Used in connection with the \code{par.names.cs} argument
#' (see below)
#' @param chain.ids Vector of chain identifiers. By default, all chains available in the \code{mcmc.list}
#' object are included
#' @param sim.dir Directory with the MCMC simulation results. Only used if \code{mcmc} is \code{NULL}
#' @param par.names Names of country-independent parameters to be included. Default names are
#' those returned by the \code{mig.parameter.names} function, which includes all country-independent
#' parameters in the BHM.
#' @param par.names.cs Names of country-specific parameters to be included. The argument \code{country}
#' is used to filter out traces that correspond to a specific country. If \code{country} is not given, 
#' traces of each parameter are given for all countries. Default names are those returned by 
#' \code{mig.parameter.names.cs()}, which includes all country-specific parameters in the BHM.
#' @param rm.const.pars Logical indicating if parameters with constant values should be removed.
#' @param burnin Number of iterations that should be removed from the beginning of each chain.
#' @param ... Other variables passed to called functions

coda.list.mcmc <- function(mcmc=NULL, country=NULL, chain.ids=NULL,
                           sim.dir=file.path(getwd(), 'bayesMig.output'), 
                           par.names=mig.parameter.names(), 
                           par.names.cs=mig.parameter.names.cs(), 
                           rm.const.pars=FALSE,
                           burnin=0, 
                           ...) {
  # return a list of mcmc objects that can be analyzed using the coda package
  mcmc.list <- mcmc
  if (is.null(mcmc.list)) {
    mcmc.list <- get.mig.mcmc(sim.dir, chain.ids=chain.ids, burnin=burnin)$mcmc.list
  } else {
    mcmc.list <- get.mcmc.list(mcmc.list)
    if (!is.null(chain.ids)) {
      mcmc.list <- mcmc.list[chain.ids]
    }
  }
  result <- list()
  i <- 1
  for(mc in mcmc.list) {
    result[[i]] <- coda.mcmc(mc, country=country, par.names=par.names, par.names.cs=par.names.cs,
                             burnin=burnin, ...)
    
    #########EDIT
    if (rm.const.pars) {
      # remove parameters that are constant for all iterations
      if (is.null(dim(result[[i]]))) { # one parameter in the chain
        if(all(result[[i]]==result[[i]][1]))  { # no parameters are kept
          result[[i]] <- NULL
          i <- i - 1
        }
      } else { # multiple parameters in the chain
        if (dim(result[[i]])[2]==1) {
          colindex <- !all(result[[i]][,1] == result[[i]][1,1])
        } else {
          colindex <- !apply(t(apply(result[[i]], 1, '==', result[[i]][1,])), 2, all)
        }
        if (sum(colindex) == 0) { # no parameters are kept
          result[[i]] <- NULL
          i <- i - 1
        } else {
          result[[i]] <- result[[i]][,colindex]
        }
      }
    }
    i <- i+1
  }
  return(mcmc.list(result))
}

get.mcmc.list.bayesMig.mcmc.set <- function(mcmc.list, ...) return(mcmc.list$mcmc.list)
get.mcmc.list.bayesMig.mcmc <- function(mcmc.list, ...) return(list(mcmc.list))
get.mcmc.list.bayesMig.prediction <- function(mcmc.list, ...) return(mcmc.list$mcmc.set$mcmc.list)


get.total.iterations <- function(mcmc.list, burnin=0) {
  # Return total number of iterations sum-up over chains after discarding burnin in each chain
  get.iter <- function(x) return(x$finished.iter - burnin)
  return(sum(sapply(mcmc.list, get.iter)))
}

get.thinned.mig.mcmc <- function(mcmc.set, thin=1, burnin=0) {
  dir.name <- file.path(mcmc.set$meta$output.dir, paste('thinned_mcmc', thin, burnin, sep='_'))
  if(file.exists(dir.name)) return(get.mig.mcmc(dir.name))
  return(NULL)
}

get.stored.mcmc.length <- function(mcmc.list, burnin=0) {
  # Return total number of iterations sum-up over chains after discarding burnin in each chain,
  # taking into account the original value of thin.
  # It should correspond to the total length of all chains stored on disk, minus burnin.
  get.iter <- function(x) return(x$length - get.thinned.burnin(x, burnin))
  return(sum(sapply(mcmc.list, get.iter)))
}

has.mig.prediction <- function(mcmc=NULL, sim.dir=NULL) {
  if (!is.null(mcmc)) sim.dir <- if(is.character(mcmc)) mcmc else mcmc$meta$output.dir
  if (is.null(sim.dir)) stop('Either mcmc or directory must be given.')
  if(file.exists(file.path(sim.dir, 'predictions', 'prediction.rda'))) return(TRUE)
  return(FALSE)
}



get.mig.parameter.traces <- function(mcmc.list, par.names=NULL, 
                                     burnin=0, thinning.index=NULL, thin=NULL) {
  # get parameter traces either from disk or from memory, if they were already loaded
  mcmc.list <- get.mcmc.list(mcmc.list)
  if(missing(par.names)) par.names <- mig.parameter.names()
  return(bayesTFR:::do.get.tfr.parameter.traces(is.cs=FALSE, mcmc.list=mcmc.list, par.names=par.names, 
                                     burnin=burnin, thinning.index=thinning.index, thin=thin))
}

get.mig.parameter.traces.cs <- function(mcmc.list, country.obj, par.names=NULL, 
                                        burnin=0, thinning.index=NULL, thin=NULL) {
  # country.obj is result of get.country.object()
  # get traces for country-specific parameters either from disk or from memory, if they were already loaded
  mcmc.list <- get.mcmc.list(mcmc.list)
  if(missing(par.names)) par.names <- mig.parameter.names.cs()
  return(bayesTFR:::do.get.tfr.parameter.traces(is.cs=TRUE, mcmc.list=mcmc.list, par.names=par.names, country.obj=country.obj,
                                     burnin=burnin, thinning.index=thinning.index, thin=thin))
}


create.thinned.mig.mcmc <- function(mcmc.set, thin=1, burnin=0, output.dir=NULL, verbose=TRUE) {
  #Return a thinned mcmc.set object with burnin removed and all chanins collapsed into one
  mcthin <- max(sapply(mcmc.set$mcmc.list, function(x) x$thin))
  thin <- max(c(thin, mcthin))
  meta <- mcmc.set$meta
  total.iter <- get.stored.mcmc.length(mcmc.set$mcmc.list, burnin=burnin)
  meta$is.thinned <- TRUE
  meta$parent.iter <- get.total.iterations(mcmc.set$mcmc.list, burnin)
  meta$parent.meta <- mcmc.set$meta
  meta$nr.chains <- 1
  
  if(verbose) cat('\nStoring thinned mcmc:')
  # store the meta object
  meta$output.dir <- file.path(
    if(is.null(output.dir)) meta$output.dir else output.dir, 
    paste('thinned_mcmc', thin, burnin, sep='_'))
  if(!file.exists(meta$output.dir)) 
    dir.create(meta$output.dir, recursive=TRUE)
  store.bayesMig.meta.object(meta, meta$output.dir)
  
  thin.index <- if(thin > mcthin) unique(round(seq(1, total.iter, by=thin/mcthin))) else 1:total.iter
  nr.points <- length(thin.index)
  
  #create one collapsed mcmc
  thinned.mcmc <- mcmc.set$mcmc.list[[1]]
  thinned.mcmc$meta <- meta
  thinned.mcmc$thin <- 1
  thinned.mcmc$id <- 1
  thinned.mcmc$traces <- 0
  thinned.mcmc$length <- nr.points
  thinned.mcmc$finished.iter <- nr.points
  thinned.mcmc$output.dir <- 'mc1'	
  outdir.thin.mcmc <- file.path(meta$output.dir, 'mc1')
  if(!file.exists(outdir.thin.mcmc)) dir.create(outdir.thin.mcmc)
  
  store.bayesMig.object(thinned.mcmc, outdir.thin.mcmc)
  
  if(verbose) cat('\nStoring country-independent parameters ...')
  for (par in mig.parameter.names()) {
    values <- get.mig.parameter.traces(mcmc.set$mcmc.list, par, burnin,
                                       thinning.index=thin.index)
    #values = as.vector(t(values))#JA: This seems to be necessary to switch from multiple rows to a single vector.
    bayesTFR:::write.values.into.file.cindep(par, values, outdir.thin.mcmc, compression.type='None')
  }
  if(verbose) cat('done.\nStoring country-specific parameters ...')
  par.names.cs <- mig.parameter.names.cs()
  for (country in mcmc.set$meta$countryIndices){
    country.obj <- get.country.object(country, mcmc.set$meta, index=TRUE)
    for (par in par.names.cs) {
      values <- get.mig.parameter.traces.cs(mcmc.set$mcmc.list, country.obj, par, 
                                            burnin=burnin, thinning.index=thin.index)
      #values = as.vector(t(values))#JA: This seems to be necessary to switch from multiple rows to a single vector.
      bayesTFR:::write.values.into.file.cdep(par, values, outdir.thin.mcmc, country.code=country.obj$code,
                                  compression.type='None')
    }
  }
  if(verbose) cat('done.\n')
  #JA: We'll assume all countries are used for estimation.
  # if (mcmc.set$meta$nr_countries > mcmc.set$meta$nr_countries_estimation) {
  #   .update.thinned.extras(mcmc.set, (mcmc.set$meta$nr_countries_estimation+1):mcmc.set$meta$nr_countries,
  #                          burnin=burnin, nr.points=nr.points, dir=outdir.thin.mcmc, verbose=verbose)
  # }
  invisible(structure(list(meta=meta, mcmc.list=list(thinned.mcmc)), class='bayesMig.mcmc.set'))
}

# .update.thinned.extras <- function (mcmc.set, country.index, burnin, nr.points, dir, verbose=TRUE) {
#   if(verbose) cat('done.\nStoring country-specific parameters for extra countries ...')
#   # thin mcmc for extra countries (they can have different length than the other countries)
#   par.names.cs <- mig.parameter.names.cs()
#   for (country in country.index){
#     country.obj <- get.country.object(country, mcmc.set$meta, index=TRUE)
#     for (par in par.names.cs) {
#       values <- get.mig.parameter.traces.cs(mcmc.set$mcmc.list, country.obj, par, 
#                                             burnin=burnin)
#       selected.simu <- get.thinning.index(nr.points, dim(values)[1])
#       if (length(selected.simu$index) < nr.points)
#         selected.simu$index <- sample(selected.simu$index, nr.points, replace=TRUE)
#       values <- values[selected.simu$index,]
#       write.values.into.file.cdep(par, values, dir, country.code=country.obj$code, 
#                                   compression.type='None')
#     }
#   }
#   if(verbose) cat('done.\n')
# }

# get.thinning.index <- function(nr.points, all.points) {
#   if (!is.null(nr.points)) {
#     nr.points <- ifelse(nr.points >= all.points, all.points, nr.points)
#   } else {
#     nr.points <- all.points
#   }
#   if (nr.points > 0) {
#     step <- all.points/nr.points
#     idx <- floor(seq(floor(step), all.points, by=step))
#   } else idx<-NULL
#   return(list(nr.points=nr.points, index=idx))
# }


summary.bayesMig.mcmc <- function(object, country = NULL, 
                                   par.names = NULL, par.names.cs = NULL, thin = 1, burnin = 0, ...) {
  res <- list()
  class(res) <- "summary.bayesTFR.mcmc"
  if(missing(par.names)) par.names <- mig.parameter.names()
  if(missing(par.names.cs)) par.names.cs <- mig.parameter.names.cs()
  if (!is.null(country)) {
    country.obj <- get.country.object(country, object$meta)
    if(is.null(country.obj$name)) stop("Country ", country, " not found.")
    res$country.name <- country.obj$name
    country <- country.obj$code
  } 
  res$results <- summary(coda.mcmc(object, country=country, par.names=par.names,
                                   par.names.cs=par.names.cs, thin=thin, burnin=burnin), ...)
  return(res)
}

summary.bayesMig.mcmc.set <- function(object, country=NULL, chain.id=NULL, 
                                       par.names = NULL, 
                                       par.names.cs = NULL, 
                                       meta.only=FALSE, thin=1, burnin=0, ...) {
  if(missing(par.names)) par.names <- mig.parameter.names()
  if(missing(par.names.cs) && !is.null(country)) par.names.cs <- mig.parameter.names.cs()
  
  res <- list(meta = summary(object$meta))
  class(res) <- "summary.bayesMig.mcmc.set"
  if(meta.only) {
    res$chain.info <- bayesTFR:::chain.info(object)
    return(res)
  }
  if (!is.null(chain.id)) {
    res$mcmc <- summary(object$mcmc.list[[chain.id]], country = country, 
                        par.names = par.names, par.names.cs = par.names.cs, 
                        thin = thin, burnin = burnin, ...)
    return(res)
  }
  if (!is.null(country)) {
    country.obj <- get.country.object(country, object$meta)
    if(is.null(country.obj$name)) stop("Country ", country, " not found.")
    res$country.name <- country.obj$name
    country <- country.obj$code
  }
  res$results <- summary(coda.list.mcmc(object, country = country, par.names = par.names,
                                        par.names.cs = par.names.cs, thin = thin, 
                                        burnin = burnin), ...)
  return(res)
}

summary.bayesMig.mcmc.meta <- function(object, ...) {
  res <- list(est.period = paste(object$start.year, object$present.year, sep = '-'),
              nr.countries = object$nr.countries,
              data.source = if(object$user.data) "user-defined" else "WPP ",
              wpp.year = if(object$user.data) NULL else object$wpp.year
  )
  class(res) <- "summary.bayesMig.mcmc.meta"
  return(res)
}

print.summary.bayesMig.mcmc.meta <- function(x, ...) {
  cat('\nNumber of countries:', x$nr.countries)
  cat('\nData source:', x$data.source, x$wpp.year)
  cat('\nInput data: migration for period', x$est.period)
  cat('\n')
}

print.summary.bayesMig.mcmc.set <- function(x, ...) {
  print(x$meta)
  if(!is.null(x$chain.info)) bayesTFR:::print.summary.chain.info(x$chain.info)
  if(!is.null(x$mcmc)) print(x$mcmc)
  if(!is.null(x$country.name)){
    cat('\nCountry:', x$country.name, '\n')
    if (is.null(x$results))
      cat('\tnot used for estimation.\n')
  }
  if(!is.null(x$results)) print(x$results)
}

print.bayesMig.mcmc <- function(x, ...) {
  print(summary(x, ...))
}

print.summary.bayesMig.mcmc <- function(x, ...) {
  if(!is.null(x$country.name)){
    cat('\nCountry:', x$country.name, '\n')
    if (is.null(x$results))
      cat('\tnot used for estimation.\n')
  }
  if(!is.null(x$results))
    print(x$results)
}

print.bayesMig.mcmc.set <- function(x, ...) {
  print(summary(x, ...))
}

print.bayesMig.mcmc.meta <- function(x, ...) {
  print(summary(x, ...))
}

print.bayesMig.prediction <- function(x, ...) {
  print(summary(x, ...))
}

summary.bayesMig.prediction <- function(object, country = NULL, compact = TRUE, ...) {
  res <- bayesTFR:::get.prediction.summary.data(object, 
                                                unchanged.pars=c('burnin', 'nr.traj'), 
                                                country=country, compact=compact)
  class(res) <- 'summary.bayesMig.prediction'
  return(bayesTFR:::.update.summary.data.by.shift(res, object, country))
}

print.summary.bayesMig.prediction <- function(x, digits = 3, ...) {
  cat('\nProjections:', length(x$projection.years), '(', x$projection.years[1], '-', 
      x$projection.years[length(x$projection.years)], ')')
  cat('\nTrajectories:', x$nr.traj)
  cat('\nBurnin:', x$burnin)
  
  if(!is.null(x$country.name)) {
    cat('\nCountry:', x$country.name, '\n')
    cat('\nProjected Migration Rate:')
    cat('\n')
    print(x$projections, digits=digits, ...)
  }
}


get.nr.countries.bayesMig.mcmc.meta <- function(meta, ...) 
  return(meta$nr.countries)

get.nrest.countries.bayesMig.mcmc.meta <- function(meta, ...) 
  return(meta$nr.countries)

country.names.bayesMig.mcmc.meta <- function(meta) {
  #JA: This assumes we produce projections for all listed countries.
  return(meta$fullCountryNameVec)
}

get.countries.table.bayesMig.mcmc.set <- function(object, ...) 
  return(bayesTFR:::get.countries.table.bayesTFR.mcmc.set(object,...))
get.countries.table.bayesMig.prediction <- function(object, ...) 
  return(bayesTFR:::get.countries.table.bayesTFR.prediction(object,...))

get.data.matrix.bayesMig.mcmc.meta <- function(meta, ...) return (t(meta$mig.rates))

get.mcmc.meta.bayesMig.mcmc.set <- function(meta, ...) return(meta$meta)
get.mcmc.meta.bayesMig.mcmc.meta <- function(meta, ...) return(meta)
get.mcmc.meta.bayesMig.mcmc <- function(meta, ...) return(meta$meta)
get.mcmc.meta.list <- function(meta, ...) return(meta[[1]]$meta)
get.mcmc.meta.bayesMig.prediction <- function(meta, ...) return(meta$mcmc.set$meta)



#########
#Test code!
#########

#load("./bayesMig.output/bayesMig.mcmc.meta.rda")
#load("./bayesMig.output/mc1/bayesMig.mcmc.rda")
#mc <- c(bayesMig.mcmc, list(meta=bayesMig.mcmc.meta))
#class(mc) <- class(bayesMig.mcmc)
#result=load.mig.parameter.traces.all(mc)

#mcmc=get.mig.mcmc()
#result=load.mig.parameter.traces.all(mcmc[[1]])

#o=get.mig.mcmc(verbose=TRUE)

#o=coda.list.mcmc()#Grabs from default directory.

# do.get.mig.parameter.traces(is.cs=TRUE, mcmc.list=mcmc.set$mcmc.list, par.names=c("mu_c","phi_c"), 
#                             country.obj=get.country.object(1,meta=mcmc$meta,index=TRUE), 
#                                         burnin=0, thinning.index=NULL, thin=NULL)
# 
# load.mig.parameter.traces.cs(mcmc=load.mcmc.set$mcmc.list[[1]],
#                              country=4,
#                              par.names=c("mu_c","phi_c"))
