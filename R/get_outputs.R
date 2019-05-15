
##################
#MIGRATION
##################

get.country.object <- function(country, meta=NULL, country.table=NULL, index=FALSE) {
  #TABLE OPTION NOT IMPLEMENTED. At the moment, everything has to be handled with a meta object.
  
  # If meta object is not given, country.table containing columns 'code' and 'name' must be given. 
  # If 'country' is numeric, 'index' determines if 'country' is an index (TRUE) or code (FALSE)
  if (!is.null(meta)) {
    codes <- meta$regions$country_code
    names <- meta$regions$country_name
  } else {
    codes <- country.table[,'code']
    names <- country.table[,'name']	
  }
  l <- length(codes)
  found <- TRUE
  if (is.numeric(country)) {#If "country" is a number
    if (index) {
      country.idx <- country
      country.code <- codes[country.idx]
    } else { 
      country.code <- country
      country.idx <- (1:l)[codes==country.code]
    }
    country.name <- as.character(names[country.idx])
    if (length(country.name) == 0) found <- FALSE
  } else {#If "country" is non-numeric, assume it's a name.
    country.name <- country
    country.idx <- (1:l)[names==country.name]
    country.code <- codes[country.idx]
    if (length(country.idx) == 0) found <- FALSE
  }
  if (!found) 
    country.name <- country.idx <- country.code <- NULL
  return(list(name=country.name, index=country.idx, code=country.code))
}

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
    th.burnin <- get.thinned.burnin(mc, burnin)
    mc$traces <- load.mig.parameter.traces.all(mc, burnin=th.burnin)
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

load.mig.parameter.traces <- function(mcmc, par.names=mig.parameter.names(), burnin=0, thinning.index=NULL) 
  return(bdem.parameter.traces(mcmc, par.names, burnin=burnin, thinning.index=thinning.index))

load.mig.parameter.traces.cs <- function(mcmc, country, par.names=NULL, burnin=0, 
                                         thinning.index=NULL) {
  #The par.names input should look something like c("mu_c","phi_c","sigma2_c")
  if(is.null(par.names)){
    stop("Need country-specific parameter names for load.mig.parameter.traces.cs")
  }
  return(bdem.parameter.traces(mcmc, par.names=paste(par.names,"_country",country,sep=""),
                               burnin=burnin, 
                               thinning.index=thinning.index))
}

"bdem.parameter.traces" <- function(mcmc, ...) UseMethod("bdem.parameter.traces")

bdem.parameter.traces.bayesMig.mcmc <- function(mcmc, par.names, ...) {
  return(.do.get.traces(mcmc, par.names=par.names, ...))
}

.do.get.traces <- function(mcmc, par.names, burnin=0, 
                           thinning.index=NULL) {  
  values <- c()
  for(name in par.names) {
    file.name <- paste(name, '.txt', sep='')
    vals <- as.matrix(read.table(file.path(mcmc$meta$output.dir, mcmc$output.dir, file.name)))
    if (burnin > 0) {
      if (burnin > dim(vals)[1]) stop('Burnin is larger than the data size.')
      vals <- vals[-seq(1, burnin),,drop=FALSE]
    }
    if(!is.null(thinning.index))
      vals <- vals[thinning.index,,drop=FALSE]
    values <- cbind(values, vals)
  } # end of loading loop
  colnames(values) <- par.names

  return(values)
}

"coda.mcmc" <- function(mcmc, ...) UseMethod("coda.mcmc")

coda.mcmc.bayesMig.mcmc <- function(mcmc, country=NULL, par.names=NULL, 
                                    par.names.cs=NULL, burnin=0, thin=1, ...
) {
  # Return a coda object for this mcmc and parameter names
  index <- NULL
  btobject <- burn.and.thin(mcmc, burnin=burnin, thin=thin)
  thin <- btobject$thin
  th.burnin <- btobject$burnin
  if(!is.null(btobject$index)) index <- btobject$index
  if(missing(par.names)) {
    par.names <- mig.parameter.names()    
  }
  if(missing(par.names.cs)){
    par.names.cs <- mig.parameter.names.cs()
  }

  if (!is.null(country)) { # for specific country
    if (burnin < mcmc$traces.burnin) {
      values <- load.mig.parameter.traces(mcmc, par.names, burnin=th.burnin, thinning.index=index)
      values <- cbind(values, load.mig.parameter.traces.cs(mcmc, 
                                                           get.country.object(country, mcmc$meta)$code, par.names.cs, 
                                                           burnin=th.burnin, thinning.index=index))
    } else {
      values <- get.burned.mig.traces(mcmc, 
                                      c(get.full.par.names(par.names, colnames(mcmc$traces)), 
                                        get.full.par.names.cs(par.names.cs, colnames(mcmc$traces), 
                                                              get.country.object(country, mcmc$meta)$code)), th.burnin,
                                      thinning.index=index)
    }
    
  } else { #no country specified
    #Again, we'll assume traces have been loaded, but we may need to get the right burnin.
    values <- get.burned.mig.traces(mcmc, c(get.full.par.names(par.names, colnames(mcmc$traces)), 
                                            get.full.par.names.cs(par.names.cs, colnames(mcmc$traces))), th.burnin,
                                    thinning.index=index)
  }
  # filter out unnecessary parameters here if needed.

  return(mcmc(values, end=mcmc$finished.iter, thin=thin, ...))
}

burn.and.thin <- function(mcmc, burnin=0, thin=1) {
  # Return thin and burnin that is consolidated with the original thin used for storing mcmc.
  # If there is need for more thinning, it returns the corresponding thinning index
  th.burnin <- get.thinned.burnin(mcmc,burnin)
  index <- NULL
  if (thin > mcmc$thin) {
    index <- unique(round(seq(1,mcmc$length-th.burnin, by=thin/mcmc$thin)))
  }
  thin <- max(mcmc$thin, thin) # thin cannot be smaller than the original thin
  return(list(thin=thin, burnin=th.burnin, index=index))
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

"get.mcmc.list" <- function(mcmc.list, ...) UseMethod("get.mcmc.list")

get.mcmc.list.bayesMig.mcmc.set <- function(mcmc.list, ...) return(mcmc.list$mcmc.list)
get.mcmc.list.bayesMig.mcmc <- function(mcmc.list, ...) return(list(mcmc.list))
get.mcmc.list.list <- function(mcmc.list, ...) return(mcmc.list)

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



get.mig.parameter.traces <- function(mcmc.list, par.names=mig.parameter.names(), 
                                     burnin=0, thinning.index=NULL, thin=NULL) {
  # get parameter traces either from disk or from memory, if they were already loaded
  return(do.get.mig.parameter.traces(is.cs=FALSE, mcmc.list=mcmc.list, par.names=par.names, 
                                     burnin=burnin, thinning.index=thinning.index, thin=thin))
}

get.mig.parameter.traces.cs <- function(mcmc.list, country.obj, par.names=mig.parameter.names.cs(), 
                                        burnin=0, thinning.index=NULL, thin=NULL) {
  # country.obj is result of get.country.object()
  # get traces for country-specific parameters either from disk or from memory, if they were already loaded
  return(do.get.mig.parameter.traces(is.cs=TRUE, mcmc.list=mcmc.list, par.names=par.names, country.obj=country.obj,
                                     burnin=burnin, thinning.index=thinning.index, thin=thin))
}

do.get.mig.parameter.traces <- function(is.cs, mcmc.list, par.names, country.obj=NULL, 
                                        burnin=0, thinning.index=NULL, thin=NULL) {
  # get parameter traces either from disk or from memory (if they were already loaded)
  # par.names are either country-independent (if is.cs is FALSE), or country-specific (if is.cs is TRUE)
  values <- c()
  if (is.null(thinning.index) && !is.null(thin) && thin > mcmc.list[[1]]$thin) {
    total.iter <- get.stored.mcmc.length(mcmc.list, burnin)
    thinning.index <- unique(round(seq(1, total.iter, by=thin/mcmc.list[[1]]$thin)))
  }
  at.iter <- 1
  for(mcmc in mcmc.list) {
    this.thinning.index <- NULL
    th.burnin <- get.thinned.burnin(mcmc, burnin)
    if(!is.null(thinning.index)) {
      this.thinning.index <- thinning.index[(thinning.index >= at.iter) & 
                                              (thinning.index < at.iter+mcmc$length-th.burnin)] - at.iter+1
      if (length(this.thinning.index) == 0) {
        at.iter <- at.iter+mcmc$length-th.burnin
        next
      }
    }
    if (no.traces.loaded(mcmc) || th.burnin < mcmc$traces.burnin) {
      traces <- if(is.cs) load.mig.parameter.traces.cs(mcmc, country.obj$code, par.names, burnin=th.burnin, 
                                                       thinning.index=this.thinning.index)
      else bdem.parameter.traces(mcmc, par.names, burnin=th.burnin, thinning.index=this.thinning.index)
    } else {
      traces <- if(is.cs) get.burned.mig.traces(mcmc, get.full.par.names.cs(par.names, colnames(mcmc$traces), 
                                                                            country=country.obj$code),#Changed index to code here!
                                                th.burnin, thinning.index=this.thinning.index)
      else get.burned.mig.traces(mcmc, par.names, th.burnin, thinning.index=this.thinning.index)
    }
    values <- rbind(values, traces)
    at.iter <- at.iter+mcmc$length-th.burnin
  }
  return(values)
}

no.traces.loaded <- function(mcmc) return((length(mcmc$traces) == 1) && mcmc$traces == 0)

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
    values = as.vector(t(values))#JA: This seems to be necessary to switch from multiple rows to a single vector.
    write.values.into.file.cindep(par, values, outdir.thin.mcmc, compression.type='None')
  }
  if(verbose) cat('done.\nStoring country-specific parameters ...')
  par.names.cs <- mig.parameter.names.cs()
  for (country in mcmc.set$meta$countryIndices){
    country.obj <- get.country.object(country, mcmc.set$meta, index=TRUE)
    for (par in par.names.cs) {
      values <- get.mig.parameter.traces.cs(mcmc.set$mcmc.list, country.obj, par, 
                                            burnin=burnin, thinning.index=thin.index)
      values = as.vector(t(values))#JA: This seems to be necessary to switch from multiple rows to a single vector.
      write.values.into.file.cdep(par, values, outdir.thin.mcmc, country.code=country.obj$code,
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

.update.thinned.extras <- function (mcmc.set, country.index, burnin, nr.points, dir, verbose=TRUE) {
  if(verbose) cat('done.\nStoring country-specific parameters for extra countries ...')
  # thin mcmc for extra countries (they can have different length than the other countries)
  par.names.cs <- mig.parameter.names.cs()
  for (country in country.index){
    country.obj <- get.country.object(country, mcmc.set$meta, index=TRUE)
    for (par in par.names.cs) {
      values <- get.mig.parameter.traces.cs(mcmc.set$mcmc.list, country.obj, par, 
                                            burnin=burnin)
      selected.simu <- get.thinning.index(nr.points, dim(values)[1])
      if (length(selected.simu$index) < nr.points)
        selected.simu$index <- sample(selected.simu$index, nr.points, replace=TRUE)
      values <- values[selected.simu$index,]
      write.values.into.file.cdep(par, values, dir, country.code=country.obj$code, 
                                  compression.type='None')
    }
  }
  if(verbose) cat('done.\n')
}

get.thinning.index <- function(nr.points, all.points) {
  if (!is.null(nr.points)) {
    nr.points <- ifelse(nr.points >= all.points, all.points, nr.points)
  } else {
    nr.points <- all.points
  }
  if (nr.points > 0) {
    step <- all.points/nr.points
    idx <- floor(seq(floor(step), all.points, by=step))
  } else idx<-NULL
  return(list(nr.points=nr.points, index=idx))
}

"get.nr.countries" <- function(meta, ...) UseMethod("get.nr.countries")

get.nr.countries.bayesMig.mcmc.meta <- function(meta, ...) 
  return(meta$nr.countries)

"get.nr.countries.est" <- function(meta, ...) UseMethod("get.nr.countries.est")

get.nr.countries.est.bayesMig.mcmc.meta <- function(meta, ...) 
  return(meta$nr.countries)

country.names <- function(meta) {
  #JA: This assumes we produce projections for all listed countries.
  return(meta$fullCountryNameVec)
}


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
