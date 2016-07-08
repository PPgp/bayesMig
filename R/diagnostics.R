###############
#MIGRATION
###############

mig.raftery.diag <- function(mcmc=NULL, 
                             sim.dir=file.path(getwd(), 'bayesMig.output'),
                             burnin=0, country=NULL,
                             par.names = mig.parameter.names(),
                             par.names.cs =mig.parameter.names.cs(),
                             country.sampling.prop=1,
                             verbose=TRUE, ...
) {
  is.error <- function(rd) {
    if(is.null(dim(rd[[1]]$resmatrix))) {
      print(rd)
      return(TRUE)
    }
    return(FALSE)
  }
  if(verbose) cat('\nComputing Raftery Diagnostics ... this can take a while ...\n')
  if (is.null(mcmc)) {
    mcmc.set <- get.mig.mcmc(sim.dir=sim.dir)
  } else {
    mcmc.set <- mcmc
  }
  gui.option.name <- paste('bDem', strsplit(class(mcmc.set), '.', fixed=TRUE)[[1]][1], 'diagnose', sep='.')
  gui.option.name.status <- paste(gui.option.name, 'status', sep='.')
  gui.options <- list()
  result.025 <- result.975 <- burnin.025 <- burnin.975 <- rd.025<- thin.ind.025 <- thin.ind.975 <- NULL
  if (!is.null(par.names)) {
    coda.mc <- coda.list.mcmc(mcmc.set, country=country, rm.const.pars=TRUE,
                              par.names=par.names, par.names.cs=NULL, burnin=burnin, ...
    )
    gui.options[[gui.option.name.status]] <- 'processing country-independent pars (q=0.025)'
    unblock.gtk(gui.option.name, gui.options)
    if(verbose) cat('\t\tProcessing raftery.diag(..., r=0.0125, q=0.025) for country-independent parameters\n')
    rd.025 <- raftery.diag(coda.mc, r=0.0125, q=0.025)
    if(is.error(rd.025)) return()
    thin.ind.025 <- diag.thin.indep(coda.mc, q=0.025)
    colnames(thin.ind.025) <- rownames(rd.025[[1]]$resmatrix)
    gui.options[[gui.option.name.status]] <- 'processing country-independent pars (q=0.975)'
    unblock.gtk(gui.option.name, gui.options)
    if(verbose) cat('\t\tProcessing raftery.diag(..., r=0.0125, q=0.975) for country-independent parameters\n')
    rd.975 <- raftery.diag(coda.mc, r=0.0125, q=0.975)
    if(is.error(rd.975)) return()
    thin.ind.975 <- diag.thin.indep(coda.mc, q=0.975)
    colnames(thin.ind.975) <- rownames(rd.025[[1]]$resmatrix)
    result.025 <- result.975 <- burnin.025 <- burnin.975 <- matrix(NA, 
                                                                   nrow=length(rd.025), ncol=dim(rd.025[[1]]$resmatrix)[1],
                                                                   dimnames=list(c(), rownames(rd.025[[1]]$resmatrix)))
    for(i in 1:length(rd.025)) {
      result.025[i,] <- rd.025[[i]]$resmatrix[,'N']
      result.975[i,] <- rd.975[[i]]$resmatrix[,'N']
      burnin.025[i,] <- rd.025[[i]]$resmatrix[,'M']
      burnin.975[i,] <- rd.975[[i]]$resmatrix[,'M']
    }
  }
  c.index <- NULL
  if(!is.null(par.names.cs)) {
    if (!is.null(country)) {
      country.obj <- get.country.object(country, mcmc.set$meta)
      c.index <- country.obj$index
    } else {
      c.index <- get.countries.index(mcmc.set$meta)
      if(country.sampling.prop < 1) 
        c.index <- sort(sample(c.index, size=round(length(c.index)*country.sampling.prop,0)))
    }
    if(verbose) 
      cat('\t\tProcessing raftery.diag for country ')
    country.counter <- 0
    status.for.gui <- paste('out of', length(c.index), 'countries.')
    for(country.idx in c.index) {
      if(getOption(gui.option.name, default=FALSE)) {
        # This is to unblock the GUI, if the run is invoked from bayesDem
        # and pass info about its status
        # In such a case the gtk libraries are already loaded
        country.counter <- country.counter + 1
        gui.options[[gui.option.name.status]] <- paste('finished', country.counter, status.for.gui)
        unblock.gtk(gui.option.name, gui.options)
      }
      country.obj <- get.country.object(country.idx, mcmc.set$meta, index=TRUE)
      coda.mc.cs <- coda.list.mcmc(mcmc.set, country=country.obj$code, 
                                   rm.const.pars=TRUE,
                                   par.names = NULL,
                                   par.names.cs=par.names.cs, burnin=burnin, ...
      )
      if(verbose) 
        cat(country.idx, ', ')
      
      #If all parameters are constant for that country (as they are for Montserrat),
      #then coda.mc.cs is an empty list and we should just move on to the next country.
      if(length(coda.mc.cs)==0){
        next
      }
      rd.025 <- raftery.diag(coda.mc.cs, r=0.0125, q=0.025)
      if(is.error(rd.025)) return()
      thin.ind.025 <- cbind(thin.ind.025, diag.thin.indep(coda.mc.cs, q=0.025))
      npar <- dim(rd.025[[1]]$resmatrix)[1]
      ncols <- ncol(thin.ind.025)
      colnames(thin.ind.025)[(ncols-npar+1):ncols] <- rownames(rd.025[[1]]$resmatrix)
      rd.975 <- raftery.diag(coda.mc.cs, r=0.0125, q=0.975)
      if(is.error(rd.975)) return()
      thin.ind.975 <- cbind(thin.ind.975, diag.thin.indep(coda.mc.cs, q=0.975))
      colnames(thin.ind.975)[(ncols-npar+1):ncols] <- rownames(rd.025[[1]]$resmatrix)
      
      m <- matrix(NA, nrow=length(rd.025), ncol=npar, 
                  dimnames=list(c(), rownames(rd.025[[1]]$resmatrix)))
      result.025 <- cbind(result.025, m)
      result.975 <- cbind(result.975, m)
      burnin.025 <- cbind(burnin.025, m)
      burnin.975 <- cbind(burnin.975, m)
      for(i in 1:length(rd.025)) {
        idx <- (ncol(result.025)-npar+1):ncol(result.025)
        result.025[i, idx] <- rd.025[[i]]$resmatrix[,'N']
        result.975[i, idx] <- rd.975[[i]]$resmatrix[,'N']
        burnin.025[i, idx] <- rd.025[[i]]$resmatrix[,'M']
        burnin.975[i, idx] <- rd.975[[i]]$resmatrix[,'M']
      }
      gc()
    }
    if(verbose) 
      cat('\n')		
  }
  if(is.null(rd.025)) {
    if(verbose) cat('\t\tNothing to be done.\n')
    return()
  }
  if(verbose) cat('\t\tProcessing results ...')
  nr.chains <- length(rd.025)
  nr.par <- ncol(result.025)
  par.names.all <- colnames(result.025)
  if (is.null(par.names.all)) { # probably an error occured in raftery.diag
    print(rd.025)
    return()
  }
  
  colidx.cind <- get.full.par.names(par.names, par.names.all, index=TRUE)
  colidx.cdep <- get.full.par.names.cs(par.names.cs, par.names.all, index=TRUE)
  
  lcolidx.cind <- length(colidx.cind)
  lcolidx.cdep <- length(colidx.cdep)
  N.country.ind <- N.country.dep <- notconv1 <- notconv2 <- notconvinchain1 <- notconvinchain2 <- NULL
  iter <- nr.chains*(mcmc.set$mcmc.list[[1]]$finished.iter - burnin)
  for(chain in 1:nr.chains) {
    N <- rbind(result.025[chain,], result.975[chain,])
    where.larger <- N > iter
    notconv1 <- rbind(notconv1, 
                      cbind(parameter.name=par.names.all[where.larger[1,]], 
                            chain.id=rep(chain, sum(where.larger[1,])),
                            N=N[1,where.larger[1,]])
    )
    notconv2 <- rbind(notconv2, 
                      cbind(parameter.name=par.names.all[where.larger[2,]], 
                            chain.id=rep(chain, sum(where.larger[2,])),
                            N=N[2,where.larger[2,]])
    )
    where.larger.inchain <- N > (mcmc.set$mcmc.list[[chain]]$finished.iter - burnin)
    notconvinchain1 <- rbind(notconvinchain1, 
                             cbind(parameter.name=par.names.all[where.larger.inchain[1,]], 
                                   chain.id=rep(chain, sum(where.larger.inchain[1,])),
                                   N=N[1,where.larger.inchain[1,]])
    )
    notconvinchain2 <- rbind(notconvinchain2, 
                             cbind(parameter.name=par.names.all[where.larger.inchain[2,]], 
                                   chain.id=rep(chain, sum(where.larger.inchain[2,])),
                                   N=N[2,where.larger.inchain[2,]])
    )
    
    N.country.ind <- rbind(N.country.ind,
                           cbind(parameter.name=par.names.all[colidx.cind],
                                 chain.id=rep(chain, lcolidx.cind),
                                 N0.025=N[1,colidx.cind],
                                 N0.975=N[2,colidx.cind]))
    N.country.dep <- rbind(N.country.dep,
                           cbind(parameter.name=par.names.all[colidx.cdep],
                                 chain.id=rep(chain, lcolidx.cdep),
                                 N0.025=N[1,colidx.cdep],
                                 N0.975=N[2,colidx.cdep]))
  }
  notconv <- list(notconv1, notconv2)
  notconvinchain <- list(notconvinchain1, notconvinchain2)
  names(notconv) <- names(notconvinchain) <- c('0.025', '0.975')
  # change data types
  for (i in 1:2) {
    if (nrow(notconv[[i]]) > 0) {
      notconv[[i]] <- data.frame(notconv[[i]], row.names=1:nrow(notconv[[i]]))
    }else{notconv[[i]] <- NA}
    if (nrow(notconvinchain[[i]]) > 0) {
      notconvinchain[[i]] <- data.frame(notconvinchain[[i]], row.names=1:nrow(notconvinchain[[i]]))
    }else{notconvinchain[[i]] <- NA}
  }
  if (!is.null(par.names)) {
    N.country.ind <- data.frame(N.country.ind, row.names=1:nrow(N.country.ind))
  } else {N.country.ind <- NULL}
  if (!is.null(par.names.cs)) {
    N.country.dep <- data.frame(N.country.dep, row.names=1:nrow(N.country.dep))
  } else {N.country.dep <- NULL}
  
  #compute medians over chains
  s <- rbind(apply(result.025, 2, median), apply(result.975, 2, median))
  bi <- rbind(apply(burnin.025, 2, median), apply(burnin.975, 2, median))
  thin.ind <- rbind(apply(thin.ind.025, 2, median), apply(thin.ind.975, 2, median))
  
  colnames(s) <- par.names.all
  
  full.par.names <- get.full.par.names(par.names, par.names.all)
  full.par.names.cs <- get.full.par.names.cs(par.names.cs, par.names.all)
  if (is.null(country) && !is.null(full.par.names.cs)) {
    short.par.names.cs <- strsplit(full.par.names.cs, "_country.*") # removes the '_country*' postfix
    short.par.names.cs <- unique(short.par.names.cs)
  } else short.par.names.cs <- full.par.names.cs
  qcs <- qbics <- matrix(NA, nrow=2, ncol=length(short.par.names.cs)+length(full.par.names),
                         dimnames=list(c(), c(short.par.names.cs, full.par.names)))
  if (is.null(country)) {
    # compute maximum for mu_c, phi_c, and sigma2_c over countries
    for (pname in short.par.names.cs) {
      colidx <- grep(paste('^', pname, '_.*[[:digit:]]$', sep=''), par.names.all)
      for (j in 1:2) {
        qcs[j,pname] <- max(s[j,colidx])
        qbics[j,pname] <- max(bi[j,colidx])
      }
    }
  } else {
    qcs[,full.par.names.cs] <- s[, full.par.names.cs]
    qbics[,full.par.names.cs] <- bi[, full.par.names.cs]
  }
  Nmed.cs <- s[, full.par.names.cs, drop=FALSE]
  qcs[,full.par.names] <- s[, full.par.names]
  qbics[,full.par.names] <- bi[, full.par.names]
  if(verbose) cat(' done.\n')
  return(list(Nmedian=round(qcs,0),
              burnin=round(qbics,0),
              not.converged.parameters=notconv,
              not.converged.inchain.parameters=notconvinchain,
              N.country.indep=N.country.ind,
              N.country.spec=N.country.dep,
              Nmedian.country.spec=round(Nmed.cs,0),
              thin.ind=list('0.025'=thin.ind.025, '0.975'=thin.ind.975, median=thin.ind),
              nr.countries=c(used=length(c.index), total=get.nr.countries.est(mcmc.set$meta))))
}

diag.thin.indep <- function(mcmc.list, q) {
  k <- matrix(NA, nrow=length(mcmc.list), ncol=ncol(mcmc.list[[1]]))
  for(imcmc in 1:length(mcmc.list)) {
    mcmc <- mcmc.list[[imcmc]]
    for (i in 1:ncol(mcmc)) 
      k[imcmc,i] <- thinindep(mcmc[,i], q)	}
  return(k)
}

thinindep <- function(x,q){
  ## find the smallest integer k that makes the thinned chain independent ## 
  ## x is the MCMC samples, and q is the quantile being estimated ##
  k=0
  bic=0
  n0=length(x)
  qx=quantile(x,q)
  while(bic>=0){
    k=k+1
    u=x[seq(1,n0,by=k)]
    n=length(u)
    ## change the continuous chain to 0-1 chain ##
    zt <- factor(u<=qx, levels=c(FALSE, TRUE))
    ## calculate the 2*2 contingency table ##
    nij <- table(zt[-n], zt[-1])
    nc=apply(nij, 2, sum)
    nr=apply(nij, 1, sum)
    ni=c(nc[1],nc[1],nc[2],nc[2])
    nj=c(nr, nr)
    nij <- c(nij)
    logn <- log(n)
    idx <- nij > 0
    lognij <- log(nij[idx])
    logni <- log(ni[idx])
    lognj <- log(nj[idx])
    ## use BIC to determine whether thinning every kth sample is enough ##
    bic <- sum(2*(nij[idx]*(logn+lognij-logni-lognj))) - logn
  }
  return(k)
}

.do.diagnose <- function(type, class.name, sim.dir, thin=80, burnin=2000, express=FALSE, 
                         country.sampling.prop=NULL, keep.thin.mcmc=FALSE, verbose=TRUE) {
  get.country.name <- function(par) {
    cindex <- strsplit(par, '_country')[[1]]
    cindex <- as.numeric(cindex[length(cindex)])
    return(get.country.object(cindex, meta=mcmc.set$meta)$name)
  }
  
  mcmc.set <- do.call(paste('get.', type, '.mcmc', sep=''), 
                      list(sim.dir=sim.dir))
  if(is.null(mcmc.set))
    stop('No valid simulation in ', sim.dir)
  iter <- get.total.iterations(mcmc.set$mcmc.list, burnin)
  if(iter <= 0) stop('0 number of iterations. Check the value of burnin.')
  if(iter/thin < 1) stop(paste('Value of thin is too high (', thin, 
                               ') given the total number of iterations (', iter, ').', sep=''))
  #run raftery.diag on country-independent parameters
  diag.procedure <- paste(type, '.raftery.diag', sep='')
  raftery.diag.res <- do.call(diag.procedure, 
                              list(mcmc.set, par.names.cs=NULL, thin=thin,
                                   burnin=burnin, verbose=verbose))
  if (is.null(raftery.diag.res)) stop(paste('Problem in', diag.procedure))
  raftery.diag.res.cs <- NULL
  if(!is.null(country.sampling.prop)) express <- FALSE
  if(!express &&((!is.null(country.sampling.prop) && (country.sampling.prop>0)) 
                 || is.null(country.sampling.prop))) {
    #run raftery.diag on country-specific parameters
    raftery.diag.res.cs <- do.call(diag.procedure, list(mcmc.set, 
                                                        par.names = NULL,
                                                        thin=thin,  burnin=burnin,
                                                        country.sampling.prop=if(is.null(country.sampling.prop)) 1 else country.sampling.prop,
                                                        verbose=verbose
    ))
    if (is.null(raftery.diag.res.cs)) stop(paste('Problem in', diag.procedure))
  }
  status <- c(red=FALSE, green=FALSE)
  res <- NULL
  lres.cind <- 0
  if (!is.null(raftery.diag.res)) {
    res <- process.not.converged.parameters(raftery.diag.res, iter)
    lres.cind <- nrow(res)
  }
  if (!is.null(raftery.diag.res.cs)) {
    res <- rbind(res, process.not.converged.parameters(raftery.diag.res.cs, iter))
  } 
  to.run <- 0
  if(nrow(res) > 0) {
    max.idx <- which.max(res[,2])
    to.run <- res[max.idx,2]
  }
  if(to.run <= 0) status['green'] <- TRUE
  else status['red'] <- TRUE
  nr.countries <- if(!is.null(raftery.diag.res.cs)) raftery.diag.res.cs$nr.countries 
  else c(used=0, total=get.nr.countries.est(mcmc.set$meta))
  use.nr.traj <- floor(iter/thin)
  thinned.mcmc <- NULL
  if(keep.thin.mcmc) {
    thinned.mcmc <- do.call(paste('get.thinned.', type, '.mcmc', sep=''), 
                            list(mcmc.set, thin=thin, burnin=burnin))
    if(is.null(thinned.mcmc) || thinned.mcmc$meta$parent.iter < iter) 
      thinned.mcmc <- do.call(paste('create.thinned.', type, '.mcmc', sep=''), 
                              list(mcmc.set, thin=thin, burnin=burnin, verbose=verbose))
  }
  diag <- structure(list(result=res,
                         lresult.country.independent=lres.cind,
                         country.independent=raftery.diag.res, 
                         country.specific=raftery.diag.res.cs, 
                         iter.needed=to.run,
                         iter.total=iter, use.nr.traj=use.nr.traj,
                         status=status,
                         mcmc.set=mcmc.set,
                         thin.mcmc=thinned.mcmc,
                         burnin=burnin,
                         thin=thin,
                         express=express, 
                         nr.countries=nr.countries),
                    class=class.name)
  if(verbose) summary(diag)
  save.dir <- file.path(mcmc.set$meta$output.dir, 'diagnostics')
  if(!file.exists(save.dir)) 
    dir.create(save.dir, recursive=TRUE)
  save.file <- do.call(paste('store.', class.name, sep=''), list(diag, thin=thin, burnin=burnin, 
                                                                 output.dir=save.dir))
  if(verbose) cat('\nConvergence diagnostics stored in', save.file, '\n')
  return(diag)
}


process.not.converged.parameters <- function(diag, iter) {
  diff <-apply(diag$Nmedian, 2, max) - iter
  posdiff <- diff[diff > 0]
  not.converged <- names(posdiff)
  lnot.converged <- length(not.converged)
  N <- matrix(0, nrow=lnot.converged, ncol=2, dimnames=list(not.converged,
                                                            c('Total iterations needed', 'Remaining iterations')))
  if(lnot.converged == 0) return(N)
  for (i in 1:2) {
    is.contained <- is.element(not.converged, 
                               colnames(diag$Nmedian))
    parnames <- not.converged[is.contained]
    for (parname in parnames) {
      N[parname,1] <- max(N[parname,1], diag$Nmedian[i, parname])
    }
  }
  N[,2] <- pmax(N[,1] - iter, 0)
  return(N)
}

mig.diagnose <- function(sim.dir, thin=80, burnin=2000, express=FALSE, 
                         country.sampling.prop=NULL, keep.thin.mcmc=FALSE, verbose=TRUE) {
  invisible(.do.diagnose(type='mig', class.name='bayesMig.convergence', 
                         sim.dir=sim.dir, thin=thin, burnin=burnin, express=express,
                         country.sampling.prop=country.sampling.prop, keep.thin.mcmc=keep.thin.mcmc,	verbose=verbose))
}


#########
#Test code
#########

#source("mcmc_ini.R")
#source("mcmc_update.R")
#source("mcmc_sampling.R")
#source("mcmc_storage.R")
#source("get_outputs.R")

#o=mig.raftery.diag(sim.dir=file.path(getwd(),'bayesMig.diagnosis'))

#Test. Something's going wrong with Montserrat because all of its parameters are constant.
#mcmc.set=get.mig.mcmc(sim.dir='bayesMig.output')
#get.country.object(175, mcmc.set$meta, index=TRUE)#Montserrat's code is 500.
#coda.mc.cs = coda.list.mcmc(mcmc.set, country=500, rm.const.pars=TRUE,
#                            par.names=NULL, par.names.cs=mig.parameter.names.cs(),
#                            burnin=0)
