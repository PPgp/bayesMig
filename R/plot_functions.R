##########
#MIGRATION
##########

get.trajectories <- function(mig.pred, country, nr.traj=NULL, base.name='traj') {
  traj.file <- file.path(mig.pred$output.directory, paste(base.name, '_country', country, '.rda', sep=''))
  if (!file.exists(traj.file)) return(list(trajectories=NULL))
  load(traj.file)
  thintraj <- get.thinning.index(nr.traj, dim(trajectories)[2]) 
  if (thintraj$nr.points == 0) return(list(trajectories=NULL))
  traj.idx <- thintraj$index

  if(!is.null(trajectories)) {
    rownames(trajectories) <- get.prediction.years(mig.pred$mcmc.set$meta, dim(trajectories)[1])
  }
  return(list(trajectories=trajectories, index=traj.idx))
}


get.quantile.from.prediction <- function(mig.pred, quantile, country.index) {
  quant.values <- mig.pred$quantiles[country.index, as.character(quantile),]
  return(quant.values)
}

get.median.from.prediction <- function(mig.pred, country.index) {
  return(get.quantile.from.prediction(mig.pred, quantile=0.5, country.index=country.index))
}

get.traj.quantiles <- function(mig.pred, country.index, country.code, trajectories=NULL, pi=80) {
  al <- (1-pi/100)/2
  quantile.values <- as.numeric(dimnames(mig.pred$quantiles)[[2]])
  alidx<-round(quantile.values,6)==round(al,6)
  cqp <- NULL
  if (any(alidx)) { # pre-saved quantiles
    alidx2 <- round(quantile.values,6)==round(1-al,6)
    cqp <- rbind(mig.pred$quantiles[country.index, alidx,], 
                 mig.pred$quantiles[country.index, alidx2,])
  } else { # non-standard quantiles
    reload <- FALSE
    if (is.null(trajectories)) {
      if(mig.pred$nr.traj > 0) reload <- TRUE
    } else { 
      if (dim(trajectories)[2] < 2000 && mig.pred$nr.traj > dim(trajectories)[2]) reload <- TRUE
    }
    if(reload) {
      #load 2000 trajectories maximum for computing quantiles
      traj.reload <- get.trajectories(mig.pred, mig.pred$mcmc.set$meta$regions$country_code[country.index], 2000)
      trajectories <- traj.reload$trajectories
    }
    if (!is.null(trajectories)) {
      cqp <- apply(trajectories, 1, 
                   quantile, c(al, 1-al), na.rm = TRUE)
    }
  }
  return(cqp)
}

#' @title Graphical output of posterior distribution of migration trajectories
#'
#' @description This function plots the posterior distribution of migration trajectories for all countries,
#' including their median and given probability intervals
#' 
#' @usage mig.trajectories.plot.all(mig.pred, output.dir=file.path(getwd(), 'migTrajectories'),
#' output.type="png", main=NULL, verbose=FALSE, ...)
#' 
#' @param mig.pred Prediction object of class \code{bayesMig.prediction}
#' @param output.dir Directory into which resulting plots are written
#' @param output.type Type of the resulting plot files. Can be "png", "pdf", "jpeg", "bmp",
#' "tiff", or "postscript"
#' @param main Main title for the plots. Any occurrence of the string "XXX" is replaced by the
#' name of the appropriate country.
#' @param verbose Logical value. Switches log messages on and off.
#' @param ... Other variables passed to called functions


mig.trajectories.plot.all <- function(mig.pred, 
                                     output.dir=file.path(getwd(), 'migTrajectories'),
                                     output.type="png", verbose=FALSE, ...) {
  
  # plots e0 trajectories for all countries
  bayesTFR:::.do.plot.all(mig.pred$mcmc.set$meta, output.dir, mig.trajectories.plot, output.type=output.type, 
                          file.prefix='Migplot', plot.type='Mig graph', verbose=verbose, mig.pred=mig.pred, ...)
}


mig.trajectories.plot <- function(mig.pred, country, pi=c(80, 95), 
                                  nr.traj=10,
                                  xlim=NULL, ylim=NULL, type='b', 
                                  xlab='Year', ylab='Migration rate', main=NULL, lwd=c(2,2,2,1), 
                                  col=c('black', 'red', 'red','#00000020'),
                                  show.legend=TRUE, add=FALSE, ...
) {
  # lwd/col is a vector of 4 line widths/colors for: 
  #	1. observed data, 2. median, 3. quantiles, 4. trajectories
  if (missing(country)) {
    stop('Argument "country" must be given.')
  }
  country <- get.country.object(country, mig.pred$mcmc.set$meta)
  mig_observed <- mig.pred$mcmc.set$meta$mig.rates[country$index,]
  #JA: Assume present year is 2015, observations cover 1950--2015 and predictions cover 2015--2100.
  lpart1 <- length(mig_observed)
  # y1.part2 <- NULL
  # lpart2 <- min(tfr.pred$mcmc.set$meta$T_end, tfr.pred$present.year.index) - T_end_c[country$index] + suppl.T
  # if (lpart2 > 0) {
  #   p2idx <- (T_end_c[country$index]+1-suppl.T):nrow(tfr_matrix_reconstructed)
  #   y1.part2 <- tfr_matrix_reconstructed[p2idx,country$index]
  #   names(y1.part2) <- rownames(tfr_matrix_reconstructed)[p2idx]
  # }

  x1 <- as.integer(names(mig_observed))
  x2 <- as.numeric(dimnames(mig.pred$quantiles)[[3]])
  trajectories <- get.trajectories(mig.pred, country$code, nr.traj=nr.traj)
  # plot historical data: observed
  if (!add) {
    if(is.null(xlim)) xlim <- c(min(x1,x2), max(x1,x2))
    if(is.null(ylim)) ylim <- c(min(trajectories$trajectories, mig_observed, mig.pred$quantiles[country$index,,]),
                                max(trajectories$trajectories, mig_observed, mig.pred$quantiles[country$index,,]))
    if(is.null(main)) main <- country$name
    plot(xlim, ylim, type='n', xlim=xlim, ylim=ylim, ylab=ylab, xlab=xlab, main=main, 
         panel.first = grid())
  }
  points.x <- x1
  points.y <- mig_observed
  points(points.x, points.y, type=type, lwd=lwd[1], col=col[1])

  # plot trajectories
  if(!is.null(trajectories$trajectories)) { 
    for (i in 1:length(trajectories$index)) {
      lines(x2, trajectories$trajectories[,trajectories$index[i]], type='l', col=col[4], lwd=lwd[4])
    }
  }
  # plot median
  mig.median <- get.median.from.prediction(mig.pred, country$index)
  lines(x2, mig.median, type='l', col=col[3], lwd=lwd[3]) 
  # plot given CIs
  lty <- 2:(length(pi)+1)
  for (i in 1:length(pi)) {
    cqp <- get.traj.quantiles(mig.pred, country$index, country$code, trajectories$trajectories, pi[i])
    if (!is.null(cqp)) {
      lines(x2, cqp[1,], type='l', col=col[3], lty=lty[i], lwd=lwd[3])
      lines(x2, cqp[2,], type='l', col=col[3], lty=lty[i], lwd=lwd[3])
    }
  }
  legend <- c()
  cols <- c()
  lwds <- c()
  lty <- c(1, lty)
  median.legend <- 'median'
  legend <- c(legend, median.legend, paste(pi, '% PI', sep=''))
  cols <- c(cols, col[2], rep(col[3], length(pi)))
  lwds <- c(lwds, lwd[2], rep(lwd[3], length(pi)))
  if(show.legend) {
    legend <- c(legend, 'observed migration')
    cols <- c(cols, col[1])
    lty <- c(lty, 1)
    pch <- c(rep(-1, length(legend)-1), 1)
    lwds <- c(lwds, lwd[1])
    legend('bottomleft', legend=legend, lty=lty, bty='n', col=cols, pch=pch, lwd=lwds)
  }
}

mig.trajectories.table <- function(mig.pred, country, pi=c(80, 95), ...) {
  return(tfr.trajectories.table(mig.pred, country=country, pi=pi, half.child.variant = FALSE, ...))
}


mig.partraces.plot <- function(mcmc.list=NULL, sim.dir=file.path(getwd(), 'bayesMig.output'), 
                               chain.ids=NULL, par.names=mig.parameter.names(), 
                               nr.points=NULL, dev.ncol=3, ...) {
  if (is.null(mcmc.list))
    mcmc.list <- get.mig.mcmc(sim.dir)
  bayesTFR:::do.plot.tfr.partraces(mcmc.list, load.mig.parameter.traces, chain.ids=chain.ids, 
                        nr.points=nr.points, par.names=par.names, dev.ncol=dev.ncol, ...)
}

mig.partraces.cs.plot <- function(country, mcmc.list=NULL, sim.dir=file.path(getwd(), 'bayesMig.output'),
                                  chain.ids=NULL, par.names=mig.parameter.names.cs(),
                                  nr.points=NULL, dev.ncol=3, ...) {
  if (is.null(mcmc.list))
    mcmc.list <- get.mig.mcmc(sim.dir)
  mcmc.list <- get.mcmc.list(mcmc.list)
  country.obj <- get.country.object(country, mcmc.list[[1]]$meta)
  if (is.null(country.obj$name))
    stop('Country ', country, ' not found.')
  bayesTFR:::do.plot.tfr.partraces(mcmc.list, load.mig.parameter.traces.cs, 
                        main.postfix=paste0('(',country.obj$name,')'), chain.ids=chain.ids, nr.points=nr.points, 
                        country=country.obj$code, par.names=par.names, dev.ncol=dev.ncol, ...)
}

