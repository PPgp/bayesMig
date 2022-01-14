##########
#MIGRATION
##########


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
  trajectories <- bayesTFR:::get.trajectories(mig.pred, country$code, nr.traj=nr.traj)
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
    cqp <- bayesTFR:::get.traj.quantiles(mig.pred, country$index, country$code, trajectories$trajectories, pi[i])
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
                               nr.points=NULL, dev.ncol=2, ...) {
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

