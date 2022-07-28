##########
#MIGRATION
##########


#' @title Output of posterior distribution of migration trajectories
#'
#' @description The functions plot/tabulate the posterior distribution of trajectories of net migration rates
#'     for a given country, or for all countries, including their median and given probability 
#'     intervals.
#' 
#' @param mig.pred Prediction object of class \code{bayesMig.prediction}.
#' @param country Name or numerical code of a country. It can also be given as ISO-2 or ISO-3 characters.
#' @param pi Probability interval (as percentage) to be included in the output. It can be a single number or a vector.
#' @param nr.traj Number of trajectories to be plotted. If \code{NULL}, all trajectories are plotted, otherwise they are thinned evenly.
#' @param xlim,ylim,type,xlab,ylab Graphical parameters passed to the \code{\link{plot}} function.
#' @param main Main title for the plot(s). In \code{mig.trajectories.plot.all} any occurrence of the string 
#'     \dQuote{XXX} is replaced by the name of the appropriate country.
#' @param lwd,col Vector of four elements giving the line width and color for: 1. observed data, 
#'     2. median, 3. quantiles, 4. trajectories.
#' @param show.legend Logical controlling whether a legend should be drawn.
#' @param add Logical controlling whether the trajectories should be plotted into a new graphic 
#'     device (\code{FALSE}) or into an existing device (\code{TRUE}). One can use this argument to plot
#'     trajectories from multiple countries into one graphics.
#' @param scale Logical. If \code{TRUE}, values are scaled to be \dQuote{per population}, i.e. 
#'     they are divided by \code{pop.denom} passed to \code{\link{mig.predict}}.
#' @param \dots Additional graphical parameters. In addition, for \code{mig.trajectories.plot.all} 
#'     any of the arguments of \code{tfr.trajectories.plot} can be passed here.
#'     
#' @details \code{mig.trajectories.plot} plots posterior distribution of trajectories of net migration
#'     rates for a given country. \code{mig.trajectories.table} gives the same output as a table. 
#'     \code{mig.trajectories.plot.all} creates a set of graphs (one per country) that are stored in 
#'     \code{output.dir}.
#'     
#'     The median and given probability intervals are computed using all available trajectories. 
#'     Thus, \code{nr.traj} does not influence those values - it is used only to control the number 
#'     of trajectories in the graphs.
#'     
#' @seealso \code{\link{mig.predict}}, \code{\link{summary.bayesMig.prediction}}
#' @examples
#' # See example in ?mig.predict
#' 
#' @export
#' @rdname plot-traj
#' 
mig.trajectories.plot <- function(mig.pred, country, pi=c(80, 95), 
                                  nr.traj=50,
                                  xlim=NULL, ylim=NULL, type='b', 
                                  xlab='Year', ylab='Migration rate', main=NULL, lwd=c(2,2,2,1), 
                                  col=c('black', 'red', 'red','#00000020'),
                                  show.legend=TRUE, add=FALSE, scale = FALSE, ...
                              ) {
  # lwd/col is a vector of 4 line widths/colors for: 
  #	1. observed data, 2. median, 3. quantiles, 4. trajectories
  if (missing(country)) {
    stop('Argument "country" must be given.')
  }
  country <- get.country.object(country, mig.pred$mcmc.set$meta)
  mig_observed <- mig.pred$mcmc.set$meta$mig.rates[country$index,]
  if(scale) mig_observed <- mig_observed / mig.pred$mcmc.set$meta$prior.scaler
  
  lpart1 <- length(mig_observed)
  # missing values imputation is not implemented yet
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
  mig.median <- get.median.from.prediction(mig.pred, country$index)
  
  if(scale) { # scale to be interpreted as "per population"
    if(!is.null(trajectories$trajectories)) trajectories$trajectories <- trajectories$trajectories / mig.pred$mcmc.set$meta$prior.scaler
    mig.pred$quantiles <- mig.pred$quantiles / mig.pred$mcmc.set$meta$prior.scaler
    mig.median <- mig.median / mig.pred$mcmc.set$meta$prior.scaler
  }
  
  # plot historical data: observed
  if (!add) {
    if(is.null(xlim)) xlim <- c(min(x1,x2), max(x1,x2))
    if(is.null(ylim)) {
      ylim <- c(min(trajectories$trajectories, mig_observed, mig.pred$quantiles[country$index,,]),
                                max(trajectories$trajectories, mig_observed, mig.pred$quantiles[country$index,,]))
    }
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

#' @param output.dir Directory into which resulting plots are written
#' @param output.type Type of the resulting plot files. Can be "png", "pdf", "jpeg", "bmp",
#' "tiff", or "postscript"
#' @param verbose Logical value. Switches log messages on and off.
#' @export
#' @rdname plot-traj

mig.trajectories.plot.all <- function(mig.pred, 
                                      output.dir=file.path(getwd(), 'migTrajectories'),
                                      output.type="png", verbose=FALSE, ...) {
  
  # plots e0 trajectories for all countries
  bayesTFR:::.do.plot.all(mig.pred$mcmc.set$meta, output.dir, mig.trajectories.plot, output.type=output.type, 
                          file.prefix='Migplot', plot.type='Mig graph', verbose=verbose, mig.pred=mig.pred, ...)
}

#' @rdname plot-traj
#' @export
mig.trajectories.table <- function(mig.pred, country, pi=c(80, 95), ...) {
  return(tfr.trajectories.table(mig.pred, country=country, pi=pi, half.child.variant = FALSE, ...))
}

#' @title Plotting MCMC Parameter Traces
#'
#' @description Functions for plotting the MCMC parameter traces from the migration model.
#' 
#' @param mcmc.list List of \code{\link{bayesMig.mcmc}} objects, or an object of class
#'     \code{\link{bayesMig.mcmc.set}} or of class \code{bayesMig.prediction}. If it is \code{NULL}, the 
#'     traces are loaded from \code{sim.dir}.
#' @param sim.dir Directory with the MCMC simulation results. It is only used if \code{mcmc.list} is \code{NULL}.
#' @param chain.ids List of MCMC identifiers to be plotted. If it is \code{NULL}, all chains found in 
#'     \code{mcmc.list} or \code{sim.dir} are plotted.
#' @param par.names Names of parameters for which traces should be plotted. By default all 
#'     country-independent parameters are plotted if used within \code{mig.partraces.plot}, or 
#'     country-specific parameters are plotted if used within \code{mig.partraces.cs.plot}.
#' @param nr.points Number of points to be plotted. If \code{NULL}, all points are plotted, 
#'     otherwise the traces are thinned evenly.
#' @param dev.ncol Number of column for the graphics device. If the number of parameters is smaller
#'     than \code{dev.ncol}, the number of columns is automatically decreased.
#' @param \dots Additional graphical parameters. 
#'     
#' @details The functions plot MCMC traces either for country-independent parameters 
#'     (\code{mig.partraces.plot} or for country-specific parameters (\code{mig.partraces.cs.plot}, 
#'     one graph per parameter.  One can restrict it to specific chains by setting 
#'     the \code{chain.ids} argument, and to specific parameters by setting the \code{par.names} 
#'     argument.
#' @export
#' @rdname plot-traces
#' 

mig.partraces.plot <- function(mcmc.list=NULL, sim.dir=file.path(getwd(), 'bayesMig.output'), 
                               chain.ids=NULL, par.names=mig.parameter.names(), 
                               nr.points=NULL, dev.ncol=2, ...) {
  if (is.null(mcmc.list))
    mcmc.list <- get.mig.mcmc(sim.dir)
  bayesTFR:::do.plot.tfr.partraces(mcmc.list, load.mig.parameter.traces, chain.ids=chain.ids, 
                        nr.points=nr.points, par.names=par.names, dev.ncol=dev.ncol, ...)
}

#' @param country Name or numerical code of a country. It can also be given as ISO-2 or ISO-3 characters.
#' @export
#' @rdname plot-traces
#' 
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

#' @title Plotting MCMC Parameter Density
#'
#' @description Functions for plotting the density of the posterior distribution of the MCMC parameters from the migration model.
#' 
#' @param mcmc.list List of \code{\link{bayesMig.mcmc}} objects, or an object of class
#'     \code{\link{bayesMig.mcmc.set}} or of class \code{bayesMig.prediction}. If it is \code{NULL}, the 
#'     values are loaded from \code{sim.dir}.
#' @param sim.dir Directory with the MCMC simulation results. It is only used if \code{mcmc.list} is \code{NULL}.
#' @param chain.ids List of MCMC identifiers to be plotted. If it is \code{NULL}, all chains found in 
#'     \code{mcmc.list} or \code{sim.dir} are plotted.
#' @param par.names Names of parameters for which density should be plotted. By default all 
#'     country-independent parameters are plotted if used within \code{mig.pardensity.plot}, or 
#'     country-specific parameters are plotted if used within \code{mig.pardensity.cs.plot}.
#' @param burnin Number of iterations to be discarded from the beginning of each chain before 
#'     computing the density.
#' @param dev.ncol Number of column for the graphics device. If the number of parameters is smaller
#'     than \code{dev.ncol}, the number of columns is automatically decreased.
#' @param low.memory Logical indicating if the processing should run in a low-memory mode. If it is 
#'     \code{FALSE}, traces of all available parameters are loaded into memory. Otherwise, parameters are
#'     loaded as they are needed.
#' @param \dots Further arguments passed to the \code{\link{density}} function.
#'     
#' @details The functions plot the density of the posterior distribution either for 
#'     country-independent parameters (\code{mig.pardensity.plot} or for country-specific 
#'     parameters (\code{mig.pardensity.cs.plot}, one graph per parameter.  
#'     One can restrict it to specific chains by setting the \code{chain.ids} argument and to specific 
#'     parameters by setting the \code{par.names} argument. 
#'     
#'     If \code{mcmc.list} is an object of class \code{\link{bayesMig.prediction}} 
#'     and if this object contains thinned traces, they are used instead of the full chains. 
#'     In such a case, \code{burnin} and \code{chain.ids} cannot be modified - their value is set 
#'     to the one used when the thinned traces were created, namely when running 
#'     \code{\link{mig.predict}}. In a situation with long MCMC chains, this approach can  
#'     significantly speed-up creation of the density plots.
#' @export
#' @rdname plot-density
#' 
mig.pardensity.plot <- function(mcmc.list=NULL, sim.dir = file.path(getwd(), 'bayesMig.output'), 
                               chain.ids = NULL, par.names = mig.parameter.names(), 
                               burnin = NULL, dev.ncol = 2, low.memory = TRUE, ...) {
  if (is.null(mcmc.list))
    mcmc.list <- get.mig.mcmc(sim.dir, low.memory = low.memory)
  bayesTFR:::do.plot.tfr.pardensity(mcmc.list, get.mig.parameter.traces, chain.ids = chain.ids, par.names = par.names,
                                    par.names.ext = par.names, burnin = burnin, dev.ncol = dev.ncol, ...)
}

#' @param country Name or numerical code of a country. It can also be given as ISO-2 or ISO-3 characters.
#' @export
#' @rdname plot-density
#' 
mig.pardensity.cs.plot <- function(country, mcmc.list = NULL, sim.dir = file.path(getwd(), 'bayesLife.output'), 
                                  chain.ids = NULL, par.names = mig.parameter.names.cs(), 
                                  burnin = NULL, dev.ncol = 3, low.memory = TRUE, ...) {
  if (is.null(mcmc.list))
    mcmc.list <- get.mig.mcmc(sim.dir, low.memory=low.memory)
  mcmc.l <- get.mcmc.list(mcmc.list)
  country.obj <- get.country.object(country, mcmc.l[[1]]$meta)
  if (is.null(country.obj$name))
    stop('Country ', country, ' not found.')
  bayesTFR:::do.plot.tfr.pardensity(mcmc.list, get.mig.parameter.traces.cs, chain.ids = chain.ids, par.names = par.names,
                                    par.names.ext = par.names,
                                    main.postfix = paste0('(',country.obj$name,')'),
                                    func.args = list(country.obj = country.obj),
                                    burnin = burnin, dev.ncol = dev.ncol, ...)
}


#' @export
.map.main.default.bayesMig.prediction <- function(pred, ...) 
  return('MIG: quantile')

#' @title World Map of Net Migration Rate
#' @description Generates a world map of the net migration rate for given quantile and 
#'     time period, which can be either projection or estimation time period. A map of 
#'     country-specific model parameters is also supported.
#' @param pred Object of class \code{\link{bayesMig.prediction}}. Note that location codes
#'     must correspond to the UN country codes in order to generate a world map.
#' @param \dots In \code{mig.map}, \dots are all arguments that can be passed 
#'     to \code{\link[bayesTFR]{tfr.map}}, such as \code{quantile}, \code{year}, 
#'     \code{projection.index}, \code{par.name}, \code{adjusted}, \code{device}, \code{main}, 
#'     \code{device.args}, and \code{data.args}. 
#'     In \code{mig.map.gvis}, \dots are all arguments that can be passed 
#'     to \code{\link[bayesTFR]{tfr.map.gvis}}. In addition, the first three functions accept 
#'     arguments passed to the \code{\link[rworldmap]{mapCountryData}} function of the \pkg{rworldmap} package.
#' @details \code{mig.map} creates a single map for the given time period and quantile. 
#'     \code{mig.map.all} generates a sequence of maps, namely one for each projection period. 
#'     If the package \pkg{fields} is installed, a color bar legend at the botom of the map is created.
#'     
#'     Function \code{get.mig.map.parameters} can be used in combination with \code{mig.map}. 
#'     (Note that \code{get.mig.map.parameters} is called from inside of \code{mig.map.all}.) 
#'     It sets breakpoints for the color scheme.
#'     
#'     Function \code{mig.map.gvis} creates an interactive map using the \pkg{googleVis} package 
#'     and opens it in an internet browser. It also generates a table of the mapped values that 
#'     can be sorted by columns interactively in the browser. 
#'     
#'     By default, both \code{mig.map} and \code{mig.map.gvis} produce maps of net migration rates. 
#'     Alternatively, the functions can be used to plot country-specific MCMC parameters into a world map. 
#'     They are given by the argument \code{par.name}. One can pass any value from 
#'     \code{\link{mig.parameter.names.cs}()}.
#'     
#' @seealso \code{\link[bayesTFR]{tfr.map}}
#' @rdname map
#' @export
#' 
mig.map <- function(pred, ...) {
  return(bayesTFR::tfr.map(pred, ...))
}

#' @export
#' @rdname map
mig.map.gvis <- function(pred, ...)
  bdem.map.gvis(pred, ...)


#' @param output.dir Directory into which resulting maps are stored.
#' @param output.type Type of the resulting files. It can be \dQuote{png}, \dQuote{pdf}, 
#'     \dQuote{jpeg}, \dQuote{bmp}, \dQuote{tiff}, or \dQuote{postscript}.
#' @param mig.range Range of the migration rate to be displayed. It is of the form 
#'     \code{c(}\var{mig.min}, \var{mig.max}\code{)}. By default, the whole available range is considered. 
#'     Note that countries with values outside of the given range will appear white.
#' @param nr.cats Number of color categories.
#' @param same.scale Logical controlling if maps for all years of this prediction object 
#'     should be on the same color scale.
#' @param quantile Quantile for which the map should be generated. It must be equal to one of the 
#'     values in \code{dimnames(pred$quantiles)[[2]]}, 
#'     i.e. 0, 0.025, 0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 0.95, 0.975, 1. 
#'     Value 0.5 corresponds to the median.
#' @param file.prefix Prefix for file names.
#' @export
#' @rdname map
mig.map.all <- function(pred, output.dir, output.type='png', mig.range=NULL, nr.cats=50, same.scale=TRUE, 
                       quantile=0.5, file.prefix='migwrldmap_', ...) {
  bayesTFR:::bdem.map.all(pred=pred, output.dir=output.dir, type='mig', output.type=output.type, range=mig.range,
                          nr.cats=nr.cats, same.scale=same.scale, quantile=quantile, file.prefix=file.prefix, ...)
}


#' @export
bdem.map.gvis.bayesMig.prediction <- function(pred, ...) {
  bayesTFR:::.do.gvis.bdem.map('mig', 'Net Migration Rate', pred, ...)
}

#' @param palette Color palette to use.
#' @export
#' @rdname map
#' 
get.mig.map.parameters <- function(pred, mig.range=NULL, nr.cats=50, same.scale=TRUE, 
                                   quantile=0.5, palette = "Blue-Red", ...) {
  map.pars <- list(pred=pred, quantile=quantile, ...)
  if (same.scale) {
    data <- pred$quantiles[,as.character(quantile),1]
    q <- if(is.null(mig.range)) c(min(data), max(data)) else mig.range
    quantiles <- seq(q[1], q[2], length=nr.cats-1)
    map.pars$catMethod <- quantiles
  } else {
    map.pars$numCats <- nr.cats
  }
  map.pars$colourPalette <- sapply(palette, hcl.colors, n = nr.cats)
  return(map.pars)
}


#' @export
par.names.for.worldmap.bayesMig.prediction <- function(pred, ...) {
  return(mig.parameter.names.cs())
}

#' @export
get.data.for.worldmap.bayesMig.prediction <- function(pred, ...)
  return(bayesTFR:::get.data.for.worldmap.bayesTFR.prediction(pred, ...))
