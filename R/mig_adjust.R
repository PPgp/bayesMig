#' @title Adjusting the Projection Medians
#' 
#' @description These functions are to be used by expert analysts. They allow to 
#'     change the projection medians either to specific values, or shift the medians 
#'     by a given constant or align one projection object with another.
#'     
#' @param sim.dir Directory containing the prediction object.
#' @param country Name or numerical code of a country.
#' @param values Vector of the new median values.
#' @param years Numeric vector giving years for which to change the median. 
#'     In \code{mig.median.set} it gives years which \code{values} correspond to.
#'     Ideally it should be of the same length as \code{values}. If it is \code{NULL}, 
#'     \code{values} are set starting from the first prediction time period. 
#'     If \code{values} correspond to consecutive years, only the first year might be given here.
#'     In \code{mig.align.predictions} it gives years for which the medians should be aligned.  
#' @param \dots Additional arguments passed to the underlying adjustment functions, such as 
#'     \code{verbose} to show/hide the progress of the adjustment.
#'     For \code{mig.shift.prediction.to.wpp} it can be \code{stat} with values \dQuote{median} (default) 
#'     or \dQuote{mean} to specify which statistics should be adjusted; 
#'     \code{wpp.year} to adjust it to if it differs from the wpp year of the simulation.
#' 
#' @details The function \code{mig.median.set} can be used to set the medians of the
#'     given country to specific values.
#'    
#' @return All functions return an updated object of class \code{\link{bayesMig.prediction}}.
#' @export
#' @rdname mig.adjust
mig.median.set <- function(sim.dir, country, values, years=NULL, ...) {
    pred <- get.mig.prediction(sim.dir)
    new.pred <- bayesTFR:::.bdem.median.set(pred, type='mig', country=country, 
                                            values=values, years=years, ...)
    store.bayesMig.prediction(new.pred)
    invisible(new.pred)
}

#' @export 
#' @keywords internal
#' @rdname internal
#' 
get.mig.shift <- function(country.code, pred) return(bayesTFR::get.tfr.shift(country.code, pred))

#' @param reset Logical. If \code{TRUE} medians in a range of \code{from} and \code{to} are
#'     reset to their original values.
#' @param shift Constant by which the medians should be offset. It is not used if \code{reset} is \code{TRUE}.
#' @param from Year from which the offset/reset should start. By default, it starts at the first prediction period.
#' @param to Year until which the offset/reset should be done. By default, it is set to the last prediction period.
#' 
#' @details Function \code{mig.median.shift} can be used to offset the medians by a specific constant, or to reset
#'     the medians to their original values. 
#' @export
#' @rdname mig.adjust
mig.median.shift <- function(sim.dir, country, reset = FALSE, shift = 0, 
                             from = NULL, to = NULL) {
    pred <- get.mig.prediction(sim.dir)
    new.pred <- bayesTFR:::.bdem.median.shift(pred, type='mig', country=country, reset=reset, 
                                              shift=shift, from=from, to=to)
    store.bayesMig.prediction(new.pred)
    invisible(new.pred)
}

#' @param countries Vector of country names or codes. If this argument is \code{NULL} (default), 
#'     the reset is done for all countries.
#'     
#' @details Function \code{mig.median.reset} resets medians of the given countries
#'     to the original values. By default it deletes adjustments for all countries.
#' @export
#' @rdname mig.adjust
mig.median.reset <- function(sim.dir, countries = NULL) {
    if(is.null(countries)) {
        pred <- get.mig.prediction(sim.dir)
        pred$median.shift <- NULL
        store.bayesMig.prediction(pred)
        cat('\nMedians for all countries reset.\n')
    } else
        for(country in countries) pred <- mig.median.shift(sim.dir, country, reset=TRUE)
        invisible(pred)
}

#' @param sim.dir1 Directory with the bayesMig prediction object to be adjusted.
#' @param sim.dir2 Directory with the bayesMig prediction object used to align the medians from \code{sim.dir1} to.
#' @param country.codes Numerical codes of countries to adjust. By default all countries 
#'     found in \code{sim.dir2} are adjusted in \code{sim.dir1}.
#' @details Function \code{mig.align.predictions} shifts medians stored in \code{sim.dir1} to match 
#'     the medians in \code{sim.dir1}.
#'     
#'     In all cases, if a median is modified, the corresponding offset is stored in the prediction object 
#'     (element \code{median.shift}). All functions write the updated prediction object back to disk. All
#'     functions in the package that use trajectories and trajectory statistics use the \code{median.shift} 
#'     values to offset the results correspondingly, i.e. trajectories are shifted the same way as the
#'     medians.
#' @export
#' @rdname mig.adjust
mig.align.predictions <- function(sim.dir1, sim.dir2, country.codes = NULL, 
                                  years = NULL, ...){
    pred1 <- get.mig.prediction(sim.dir1)
    pred2 <- get.mig.prediction(sim.dir2)
    cntries1 <- if(is.null(country.codes)) get.countries.table(pred1)$code else country.codes
    cntries2 <- get.countries.table(pred2)$code
    avail.years <- dimnames(pred2$quantiles)[[3]]
    if(is.null(years)) 
        years <- avail.years
    else years <- intersect(as.character(years), avail.years)
    for(cntry in cntries1){
        cidx <- which(cntries2 == cntry)
        if(length(cidx) == 0){
            warning("Country ", cntry, " not found in sim.dir2. No adjustment for this country.")
            next
        }
        adjust.to <- get.median.from.prediction(pred2, cidx, cntry)[years]
        pred <- mig.median.set(sim.dir1, cntry, values = adjust.to, 
                               years = as.integer(names(adjust.to)), ...)
    }
    invisible(pred)
}

.do.mig.shift.prediction.to.wpp <- function(pred, stat = "median", wpp.year = NULL, verbose = TRUE){
    country_code <- year <- NULL # for CRAN check not to complain
    meta <- pred$mcmc.set$meta
    wpp.year <- if(!is.null(wpp.year)) wpp.year else meta$wpp.year
    countries <- get.countries.table(pred$mcmc.set)$code
    if(wpp.year < 2024) stop("Function not implemented for WPP < 2024.")
    n <- if(meta$annual.simulation) "1" else "5"
    migcounts.wpp <- bayesTFR:::load.bdem.dataset(paste0("migproj", n, "dt"), wpp.year = wpp.year, verbose = verbose)
    popcounts.wpp <- bayesTFR:::load.bdem.dataset(paste0("popproj", n, "dt"), wpp.year = wpp.year, verbose = verbose)
    if(!meta$annual.simulation) migcounts.wpp[, year := year + 2] # align migration and pop years
    
    migrates.wpp <- merge(migcounts.wpp[, c("country_code", "name", "year", "mig"), with = FALSE],
                          popcounts.wpp[, c("country_code", "name", "year", "pop"), with = FALSE],
                      by = c("country_code", "name", "year")
                      )
    migrates.wpp$wpp <- migrates.wpp$mig / (migrates.wpp$pop - migrates.wpp$mig)
    if(!meta$annual.simulation) migrates.wpp[, year := year - 2]
    
    pred.years <- as.numeric(dimnames(pred$quantiles)[[3]])
    pred$median.shift <- NULL
    
    if(verbose) cat("\n")
    for(icntry in seq_along(countries)) {
        if (verbose) {
            if(interactive()) cat("\rAdjusting countries' prediction to wpp", wpp.year, " ... ", round(icntry/length(countries) * 100), ' %')
            else {
                if (icntry == 1)
                    cat("Adjusting countries' prediction to wpp", wpp.year, " ... ")
                cat(icntry, ", ")
            }
        }
        cntry <- countries[icntry]
        if(!stat %in% c("median", "mean")) stop("Argument 'stat' must be 'median' or 'mean', but is ", stat, ".")
        to.match <- merge(data.table::data.table(year = pred.years, median = pred$quantiles[icntry, "0.5", ]), 
                          migrates.wpp[country_code == cntry, c("year", "wpp"), with = FALSE], 
                          by = "year", all.x = TRUE)
        if(stat == "mean") to.match[["median"]] <- to.match[["median"]] + pred$traj.mean.sd[icntry, 1, ] - to.match[["median"]] # difference between the mean and median
        to.match$wpp[is.na(to.match$wpp)] <- to.match$median[is.na(to.match$wpp)] # no shift for years that don't match
        to.match$shift <- to.match$wpp - to.match$median
        if(sum(to.match$shift) != 0)
            pred$median.shift[[as.character(cntry)]] <- to.match$shift
    }
    if(verbose) cat("\n")
    return(pred)
}

#' @details Function \code{mig.shift.prediction.to.wpp} shifts the projected medians or means 
#'     (if \code{stat} is \dQuote{mean}), so that they correspond to the values found in the \code{migproj1dt} or \code{migproj5dt}
#'     datasets of the \pkg{wpp} package that either corresponds to the package used for the simulation itself 
#'     or is given by the \code{wpp.year} argument. Currently, the function only works for \pkg{wpp2024}.
#'     Note that regardless if it is an adjustment of the median or mean, the corresponding offset is always 
#'     converted to a shift of the median.
#' @export
#' @rdname mig.adjust
mig.shift.prediction.to.wpp <- function(sim.dir, ...){
    pred <- get.mig.prediction(sim.dir)
    new.pred <- .do.mig.shift.prediction.to.wpp(pred, ...)
    store.bayesMig.prediction(new.pred)
    invisible(new.pred)
}
