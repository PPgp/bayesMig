get.wpp.mig.data <- function(start.year = 1950, present.year = 2020, 
                            wpp.year = 2019, my.mig.file = NULL, 
                            annual = FALSE, exclude.from.world = NULL, 
                            ignore.last.observed = FALSE, 
                            use.wpp.data = TRUE, verbose = FALSE) {

    ########################################
    # set data and match with areas
    ########################################
    if(!is.null(my.mig.file)){
        migdata <- bayesTFR:::do.read.subnat.file(my.mig.file, present.year = present.year)
        if(! "code" %in% colnames(migdata) && ! "country_code" %in% colnames(migdata))
           stop("Columns country_code or code must be present in the data file.")
         if("code" %in% colnames(migdata)) colnames(migdata)[colnames(migdata) == "code"] <- "country_code" # rename "code" to "country_code"
         if("country" %in% colnames(migdata)) colnames(migdata)[colnames(migdata) == "country"] <- "name" # rename country column to "name"
         if(! "name" %in% colnames(migdata)) migdata$name <- migdata$country_code
         locations <- bayesTFR:::create.sublocation.dataset(migdata)
    } else {
        migdata <- read.UNmig(wpp.year=wpp.year, 
                           present.year=present.year, annual = annual,
                           use.wpp.data = use.wpp.data, 
                           verbose=verbose)

        # get region and area data
        locations <- bayesTFR:::create.sublocation.dataset(migdata)
    }
    loc_data <- locations$loc_data
    include <- locations$include & ! (loc_data$country_code %in% exclude.from.world)
    prediction.only <- locations$prediction.only | loc_data$country_code %in% exclude.from.world
    
    data_incl <- migdata[include,]
    nr_countries_estimation <- nrow(data_incl)
    if(any(!is.na(prediction.only))) { # move prediction countries at the end of data
        data_prediction <- migdata[prediction.only,]
        data_incl <- rbind(data_incl, data_prediction)
    }

    MIGmatrix.regions <- bayesTFR:::get.observed.time.matrix.and.regions(
        data_incl, loc_data, 
        start.year = start.year, 
        present.year = present.year, annual = annual, 
        datacolnames=c(country.code='country_code', country.name='name', reg.name='reg_name',
                       reg.code='reg_code', area.name='area_name', area.code='area_code'),
        interpolate = wpp.year < 2022 && annual && is.null(my.mig.file),
        ignore.last.observed = ignore.last.observed)
    
    return(list(mig.matrix = MIGmatrix.regions$obs_matrix, 
                mig.matrix.all = MIGmatrix.regions$obs_matrix_all, 
                regions = MIGmatrix.regions$regions, 
                nr.countries.estimation = nr_countries_estimation
                )
           )
}

read.UNmig <- function(wpp.year, annual = FALSE, ...) {
    migration <- bayesTFR:::do.read.un.file("migration", wpp.year, annual = annual, ...)$data
    pop <- bayesTFR:::do.read.un.file("pop", wpp.year, annual = annual, ...)$data

    #List of all possible countries
    UNlocations <- bayesTFR:::load.bdem.dataset('UNlocations', wpp.year=wpp.year)
    fullCountryCodeVec <- UNlocations$country_code[UNlocations$location_type==4]
    fullCountryNameVec <- UNlocations$name[UNlocations$location_type==4]
    
    #Figure out the countries of overlap
    fullDataIndices=(fullCountryCodeVec %in% migration$country_code & fullCountryCodeVec %in% pop$country_code)
    fullCountryCodeVec=fullCountryCodeVec[fullDataIndices]
    fullCountryNameVec=as.character(fullCountryNameVec[fullDataIndices])

    nC <- length(fullCountryCodeVec)
    
    #Construct a matrix of initial populations
    exclude.columns <- c("country_code", "name", "country", "last.observed", "include_code")
    initialPopMat <- merge(data.frame(country_code=fullCountryCodeVec), pop, sort=FALSE)
    numcols <- as.numeric(setdiff(colnames(initialPopMat), exclude.columns))
    initialPopMat <- initialPopMat[,-which(colnames(pop) %in% c(exclude.columns, as.character(numcols[numcols <=1950])))] # need end-period pop
    rownames(initialPopMat) <- fullCountryCodeVec
    
    #Construct a matrix of total migration counts
    migCountMat <- merge(data.frame(country_code=fullCountryCodeVec), migration, sort=FALSE)
    migCountMat <- migCountMat[,-which(colnames(migCountMat) %in% exclude.columns)]
    if(!annual || (wpp.year < 2022 && annual))
        migCountMat <- migCountMat[, substr(colnames(migCountMat), 6, 9) %in% colnames(initialPopMat)]
    else migCountMat <- migCountMat[, colnames(migCountMat) %in% colnames(initialPopMat)]
    rownames(migCountMat) <- fullCountryCodeVec
    
    #Convert migration counts and initial populations to a matrix of migration "rates"
    # as count/(end pop - mig)
    migdata <- as.matrix(migCountMat/(initialPopMat - migCountMat))
    migdata <- cbind(data.frame(country_code=fullCountryCodeVec, name = fullCountryNameVec),
                     migdata)
    return(migdata)
}