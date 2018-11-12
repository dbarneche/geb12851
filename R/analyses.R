####################
# GENERAL STAN SPECS
####################
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

######################
# AUXILLIARY FUNCTIONS
######################
readFile  <-  function (filepath, ...) {
    read.csv(filepath, header = TRUE, stringsAsFactors = FALSE, na.strings = c('', 'NA'), strip.white = TRUE, ...)
}

lenUn  <-  function (x) {
    length(unique(x))
}

################
# DATA CRUNCHING
################
removeSmallestSizeClass  <-  function (data) {
    data[data$size_class != 1, ]
}

keepQuantClassData  <-  function (data) {
    data[data$type == 'quant_class', ]
}

############################
# COVERAGE-BASED RAREFACTION
############################
incidence_freq  <-  function (data) {
    c(length(unique(data$transect_id)), unname(sort(tapply(data$transect_id, data$species, lenUn), decreasing = TRUE)))
}

outputList  <-  function (placeList, placeCov) {
    list(list = placeList,
         cov  = placeCov,
         rich = data.frame(.id = as.character(placeCov[['site']]), Coverage = placeCov[['SC']], Asym = placeCov[['q = 0']], stringsAsFactors = FALSE))
}

extractCoverageRarefactionAllDataList  <-  function (dat, dataTypeFunction) {
    siteList   <-  plyr::dlply(dat, .(sites), get(dataTypeFunction))
    locList    <-  plyr::dlply(dat, .(locality), get(dataTypeFunction))
    regList    <-  plyr::dlply(dat, .(region), get(dataTypeFunction))
    
    sitesCov  <-  iNEXT::estimateD(siteList, dataTypeFunction, base = 'coverage', level = 0.83, conf = NULL)
    locCov    <-  iNEXT::estimateD(locList,  dataTypeFunction, base = 'coverage', level = 0.89, conf = NULL)
    regCov    <-  iNEXT::estimateD(regList,  dataTypeFunction, base = 'coverage', level = 0.98, conf = NULL)
        
    list(sites    = outputList(siteList, sitesCov),
         locality = outputList(locList, locCov),
         region   = outputList(regList, regCov)
    )
}

extractCoverageRarefactionAllDataListSizeClass  <-  function (data, dataTypeFunction = 'incidence_freq') {
    sizeClasses  <-  sort(unique(data$size_class))
    sizeList     <-  vector(mode = 'list', length = length(sizeClasses))
    for (i in seq_along(sizeClasses)) {
        dat      <-  data[data$size_class == sizeClasses[i], ]
        scales   <-  c('sites', 'locality', 'region')
        covLev   <-  c(0.65, 0.85, 0.95)
        outList  <-  vector(mode = 'list', length = length(scales))
        for (j in seq_along(scales)) {
            places     <-  unique(dat[[scales[j]]])
            placeList  <-  plyr::dlply(dat, c(scales[j]), get(dataTypeFunction))
            placeCov   <-  data.frame()
            for (k in seq_along(placeList)) {
                covRar   <-  try(iNEXT::estimateD(placeList[k], dataTypeFunction, base = 'coverage', level = covLev[j], conf = NULL), silent = TRUE)
                if (!inherits(covRar, 'try-error')) {
                    covRar$site    <-  as.character(covRar$site)
                    covRar$method  <-  as.character(covRar$method)
                    placeCov       <-  rbind(placeCov, covRar, stringsAsFactors = FALSE)
                }
            }
            placeCov      <-  placeCov[placeCov[['q = 0']] >= 1, ]
            placeList     <-  placeList[names(placeList) %in% placeCov$site]
            outList[[j]]  <-  outputList(placeList, placeCov)
        }
        names(outList)  <-  scales
        sizeList[[i]]   <-  outList
    }
    names(sizeList)  <-  sizeClasses
    sizeList
}

######################
# CALCULATE MODAL SIZE
######################
modalSizeQuantClass  <-  function (data) {
    spAvSizes  <-  plyr::ddply(data, .(species), function (x) {
        data.frame(species = unique(x$species), av_size = mean(unlist(mapply(function(x, y)rep(x, y), x = x$size_cm, y = x$abun))), stringsAsFactors = FALSE)
    })
    modeSizes     <-  unique(spAvSizes[, c('species', 'av_size')])
    sizeDens      <-  density(modeSizes$av_size, from = min(modeSizes$av_size), to = max(modeSizes$av_size))
    modeSize      <-  sizeDens$x[which.max(sizeDens$y)]
    data.frame(avModeSize = modeSize, stringsAsFactors = FALSE)
}

modalSizeAllMaxSize  <-  function (data) {
    modeSizes     <-  unique(data[, c('species', 'max_size')])
    sizeDens      <-  density(modeSizes$max_size, from = min(modeSizes$max_size), to = max(modeSizes$max_size))
    modeSize      <-  sizeDens$x[which.max(sizeDens$y)]
    data.frame(avModeSize = modeSize, stringsAsFactors = FALSE)
}

####################
# STATISTICAL MODELS
# IN STAN
####################
stanFitSites  <-  function (coverageRarefactionDataList, data, humanPop, geogrVar, avSst, sizeFun) {
    scaleData     <-  coverageRarefactionDataList$sites$rich
    thereIsAZero  <-  scaleData$Coverage == 0
    
    if (any(thereIsAZero)) {
        message('The following site(s) were removed due to 0 coverage:\n', paste0(scaleData$.id[thereIsAZero], collapse = '\n'), '\n')
        data       <-  data[!(data$sites %in% scaleData$.id[thereIsAZero]), ]
        scaleData  <-  scaleData[!thereIsAZero, ]
    }
    
    # response
    scaleData$lnAsym            <-  log(scaleData$Asym)
    # fixed-effect covariates
    sizeFun                     <-  get(sizeFun)
    avModeSize                  <-  plyr::ddply(data, .(sites), sizeFun)
    scaleData$lnAvModeSize      <-  log(avModeSize$avModeSize[match(scaleData$.id, avModeSize$sites)])
    scaleData$lnGravNpop        <-  log(humanPop$grav_Npop[match(scaleData$.id, humanPop$sites)])
    scaleData$lnGravMark        <-  log(humanPop$grav_Nmarket[match(scaleData$.id, humanPop$sites)])
    scaleData$lnReefArea        <-  log(geogrVar$reefArea_km2[match(scaleData$.id, geogrVar$sites)] + 1)
    scaleData$lnIsolationReef   <-  log(geogrVar$distClosesReef_km[match(scaleData$.id, geogrVar$sites)])
    scaleData$lnIsolationCoast  <-  log(geogrVar$distCoast_km[match(scaleData$.id, geogrVar$sites)])
    scaleData$avSst             <-  avSst$mean[match(scaleData$.id, avSst$sites)]
    # random effects
    scaleData$locality          <-  data$locality[match(scaleData$.id, data$sites)]
    scaleData$region            <-  data$region[match(scaleData$.id, data$sites)]
    scaleData$id                <-  as.numeric(as.factor(scaleData$locality))
    
    set.seed(1)
    stanFit  <-  brms::brm(lnAsym ~ lnAvModeSize + lnGravNpop + lnGravMark + lnReefArea + lnIsolationReef + lnIsolationCoast + avSst + (1 | region / id), data = scaleData, family = gaussian(), prior = c(prior(normal(0, 10), 'b'), prior(normal(0, 50), 'Intercept'), prior(student_t(3, 0, 20), 'sd'), prior(student_t(3, 0, 20), 'sigma')), sample_prior = TRUE, chains = 3, cores = 3, iter = 2000, warmup = 1000)
    list(data = scaleData, modelFit = stanFit)
}

moransTest  <-  function (data, stanFitList) {
    # Moran's I test of spatial autocorrelation
    coordinates       <-  unique(data[, c('sites', 'lon', 'lat')])
    scaleData         <-  stanFitList$data
    scaleData$lon     <-  coordinates$lon[match(scaleData$.id, coordinates$sites)]
    scaleData$lat     <-  coordinates$lat[match(scaleData$.id, coordinates$sites)]
    distMat           <-  as.matrix(dist(cbind(scaleData$lon, scaleData$lat)))
    invDistMat        <-  1 / distMat
    diag(invDistMat)  <-  0 
    set.seed(10)
    ape::Moran.I(residuals(stanFitList$modelFit)[, 'Estimate'], invDistMat)
}

stanFitLocalities  <-  function (coverageRarefactionDataList, data, humanPop, geogrVar, avSst, sizeFun) {
    
    scaleData     <-  coverageRarefactionDataList$locality$rich
    thereIsAZero  <-  scaleData$Coverage == 0
    
    if (any(thereIsAZero)) {
        message('The following site(s) were removed due to 0 coverage:\n', paste0(scaleData$.id[thereIsAZero], collapse = '\n'), '\n')
        data       <-  data[!(data$locality %in% scaleData$.id[thereIsAZero]), ]
        scaleData  <-  scaleData[!thereIsAZero, ]
    }

    # response
    scaleData$lnAsym        <-  log(scaleData$Asym)
    # fixed-effect covariates
    sizeFun                 <-  get(sizeFun)
    avModeSize              <-  plyr::ddply(data, .(locality), sizeFun)
    scaleData$lnAvModeSize  <-  log(avModeSize$avModeSize[match(scaleData$.id, avModeSize$locality)])
    humanPop$locality       <-  data$locality[match(humanPop$sites, data$sites)]
    geogrVar$locality       <-  data$locality[match(geogrVar$sites, data$sites)]
    avSst$locality          <-  data$locality[match(avSst$sites, data$sites)]
    
    humanPop  <-  humanPop[!is.na(humanPop$locality), ]
    geogrVar  <-  geogrVar[!is.na(geogrVar$locality), ]
    avSst     <-  avSst[!is.na(avSst$locality), ]
    
    scaleData$lnGravNpop  <-  scaleData$lnGravMark  <-  scaleData$lnReefArea  <-  scaleData$lnIsolationReef  <-  scaleData$lnIsolationCoast  <-  scaleData$avSst  <-  NA
    for (i in 1:nrow(scaleData)) {
        scaleData$lnGravNpop[i]        <-  mean(log(humanPop$grav_Npop[humanPop$locality == scaleData$.id[i]]))
        scaleData$lnGravMark[i]        <-  mean(log(humanPop$grav_Nmarket[humanPop$locality == scaleData$.id[i]]))
        scaleData$lnReefArea[i]        <-  log(sum(geogrVar$reefArea_km2[geogrVar$locality == scaleData$.id[i]]) + 1)
        scaleData$lnIsolationReef[i]   <-  mean(log(geogrVar$distClosesReef_km[geogrVar$locality == scaleData$.id[i]]))
        scaleData$lnIsolationCoast[i]  <-  mean(log(geogrVar$distCoast_km[geogrVar$locality == scaleData$.id[i]]))
        scaleData$avSst[i]             <-  mean(avSst$mean[avSst$locality == scaleData$.id[i]])
    }
    # random effects
    scaleData$region  <-  data$region[match(scaleData$.id, data$locality)]
    scaleData$id      <-  as.numeric(as.factor(scaleData$region))

    set.seed(1)
    stanFit     <-  brms::brm(lnAsym ~ lnAvModeSize + lnGravNpop + lnGravMark + lnReefArea + lnIsolationReef + lnIsolationCoast + avSst + (1 | id), data = scaleData, family = gaussian(), prior = c(prior(normal(0, 10), 'b'), prior(normal(0, 50), 'Intercept'), prior(student_t(3, 0, 20), 'sd'), prior(student_t(3, 0, 20), 'sigma')), sample_prior = TRUE, chains = 3, cores = 3, iter = 2000, warmup = 1000)
    list(data = scaleData, modelFit = stanFit)
}

stanFitRegions <-  function (coverageRarefactionDataList, data, geogrVar, avSst, sizeFun) {
    scaleData     <-  coverageRarefactionDataList$region$rich
    thereIsAZero  <-  scaleData$Coverage == 0
    
    if (any(thereIsAZero)) {
        message('The following site(s) were removed due to 0 coverage:\n', paste0(scaleData$.id[thereIsAZero], collapse = '\n'), '\n')
        data       <-  data[!(data$region %in% scaleData$.id[thereIsAZero]), ]
        scaleData  <-  scaleData[!thereIsAZero, ]
    }

    # response
    scaleData$lnAsym        <-  log(scaleData$Asym)
    # fixed-effect covariates
    sizeFun                 <-  get(sizeFun)
    avModeSize              <-  plyr::ddply(data, .(region), sizeFun)
    scaleData$lnAvModeSize  <-  log(avModeSize$avModeSize[match(scaleData$.id, avModeSize$region)])

    geogrVar$region       <-  data$region[match(geogrVar$sites, data$sites)]
    avSst$region          <-  data$region[match(avSst$sites, data$sites)]
    
    geogrVar  <-  geogrVar[!is.na(geogrVar$region), ]
    avSst     <-  avSst[!is.na(avSst$region), ]
    scaleData$lnReefArea  <-  scaleData$avSst  <-  NA
    for (i in 1:nrow(scaleData)) {
        scaleData$lnReefArea[i]        <-  log(sum(geogrVar$reefArea_km2[geogrVar$region == scaleData$.id[i]]) + 1)
        scaleData$avSst[i]             <-  mean(avSst$mean[avSst$region == scaleData$.id[i]])
    }
    
    set.seed(1)
    stanFit     <-  brms::brm(lnAsym ~ lnAvModeSize, data = scaleData, family = gaussian(), prior = c(prior(normal(0, 10), 'b'), prior(normal(0, 50), 'Intercept'), prior(student_t(3, 0, 20), 'sigma')), sample_prior = TRUE, chains = 3, cores = 3, iter = 2000, warmup = 1000)
    stanFit2     <-  brms::brm(lnAsym ~ lnReefArea, data = scaleData, family = gaussian(), prior = c(prior(normal(0, 10), 'b'), prior(normal(0, 50), 'Intercept'), prior(student_t(3, 0, 20), 'sigma')), sample_prior = TRUE, chains = 3, cores = 3, iter = 2000, warmup = 1000)
    stanFit3     <-  brms::brm(lnAsym ~ avSst, data = scaleData, family = gaussian(), prior = c(prior(normal(0, 10), 'b'), prior(normal(0, 50), 'Intercept'), prior(student_t(3, 0, 20), 'sigma')), sample_prior = TRUE, chains = 3, cores = 3, iter = 2000, warmup = 1000)
    list(data = scaleData, modelFit = stanFit, modelFitArea = stanFit2, modelFitSst = stanFit3)
}

stanFitSitesWithSampArea  <-  function (coverageRarefactionDataList, data, humanPop, geogrVar, avSst, sizeFun) {
    scaleData     <-  coverageRarefactionDataList$sites$rich
    thereIsAZero  <-  scaleData$Coverage == 0
    
    if (any(thereIsAZero)) {
        message('The following site(s) were removed due to 0 coverage:\n', paste0(scaleData$.id[thereIsAZero], collapse = '\n'), '\n')
        data       <-  data[!(data$sites %in% scaleData$.id[thereIsAZero]), ]
        scaleData  <-  scaleData[!thereIsAZero, ]
    }
    
    # response
    scaleData$lnAsym            <-  log(scaleData$Asym)
    # fixed-effect covariates
    effortArea                  <-  tapply((data$area * data$transect_id), data$sites, sum)
    sizeFun                     <-  get(sizeFun)
    avModeSize                  <-  plyr::ddply(data, .(sites), sizeFun)
    scaleData$lnEffortArea      <-  log(effortArea[match(scaleData$.id, names(effortArea))])
    scaleData$lnAvModeSize      <-  log(avModeSize$avModeSize[match(scaleData$.id, avModeSize$sites)])
    scaleData$lnGravNpop        <-  log(humanPop$grav_Npop[match(scaleData$.id, humanPop$sites)])
    scaleData$lnGravMark        <-  log(humanPop$grav_Nmarket[match(scaleData$.id, humanPop$sites)])
    scaleData$lnReefArea        <-  log(geogrVar$reefArea_km2[match(scaleData$.id, geogrVar$sites)] + 1)
    scaleData$lnIsolationReef   <-  log(geogrVar$distClosesReef_km[match(scaleData$.id, geogrVar$sites)])
    scaleData$lnIsolationCoast  <-  log(geogrVar$distCoast_km[match(scaleData$.id, geogrVar$sites)])
    scaleData$avSst             <-  avSst$mean[match(scaleData$.id, avSst$sites)]
    # random effects
    scaleData$locality          <-  data$locality[match(scaleData$.id, data$sites)]
    scaleData$region            <-  data$region[match(scaleData$.id, data$sites)]
    scaleData$id                <-  as.numeric(as.factor(scaleData$locality))
    
    set.seed(1)
    stanFit  <-  brms::brm(lnAsym ~ lnEffortArea + lnAvModeSize + lnGravNpop + lnGravMark + lnReefArea + lnIsolationReef + lnIsolationCoast + avSst + (1 | region / id), data = scaleData, family = gaussian(), prior = c(prior(normal(0, 10), 'b'), prior(normal(0, 50), 'Intercept'), prior(student_t(3, 0, 20), 'sd'), prior(student_t(3, 0, 20), 'sigma')), sample_prior = TRUE, chains = 3, cores = 3, iter = 2000, warmup = 1000)
    list(data = scaleData, modelFit = stanFit)
}

stanFitSitesPerOcean  <-  function (coverageRarefactionDataList, data, humanPop, geogrVar, avSst, sizeFun, ocean) {
    scaleData     <-  coverageRarefactionDataList$sites$rich
    if (ocean == 'atlantic') {
        scaleData     <-  scaleData[scaleData$.id %in% unique(data$sites[data$region %in% c('caribbean', 'sw_atlantic', 'offshore_sw_atlantic', 'eastern_atlantic', 'offshore_tep', 'continental_tep')]), ]
    } else if (ocean == 'pacific') {
        scaleData     <-  scaleData[scaleData$.id %in% unique(data$sites[data$region %in% c('central_iwp', 'central_pacific', 'easter', 'polynesian', 'hawaiian', 'sw_pacific', 'western_indian', 'nw_indian')]), ]
    }
    thereIsAZero  <-  scaleData$Coverage == 0
    
    if (any(thereIsAZero)) {
        message('The following site(s) were removed due to 0 coverage:\n', paste0(scaleData$.id[thereIsAZero], collapse = '\n'), '\n')
        data       <-  data[!(data$sites %in% scaleData$.id[thereIsAZero]), ]
        scaleData  <-  scaleData[!thereIsAZero, ]
    }
    
    # response
    scaleData$lnAsym            <-  log(scaleData$Asym)
    # fixed-effect covariates
    sizeFun                     <-  get(sizeFun)
    avModeSize                  <-  plyr::ddply(data, .(sites), sizeFun)
    scaleData$lnAvModeSize      <-  log(avModeSize$avModeSize[match(scaleData$.id, avModeSize$sites)])
    scaleData$lnGravNpop        <-  log(humanPop$grav_Npop[match(scaleData$.id, humanPop$sites)])
    scaleData$lnGravMark        <-  log(humanPop$grav_Nmarket[match(scaleData$.id, humanPop$sites)])
    scaleData$lnReefArea        <-  log(geogrVar$reefArea_km2[match(scaleData$.id, geogrVar$sites)] + 1)
    scaleData$lnIsolationReef   <-  log(geogrVar$distClosesReef_km[match(scaleData$.id, geogrVar$sites)])
    scaleData$lnIsolationCoast  <-  log(geogrVar$distCoast_km[match(scaleData$.id, geogrVar$sites)])
    scaleData$avSst             <-  avSst$mean[match(scaleData$.id, avSst$sites)]
    # random effects
    scaleData$locality          <-  data$locality[match(scaleData$.id, data$sites)]
    scaleData$region            <-  data$region[match(scaleData$.id, data$sites)]
    scaleData$id                <-  as.numeric(as.factor(scaleData$locality))
    
    set.seed(1)
    stanFit  <-  brms::brm(lnAsym ~ lnAvModeSize + lnReefArea + lnIsolationReef + avSst + (1 | region / id), data = scaleData, family = gaussian(), prior = c(prior(normal(0, 10), 'b'), prior(normal(0, 50), 'Intercept'), prior(student_t(3, 0, 20), 'sd'), prior(student_t(3, 0, 20), 'sigma')), sample_prior = TRUE, chains = 3, cores = 3, iter = 2000, warmup = 1000)
    list(data = scaleData, modelFit = stanFit)
}

###########################
# SAMPLE-BASED RAREFACTION
# FOR ILLUSTRATION PURPOSES
###########################
runRegionSampleBasedRarefaction  <-  function (data) {
    plyr::dlply(data, .(region), sampleBasedRarefaction)
}

sampleBasedRarefaction  <-  function (data) {
    communityTable  <-  prepareRarefactionTable(data)
    vegan::specaccum(communityTable, method = 'random', permutations = 500)
}

prepareRarefactionTable  <-  function (data) {
    tab  <-  tapply(data$abun, list(data$transect_id, data$species), sum)
    changeNAforZeros(tab)
}

changeNAforZeros  <-  function (data) {
    data[is.na(data)]  <-  0
    data
}

packingRegionRegression  <-  function (data, regionsFitAll) {
    packingData    <-  plyr::ddply(data, .(region, locality, sites, transect_id), function (x) {
        data.frame(packing = length(unique(x$species)), area = unique(x$area), stringsAsFactors = FALSE)
    })

    packingBrms   <-  brms::brm(log(packing) ~ log(area) + (1 | region / locality / sites), data = packingData, family = gaussian(), prior = c(prior(normal(0, 10), 'b'), prior(normal(0, 50), 'Intercept'), prior(student_t(3, 0, 20), 'sd'), prior(student_t(3, 0, 20), 'sigma')), sample_prior = TRUE, chains = 3, cores = 3, iter = 2000, warmup = 1000)
    packingsBrms  <-  coef(packingBrms)$region[, , 'Intercept'][, 'Estimate']
    pred          <-  regionsFitAll$data$lnAsym
    meanModel     <-  summary(lm(packingsBrms ~ pred))
    meanSlope     <-  coef(meanModel)['pred', 'Estimate']
    meanR2        <-  meanModel$r.squared

    allMcMc  <-  brms::posterior_samples(packingBrms)
    regions  <-  regionsFitAll$data$.id
    columns  <-  paste0('r_region[', regions, ',Intercept]')
    r2s      <-  numeric(length = nrow(allMcMc))
    slopes   <-  numeric(length = nrow(allMcMc))
    for (k in 1:nrow(allMcMc)) {
        resp       <-  unlist(allMcMc[k, 'b_Intercept'] + allMcMc[k, columns])
        postModel  <-  summary(lm(resp ~ pred))
        slopes[k]  <-  coef(postModel)['pred', 'Estimate']
        r2s[k]     <-  postModel$r.squared
    }
    cisSlope  <-  c('mean' = mean(slopes), 'lower95ci' = quantile(slopes, probs = 0.025), 'upper95ci' = quantile(slopes, probs = 0.975))
    cisR2  <-  c('mean' = mean(r2s), 'lower95ci' = quantile(r2s, probs = 0.025), 'upper95ci' = quantile(r2s, probs = 0.975))

    list('packingData' = packingData,
         'packingBrms' = packingBrms,
         'meanModel' = meanModel,
         'meanSlope' = meanSlope,
         'posteriorSlope' = slopes,
         '95ciSlope' = cisSlope,
         'posteriorR2' = r2s,
         '95ciR2' = cisR2)
}

###################
# KNITTING COMMANDS
# FOR RMD FILE
###################
calcDocStats  <-  function (data, figClist, sitesFitAll, localitiesFitAll, regionsFitAll, sitesFitAllNoSmall, localitiesFitAllNoSmall, regionsFitAllNoSmall) {
    dataRichness       <-  tapply(data$species, data$sites, lenUn)
    foldDiffSRichness  <-  round(max(dataRichness) / min(dataRichness))
    samplingTable      <-  plyr::ddply(data, .(region, area, type), function (x) {
        data.frame(nOfSites = length(unique(x$sites)), stringsAsFactors = FALSE)
    })        
    fitParsSites       <-  brms::fixef(sitesFitAll$modelFit)
    fitParsLocs        <-  brms::fixef(localitiesFitAll$modelFit)
    fitParsRegion      <-  brms::fixef(regionsFitAll$modelFit)
    fitParsRegionArea  <-  brms::fixef(regionsFitAll$modelFitArea)
    fitParsRegionSst   <-  brms::fixef(regionsFitAll$modelFitSst)
    rndParsSites       <-  summary(sitesFitAll$modelFit)$random$`region:id`[1, 1] # sub-provinces
    rndParsLocs        <-  summary(localitiesFitAll$modelFit)$random$id[1, 1]

    fitParsSitesNoSmall    <-  brms::fixef(sitesFitAllNoSmall$modelFit)
    fitParsLocsNoSmall     <-  brms::fixef(localitiesFitAllNoSmall$modelFit)
    fitParsRegionNoSmall   <-  brms::fixef(regionsFitAllNoSmall$modelFit)

    # range in species richness at different scales
    richnessRange  <-  list('site' = round(range(exp(sitesFitAll$data$lnAsym))),
                            'location' = round(range(exp(localitiesFitAll$data$lnAsym))),
                            'province' = round(range(exp(regionsFitAll$data$lnAsym))))

    # r2 values
    rSquared  <-  list('site' = LoLinR::rounded(brms::bayes_R2(sitesFitAll$modelFit)['R2', ] * 100),
                       'location' = LoLinR::rounded(brms::bayes_R2(localitiesFitAll$modelFit)['R2', ] * 100),
                       'province' = list('size' = LoLinR::rounded(brms::bayes_R2(regionsFitAll$modelFit)['R2', ] * 100),
                                         'area' = LoLinR::rounded(brms::bayes_R2(regionsFitAll$modelFitArea)['R2', ] * 100),
                                         'sst'  = LoLinR::rounded(brms::bayes_R2(regionsFitAll$modelFitSst)['R2', ] * 100)
                                         )
                       )

    # fold changes at the site level
    vars           <-  c('lnAvModeSize', 'lnReefArea', 'lnIsolationReef', 'avSst')
    fcts           <-  c(exp, exp, exp, I)
    folChangeSite  <-  vector(mode = 'list', length = length(vars))
    intercept      <-  fitParsSites['Intercept', 'Estimate']
    names(folChangeSite)  <-  vars
    for (i in seq_along(vars)) {
        subData      <-  sitesFitAll$data[[vars[i]]]
        rangeOfPred  <-  range(subData)
        parEff       <-  fitParsSites[vars[i], 'Estimate']
        vals         <-  exp(rangeOfPred * parEff)
        if (parEff < 0) {
            fChange  <-  ((vals[1] - vals[2]) / vals[1]) * 100
        } else {
            fChange  <-  ((vals[2] - vals[1]) / vals[1]) * 100
        }
        folChangeSite[[i]]  <-  list(folChange = LoLinR::rounded(fChange, 1), rangeOfValues = round(fcts[[i]](rangeOfPred), 2))
    }

    # analysis between abundance and average size
    dataQuant    <-  data[data$type == 'quant_class', ]
    abuSize      <-  plyr::ddply(dataQuant, .(sites, species), function (x) {
        data.frame(region = unique(x$region),
                   lnAbun = mean(log(x$abun)),
                   lnSize =  mean(log(x$size_cm)), stringsAsFactors = FALSE)
    })
    
    abuSizeOut    <-  summary(lme4::lmer(lnAbun ~ lnSize + (1 | region), data = abuSize))
    tValAbuSize   <-  LoLinR::rounded(coef(abuSizeOut)['lnSize', 't value'], 2)
    slopeAbuSize  <-  LoLinR::rounded(coef(abuSizeOut)['lnSize', 'Estimate'], 2)
    
    specRegions  <-  plyr::ddply(data, .(species), function (x) {
        data.frame(nOfRegions = length(unique(x$region)),
                   maxSize = unique(x$max_size), stringsAsFactors = FALSE)
    })
    modelSpecRegions  <-  cor.test(specRegions$nOfRegions, specRegions$maxSize)
    corVal            <-  rounded(modelSpecRegions$estimate, 2)
    corP              <-  ifelse(modelSpecRegions$p.value < 0.0001, '< 0.0001', LoLinR::rounded(modelSpecRegions$p.value, 4))
    
    clistReg  <-  plyr::ddply(figClist, .(Gaspar_code), function (x) {
        data.frame(nOfRegions = length(unique(x$Realm)),
                   maxSize = unique(x$SpSize), stringsAsFactors = FALSE)
    })
    modelSpecRegionsClist  <-  cor.test(clistReg$nOfRegions, clistReg$maxSize)
    corValClist            <-  rounded(modelSpecRegionsClist$estimate, 2)
    corPClist              <-  ifelse(modelSpecRegionsClist$p.value < 0.0001, '< 0.0001', LoLinR::rounded(modelSpecRegionsClist$p.value, 4))
    
    list(dataRichness = dataRichness, foldDiffSRichness = foldDiffSRichness, samplingTable = samplingTable, fitParsSites = fitParsSites, fitParsLocs = fitParsLocs, fitParsRegion = fitParsRegion, fitParsRegionArea = fitParsRegionArea, fitParsRegionSst = fitParsRegionSst, rndParsSites = rndParsSites, rndParsLocs = rndParsLocs, fitParsSitesNoSmall = fitParsSitesNoSmall, fitParsLocsNoSmall = fitParsLocsNoSmall, fitParsRegionNoSmall = fitParsRegionNoSmall, richnessRange = richnessRange, rSquared = rSquared, folChangeSite = folChangeSite, dataQuant = dataQuant, abuSize = abuSize, abuSizeOut = abuSizeOut, tValAbuSize = tValAbuSize, slopeAbuSize = slopeAbuSize, specRegions = specRegions, modelSpecRegions = modelSpecRegions, corVal = corVal, corP = corP, clistReg = clistReg, modelSpecRegionsClist = modelSpecRegionsClist, corValClist = corValClist, corPClist = corPClist)
}
