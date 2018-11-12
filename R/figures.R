######################
# AUXILLIARY FUNCTIONS
######################
extrafont::loadfonts(quiet = TRUE)

toDev  <-  function (expr, dev, filename, ..., verbose = TRUE) {
    if (verbose) {
        cat(sprintf('Creating %s\n', filename))
    }
    dev(filename, family = 'CM Roman', ...)
    on.exit(dev.off())
    eval.parent(substitute(expr))
}

toPdf <- function (expr, filename, ...) {
    toDev(expr, pdf, filename, ...)
}

toRadians  <-  function (grades) {
    (grades / 180 * pi)
} 

niceLimits  <-  function (range_x, gap = 0.05) {
    c(range_x[1] - (range_x[1] * gap), range_x[2] + (range_x[2] * gap))
}

linearRescale   <-  function (x, rOut) {
    p  <-  (x - min(x)) / (max(x) - min(x))
    rOut[[1]] + p * (rOut[[2]] - rOut[[1]])
}

firstElement  <-  function (x) {
    x[1]
}

drawCircleLoglog  <-  function (x, y, radius, nv = 100, border = NULL, col = NA, lty = 1, lwd = 1) {
    # from https://stackoverflow.com/questions/15921359/how-to-draw-a-circle-in-a-log-log-plot-in-r
    xylim     <-  par('usr')
    plotdim   <-  par('pin')
    ymult     <-  (xylim[4] - xylim[3]) / (xylim[2] - xylim[1]) * plotdim[1] / plotdim[2]
    angleInc  <-  2 * pi / nv
    angles    <-  seq(0, 2 * pi - angleInc, by = angleInc)
    if (length(col) < length(radius)) {
        col  <-  rep(col, length.out = length(radius))
    }
    for (circle in 1:length(radius)) {
        xv  <-  exp(cos(angles) * log(radius[circle])) * x[circle]
        yv  <-  exp(sin(angles) * ymult * log(radius[circle])) * y[circle]
        polygon(xv, yv, border = border, col = col[circle], lty = lty, lwd = lwd)
    }
    invisible(list(x = xv, y = yv))
}

######################
# AUXILLIARY FUNCTIONS
######################
getTransectsAtRandom  <-  function (data, nOfTransects, iter = 1000) {
    set.seed(1)
    chosenTrs  <-  sample(unique(data$transect_id), nOfTransects * iter, replace = TRUE)
    split(chosenTrs, rep(seq_len(iter), each = nOfTransects))
}

avBodySize  <-  function (data, ...) {
    chosenTrs  <-  getTransectsAtRandom(data, ...)
    modalSize  <-  plyr::ldply(chosenTrs, function (x, data) {
        data   <-  data[data$transect_id %in% x, ]
        onlyOneSpecies  <-  lenUn(data$species) == 1
        if (onlyOneSpecies) {
            data.frame(avModeSize = unique(data$max_size), stringsAsFactors = FALSE)
        } else {
            modalSizeAllMaxSize(data)
        }
    }, data = data)
    data.frame(avMaxSize = mean(modalSize$avModeSize), varMaxSize = var(modalSize$avModeSize), stringsAsFactors = FALSE)
}

avMaxSizeAllTrs  <-  function (data) {
    allTrs  <-  plyr::ddply(data, .(transect_id), function (x) {
        onlyOneSpecies  <-  lenUn(x$species) == 1
        if (onlyOneSpecies) {
            data.frame(avModeSize = unique(x$max_size), stringsAsFactors = FALSE)
        } else {
            modalSizeAllMaxSize(x)
        }
    })
    data.frame(avMaxSize = mean(allTrs$avModeSize), varMaxSize = var(allTrs$avModeSize), stringsAsFactors = FALSE)
}

brmsLikeScatterPPC  <-  function (object) {
    object    <-  restructure(object)
    newdArgs  <-  brms:::nlist(object, NULL, NA, TRUE, list, check_response = TRUE)
    
    newdArgs$newdata        <-  NULL
    newdArgs$internal       <-  TRUE
    newdArgs$only_response  <-  TRUE
    
    sdata     <-  do.call(brms::standata, newdArgs)
    predArgs  <-  brms:::nlist(object, NULL, NULL, FALSE, NULL, 'uncertainty', list, TRUE, NULL, NULL, NULL, 5, sort = FALSE, summary = FALSE)
    list(observed  = as.vector(sdata[[paste0('Y', brms:::usc(NULL))]]), 
         predicted = do.call('predict', predArgs))
}

individualPPcheckPlot  <-  function (brmsFitList, xlab = NA, ylab = NA, main, matchCol, colorScheme) {
    ppcArgs     <-  brmsLikeScatterPPC(brmsFitList$modelFit)
    data        <-  brmsFitList$data
    data$color  <-  colorScheme$cpoint[match(data[[matchCol]], colorScheme$region)]
    data$shape  <-  colorScheme$ppoint[match(data[[matchCol]], colorScheme$region)]
    x           <-  colMeans(ppcArgs$predicted)
    y           <-  ppcArgs$observed
    plot(x, y, xlab = xlab, ylab = ylab, axes = FALSE, type = 'n', xpd = NA)
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col = 'grey90', border = NA)
    LoLinR::whiteGrid()
    box()
    axis(1)
    axis(2, las = 1)
    LoLinR::proportionalLabel(0.5, 1.1, main, xpd = NA, cex = 1.4, font = 3, adj = c(0.5, 0.5))
    points(x, y, pch = data$shape, bg = data$color, cex = 1.1)
    abline(0, 1, lty = 2, lwd = 1.5)
}

amplitude  <-  function (vec) {
    max(vec) - min(vec)
}

individualPlotFig2  <-  function (brmsFitList, modelName = 'modelFit', xlab = NA, ylab = NA, main, matchCol, variable, r2 = FALSE, r2pos = 'right', covLeg = FALSE, leg, colorScheme, ...) {
    modelFit    <-  brmsFitList[[modelName]]
    data        <-  brmsFitList$data
    data$color  <-  colorScheme$cpoint[match(data[[matchCol]], colorScheme$region)]
    data$shape  <-  colorScheme$ppoint[match(data[[matchCol]], colorScheme$region)]

    fixefs  <-  brms::fixef(modelFit)
    ranefs  <-  try(brms::ranef(modelFit), silent = TRUE)
    rnfCor  <-  rep(0, nrow(data))
    if (!inherits(ranefs, 'try-error')) {
        for(k in seq_along(ranefs)) {
            matchRnd  <-  strsplit(names(ranefs)[k], ':')[[1]]
            if (length(matchRnd) == 2) {
                matchDataVec  <-  paste0(data[[matchRnd[1]]], '_', data[[matchRnd[2]]])
            } else {
                matchDataVec  <-  data[[matchRnd]]
            }
            rnfCor  <-  rnfCor + ranefs[[k]][match(matchDataVec, rownames(ranefs[[k]]))]
        }
    }
    
    fefCor  <-  rep(0, nrow(data))
    for (j in setdiff(rownames(fixefs), c('Intercept', variable))) {
        fefCor  <-  fefCor - fixefs[j, 'Estimate'] * data[[j]] + fixefs[j, 'Estimate'] * mean(data[[j]])
    }
        
    x   <-  data[[variable]]
    y   <-  data$lnAsym - rnfCor + fefCor
    plot(x, y, xlab = xlab, ylab = ylab, ylim = c(1.8, 7.2), axes = FALSE, type = 'n', xpd = NA, yaxs = 'i', ...)
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col = 'grey90', border = NA)
    LoLinR::whiteGrid()
    box()
    axis(1)
    axis(2, las = 1)
    LoLinR::proportionalLabel(0.5, 1.1, main, xpd = NA, cex = 1.4, font = 3, adj = c(0.5, 0.5))
    points(x, y, pch = data$shape, bg = data$color, cex = 1.4, lwd = 0.5)
    avFef  <-  rep(0, length(range(x)))
    for (j in setdiff(rownames(fixefs), c('Intercept', variable))) {
        avFef  <-  avFef + fixefs[j, 'Estimate'] * mean(data[[j]])
    }
    lines(range(x), fixefs['Intercept', 'Estimate'] + fixefs[variable, 'Estimate'] * range(x) + avFef, lty = 2, lwd = 1.8)
    
    mcmcSamp    <-  brms::posterior_samples(modelFit, pars = paste0('b_', rownames(fixefs)), exact_match = TRUE)
    sizeVec     <-  seq(min(x), max(x), length.out = 50)
    mcmcFefCor  <-  matrix(0, nrow(mcmcSamp), length(sizeVec))
    for (a in 1:nrow(mcmcFefCor)) {
        mcmcFefCor[a, ]  <-  mcmcFefCor[a, ] + mcmcSamp[a, 'b_Intercept'] + mcmcSamp[a, paste0('b_', variable)] * sizeVec
        for (j in setdiff(rownames(fixefs), c('Intercept', variable))) {
            mcmcFefCor[a, ]  <-  mcmcFefCor[a, ] + mcmcSamp[a, paste0('b_', j)] * mean(data[[j]])
        }
    }
    lines(sizeVec, apply(mcmcFefCor, 2, quantile, probs = 0.025), lty = 2, lwd = 1.3, col = 'grey30')
    lines(sizeVec, apply(mcmcFefCor, 2, quantile, probs = 0.975), lty = 2, lwd = 1.3, col = 'grey30')

    LoLinR::proportionalLabel(0.02, 0.95, leg, xpd = NA, font = 3, cex = 1.5)
    if (r2) {
        r2     <-  brms::bayes_R2(modelFit)
        if (r2pos == 'right') {
            xadj  <-  1
            xpos  <-  0.96
            ypos  <-  0.95
        } else if (r2pos == 'left') {
            xadj  <-  0
            xpos  <-  0.03
            ypos  <-  0.85
        }
        LoLinR::proportionalLabel(xpos, ypos, substitute(italic('Bayesian R'^2) == a, list(a = LoLinR::rounded(r2[1, 'Estimate'], 2))), adj = c(xadj, 0.5), cex = 1.25)
        LoLinR::proportionalLabel(xpos, ypos - 0.08, substitute('95% C.I.: '* a - b, list(a = LoLinR::rounded(r2[1, 'Q2.5'], 2), b = LoLinR::rounded(r2[1, 'Q97.5'], 2))), adj = c(xadj, 0.5), cex = 1.25)
    }
    if (covLeg) {
        LoLinR::proportionalLabel(0.03, 0.1, substitute('Cov.: '* a %+-% b, list(a = LoLinR::rounded(mean(data$Coverage), 2), b = LoLinR::rounded(sd(data$Coverage), 2))), adj = c(0, 0.5), cex = 1.25)
    }
}

######################
# TEXT FIGURES
######################
makeFigure1  <-  function (dest, ...) {
    toPdf(fig1(...), dest, width = 7, height = 9)
    extrafont::embed_fonts(dest)
}

fig1  <-  function (figClist, regionsInfo, coverageRarefactionAllDataList, data, colorScheme) {
    regionsInfo[, c('cpoint', 'ppoint')]  <-  colorScheme[match(regionsInfo$gasparProvince, colorScheme$region), c('cpoint', 'ppoint')]
    
    sites     <-  coverageRarefactionAllDataList$sites$rich
    locality  <-  coverageRarefactionAllDataList$locality$rich
    region    <-  coverageRarefactionAllDataList$region$rich
    
    coord  <-  unique(data[, c('sites', 'locality', 'region', 'lon', 'lat')])
    coord  <-  coord[coord$sites %in% sites$.id, c('sites', 'locality', 'region', 'lat', 'lon')]
    coord[, c('color', 'shape')]  <-  colorScheme[match(coord$region, colorScheme$region), c('cpoint', 'ppoint')]
    coord  <-  coord[match(sites$.id, coord$sites), ]
    
    lon          <-  mapproj::mapproject(coord$lon, coord$lat, projection = 'gilbert', orientation = c(90, 0, 200))$x
    lat          <-  mapproj::mapproject(coord$lon, coord$lat, projection = 'gilbert', orientation = c(90, 0, 200))$y
    cpoint       <-  'black'
    cpointT      <-  coord$color
    ppoint       <-  coord$shape
    
    layout(matrix(c(rep(1, 10), rep(c(2, 3), 3), rep(c(4, 5), 3)), 11, 2, byrow = TRUE))
    par(mai = c(0.4, 0, 0.2, 1.2), omi = c(0.6, 0.2, 0, 0.2), cex = 1)
    plot(NA, xlim = c(-1.0733908, 1.0782939), ylim = c(-0.6092952, 0.6058435), xlab = '', ylab = '', axes = FALSE, xaxs = 'i', yaxs = 'i')
    curveLeftX   <-  cos(seq(toRadians(160.85), toRadians(206.6), length.out = 50))
    curveLeftY   <-  sin(seq(toRadians(160.85), toRadians(206.6), length.out = 50))
    curveRightX  <-  cos(seq(toRadians(19.15), toRadians(-26.6), length.out = 50))
    curveRightY  <-  sin(seq(toRadians(19.15), toRadians(-26.6), length.out = 50))
    xpol  <-  c(curveLeftX, rev(curveRightX), curveLeftX[1])
    ypol  <-  c(curveLeftY, rev(curveRightY), curveRightY[1])
    polygon(xpol, ypol, col = 'grey90', border = 'black')
    maps::map('world', projection = 'gilbert', col = 'white', ylim = c(-50, 35), orientation = c(90, 0, 200), bg = NA, lwd = 0.5, fill = TRUE, resolution = 0, wrap = TRUE, border = 'white', add = TRUE)
    maps::map('world', projection = 'gilbert', col = 'black', ylim = c(-50, 35), orientation = c(90, 0, 200), lwd = 0.2, boundary = TRUE, interior = FALSE, fill = FALSE, add = TRUE, resolution = 0, wrap = TRUE)
    par(lwd = 0.7)
    points(lon, lat, pch = ppoint, col = cpoint, bg = cpointT, cex = 1.1)
    polygon(c(-1, 1, 1, -1), c(0.33, 0.33, 0.6, 0.6), col = 'white', border = FALSE)
    polygon(c(-1, 1, 1, -1), c(-0.45, -0.45, -0.55, -0.55), col = 'white', border = FALSE)
    polygon(xpol, ypol, col = NA, border = 'black')    
    LoLinR::proportionalLabel(0.07, 0.76, '(a)', xpd = NA, font = 3)
    
    orderedRegions  <-  c('nw_indian', 'western_indian', 'central_iwp', 'sw_pacific', 'central_pacific', 'polynesian', 'hawaiian', 'easter', 'offshore_tep', 'continental_tep', 'caribbean', 'sw_atlantic', 'offshore_sw_atlantic', 'eastern_atlantic')

    for (z in seq_along(orderedRegions)) {
        regionLab  <-  strsplit(orderedRegions[z], '_')[[1]]
        if (length(regionLab) > 1) {
            for (j in seq_along(regionLab)) {
                chars  <-  nchar(regionLab[j])
                regionLab[j]  <-  ifelse(chars <= 3, toupper(regionLab[j]), Hmisc::capitalize(regionLab[j]))
            }
            regionLab  <-  paste0(regionLab, collapse = ' ')
        } else {
            regionLab  <-  Hmisc::capitalize(regionLab)
        }
        LoLinR::proportionalLabel(1, 0.78 - (z - 1) * 0.05, text = FALSE, pch = unique(coord$shape[coord$region == orderedRegions[z]]), bg = unique(coord$color[coord$region == orderedRegions[z]]), cex = 0.9, xpd = NA) 
        LoLinR::proportionalLabel(1.02, 0.78 - (z - 1) * 0.05, regionLab, adj = c(0, 0.5), xpd = NA, cex = 0.7)
    }

    # splines to separate major regions
    lines(curveLeftX + 0.13, curveLeftY, lty = 2)
    sx1 <- c(-0.3925947, -0.2092355, -0.2101806, -0.2192736, -0.2306370, -0.2462992, -0.2660296, -0.2868434, -0.3020100, -0.3047068)
    sy1 <- c(0.3283223, -0.2519616, -0.2887141, -0.3295579, -0.3582546, -0.3781924, -0.4025558, -0.4232428, -0.4408411, -0.4448748)
    lines(sx1, sy1, lty = 2)
    sx2 <- c(0.4337201, 0.6952173, 0.7309902, 0.7501213, 0.7440247, 0.7254352, 0.7005993, 0.5939791)
    sy2 <- c(0.32741187, 0.02175180, -0.04913709, -0.11772102, -0.18249025, -0.22758669, -0.26900673, -0.44886236)
    lines(sx2 - 0.25, sy2, lty = 2)
    lines(sx2, sy2, lty = 2)

    # and now labels for regions
    LoLinR::proportionalLabel(0.09, 0.07, 'Western\nIndian', xpd = NA, font = 3, adj = c(0.5, 0.5))
    LoLinR::proportionalLabel(0.22, 0.85, 'Central\nIndo-Pacific', xpd = NA, font = 3, adj = c(0.5, 0.5))
    LoLinR::proportionalLabel(0.52, 0.07, 'Eastern\nIndo-Pacific', xpd = NA, font = 3, adj = c(0.5, 0.5))
    LoLinR::proportionalLabel(0.65, 0.85, 'Eastern\nPacific', xpd = NA, font = 3, adj = c(0.5, 0.5))
    LoLinR::proportionalLabel(0.85, 0.07, 'Atlantic', xpd = NA, font = 3, adj = c(0.5, 0.5))

    par(mai = c(0.5, 0.9, 0, 0.2), cex = 1, mgp = c(3, 0.5, 0), tcl = -0.3)
    firstElement  <-  function (x) {x[1]}
    regRich   <-  tapply(figClist$Gaspar_code, figClist$Site_name, lenUn)
    slon2     <-  tapply(figClist$Longitude, figClist$Site_name, firstElement)
    slat2     <-  tapply(figClist$Latitude, figClist$Site_name, firstElement)
    slon      <-  mapproj::mapproject(slon2, slat2, projection = 'gilbert', orientation = c(90, 0, 200))$x
    crPoint   <-  'black'
    crPointT  <-  regionsInfo[match(names(regRich), regionsInfo$site_name), 'cpoint']
    srPoint   <-  regionsInfo[match(names(regRich), regionsInfo$site_name), 'ppoint']
    xli       <-  mapproj::mapproject(c(30, 140, -159, -100, 7), rep(0, 5), projection = 'gilbert', orientation = c(90, 0, 200))$x
    plot(NA, xlab = '', ylab = '', xlim = range(xli), ylim = c(0, 2000), xaxt = 'n', las = 1, cex.axis = 0.8, xpd = NA)
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col = 'grey90', border = NA)
    LoLinR::whiteGrid()
    box()
    axis(1, at = xli, labels = c(30, 140, -159, -100, 7), cex.axis = 0.8)
    points(slon, regRich, pch = srPoint, col = 'black', bg = crPointT, cex = 1.1)
    LoLinR::proportionalLabel(0.02, 0.95, '(b)', xpd = NA, cex = 1.1, font = 3)
    LoLinR::proportionalLabel(-0.24, 0.5, 'Species richness', adj = c(0.5, 0.5), cex = 1.2, xpd = NA, srt = 90)
    LoLinR::proportionalLabel(0.95, 0.9, 'Checklists', cex = 0.9, font = 3, adj = c(1, 0.5))
    LoLinR::proportionalLabel(0.95, 0.8, substitute(italic(n) == a, list(a = length(regRich))), cex = 0.9, font = 3, adj = c(1, 0.5))
    
    lonReg   <-  tapply(lon, coord$region, firstElement)
    lonReg2  <-  tapply(coord$lon, coord$region, firstElement)
    ppReg    <-  tapply(ppoint, coord$region, firstElement)
    cpReg    <-  tapply(cpointT, coord$region, firstElement)
    
    plot(NA, xlab = '', ylab = '', xlim = range(xli), ylim = c(7,900), xaxt = 'n', las = 1,  cex.axis = 0.8, xpd = NA)
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col = 'grey90', border = NA)
    LoLinR::whiteGrid()
    box()
    axis(1, at = xli, labels = c(30, 140, -159, -100, 7), cex.axis = 0.8)
    points(lonReg, region$Asym, pch = ppReg, col = cpoint, bg = cpReg, cex = 1.1)
    LoLinR::proportionalLabel(0.95, 0.8, substitute('Cov.: ' * a %+-% b, list(a = LoLinR::rounded(mean(region$Coverage), 2), b = LoLinR::rounded(sd(region$Coverage), 2))), adj = c(1, 0.5), font = 3, cex = 0.9)
    LoLinR::proportionalLabel(0.02, 0.95, '(c)', xpd = NA, cex = 1.1, font = 3)
    LoLinR::proportionalLabel(-0.24, 0.5, 'Species richness', adj = c(0.5, 0.5), cex = 1.2, xpd = NA, srt = 90)
    LoLinR::proportionalLabel(0.95, 0.9, 'Province scale', cex = 0.9, font = 3, adj = c(1, 0.5))
    LoLinR::proportionalLabel(0.95, 0.7, substitute(italic(n) == a, list(a = length(lonReg))), cex = 0.9, font = 3, adj = c(1, 0.5))

    lonLoc   <-  tapply(lon, coord$locality, firstElement)
    lonLoc2  <-  tapply(coord$lon, coord$locality, firstElement)
    ppLoc    <-  tapply(ppoint, coord$locality, firstElement)
    cpLoc    <-  tapply(cpointT, coord$locality, firstElement)
    
    plot(NA, xlab = '', ylab = '', xlim = range(xli), ylim = c(0, 400), xaxt = 'n', las = 1,  cex.axis = 0.8, xpd = NA)
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col = 'grey90', border = NA)
    LoLinR::whiteGrid()
    box()
    axis(1, at = xli, labels = c(30, 140, -159, -100, 7), cex.axis = 0.8)
    points(lonLoc, locality$Asym, pch = ppLoc, col = cpoint, bg = cpLoc, cex = 1.1)
    LoLinR::proportionalLabel(0.95, 0.8, substitute('Cov.: ' * a %+-% b, list(a = LoLinR::rounded(mean(locality$Coverage), 2), b = LoLinR::rounded(sd(locality$Coverage), 2))), adj = c(1, 0.5), font = 3, cex = 0.9)
    LoLinR::proportionalLabel(0.02, 0.95, '(d)', xpd = NA, cex = 1.1, font = 3)
    LoLinR::proportionalLabel(-0.24, 0.5, 'Species richness', adj = c(0.5, 0.5), cex = 1.2, xpd = NA, srt = 90)
    LoLinR::proportionalLabel(0.5, -0.26, 'Longitude', adj = c(0.5, 0.5), cex = 1.2, xpd = NA)
    LoLinR::proportionalLabel(0.95, 0.9, 'Sub-province scale', cex = 0.9, font = 3, adj = c(1, 0.5))
    LoLinR::proportionalLabel(0.95, 0.7, substitute(italic(n) == a, list(a = length(lonLoc))), cex = 0.9, font = 3, adj = c(1, 0.5))

    plot(NA, xlab = '', ylab = '', xlim = range(xli), ylim = c(7, 300), xaxt = 'n', las = 1,  cex.axis = 0.8, xpd = NA)
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col = 'grey90', border = NA)
    LoLinR::whiteGrid()
    box()
    axis(1, at = xli, labels = c(30, 140, -159, -100, 7), cex.axis = 0.8)
    points(lon, sites$Asym, pch = ppoint, col = cpoint, bg = cpointT,cex = 1.1)
    LoLinR::proportionalLabel(0.95, 0.8, substitute('Cov.: ' * a %+-% b, list(a = LoLinR::rounded(mean(sites$Coverage), 2), b = LoLinR::rounded(sd(sites$Coverage), 2))), adj = c(1, 0.5), font = 3, cex = 0.9)
    LoLinR::proportionalLabel(0.02, 0.95, '(e)', xpd = NA, cex = 1.1, font = 3)
    LoLinR::proportionalLabel(-0.24, 0.5, 'Species richness', adj = c(0.5, 0.5), cex = 1.2, xpd = NA, srt = 90)
    LoLinR::proportionalLabel(0.5, -0.26, 'Longitude', adj = c(0.5, 0.5), cex = 1.2, xpd = NA)
    LoLinR::proportionalLabel(0.95, 0.9, 'Site scale', cex = 0.9, font = 3, adj = c(1, 0.5))
    LoLinR::proportionalLabel(0.95, 0.7, substitute(italic(n) == a, list(a = length(lon))), cex = 0.9, font = 3, adj = c(1, 0.5))
}

makeFigure2  <-  function (dest, ...) {
    toPdf(fig2(...), dest, width = 14, height = 14)
    extrafont::embed_fonts(dest)
}

fig2  <-  function (sitesFitAll, localitiesFitAll, regionsFitAll, colorScheme) {
    
    par(mfrow = c(3, 3), mai = c(1.02, 0.62, 0.32, 0.3),  omi = c(0, 0.5, 0.5, 0), cex = 1, cex.lab = 1.4, cex.axis = 1.2, mgp = c(3, 0.5, 0), tck = -0.03)
    individualPlotFig2(brmsFitList = sitesFitAll, xlab = 'ln Maximum body size mode (cm)', ylab = 'ln Species richness', main = 'Site', matchCol = 'region', variable = 'lnAvModeSize', r2 = TRUE, covLeg = TRUE, leg = '(a)', colorScheme = colorScheme, xlim = c(2.2, 4.25))
    individualPlotFig2(brmsFitList = localitiesFitAll, xlab = 'ln Maximum body size mode (cm)', main = 'Sub-province', matchCol = 'region', variable = 'lnAvModeSize', r2 = TRUE, covLeg = TRUE, leg = '(b)', colorScheme = colorScheme, xlim = c(2.2, 4.25))
    individualPlotFig2(brmsFitList = regionsFitAll, xlab = 'ln Maximum body size mode (cm)', main = 'Province', matchCol = '.id', variable = 'lnAvModeSize', r2 = TRUE, r2pos = 'right', covLeg = TRUE, leg = '(c)', colorScheme = colorScheme, xlim = c(2.2, 4.25))

    individualPlotFig2(brmsFitList = sitesFitAll, xlab = substitute('ln (Reef Area + 1) (km'^2 * ')'), ylab = 'ln Species richness', main = NA, matchCol = 'region', variable = 'lnReefArea', leg = '(d)', colorScheme = colorScheme, xlim = c(-0.2, 8.2))
    individualPlotFig2(brmsFitList = localitiesFitAll, xlab = substitute('ln (Reef Area + 1) (km'^2 * ')'), main = NA, matchCol = 'region', variable = 'lnReefArea', leg = '(e)', colorScheme = colorScheme, xlim = c(-0.2, 8.2))
    individualPlotFig2(brmsFitList = regionsFitAll, xlab = substitute('ln (Reef Area + 1) (km'^2 * ')'), modelName = 'modelFitArea', main = NA, matchCol = '.id', variable = 'lnReefArea', r2 = TRUE, r2pos = 'left', leg = '(f)', colorScheme = colorScheme, xlim = c(-0.2, 8.2))

    individualPlotFig2(brmsFitList = sitesFitAll, xlab = substitute('SST ('*degree*'C)'), ylab = 'ln Species richness', main = NA, matchCol = 'region', variable = 'avSst', leg = '(g)', colorScheme = colorScheme, xlim = c(21, 30))
    individualPlotFig2(brmsFitList = localitiesFitAll, xlab = substitute('SST ('*degree*'C)'), main = NA, matchCol = 'region', variable = 'avSst', leg = '(h)', colorScheme = colorScheme, xlim = c(21, 30))
    individualPlotFig2(brmsFitList = regionsFitAll, xlab = substitute('SST ('*degree*'C)'), modelName = 'modelFitSst', main = NA, matchCol = '.id', variable = 'avSst', r2 = TRUE, r2pos = 'left', leg = '(i)', colorScheme = colorScheme, xlim = c(21, 30))
}

makeFigure3  <-  function (dest, ...) {
    toPdf(fig3(...), dest, width = 6, height = 8.8)
    extrafont::embed_fonts(dest)
}

fig3  <-  function (data, coverageRarefactionAllDataListSize, scale, colorScheme) {
    ######################################
    # AVERAGE ASYMPTOTE AT THE LOCAL SCALE
    ######################################
    coord    <-  unique(data[, c('sites', 'locality', 'region', 'lon', 'lat')])
    classes  <-  names(coverageRarefactionAllDataListSize)
    sizes    <-  paste0(c('0-7', '7-15', '15-30', '30-50', '50-80', '> 80'), ' cm')
    coord[, c('color', 'shape')]  <-  colorScheme[match(coord$region, colorScheme$region), c('cpoint', 'ppoint')]
    
    par(mar = c(0.6, 0.6, 0.2, 0.3), mfrow = c(3, 2), omi = c(1.8, 0.9, 0.6, 0.1), cex = 1, cex.axis = 1)
    
    for (k in seq_along(classes)) {
        dat          <-  coverageRarefactionAllDataListSize[[k]][[scale]]$rich
        dat          <-  dat[dat$Coverage != 0, ]
        subCoord     <-  coord[coord[[scale]] %in% dat$.id, c(scale, 'color', 'shape')]
        subCoord     <-  subCoord[match(dat$.id, subCoord[[scale]]), ]
        cpoint       <-  'black'
        cpointT      <-  subCoord$color
        ppoint       <-  subCoord$shape
        regions      <-  data$region[match(tolower(subCoord[[scale]]), data[[scale]])]
        dat$regions  <-  regions
        
        if (scale != 'region') {
            ano         <-  try(aov(dat$Asym ~ dat$regions))
            isTryError  <-  any(class(ano) %in% 'try-error')
            if (!isTryError) {
                ano  <-  TukeyHSD(ano)
            }
        }
        
        if (scale == 'sites') {
            ylim  <-  c(-5, 52)
        } else if (scale == 'locality') {
            ylim  <-  c(-10, 100)
        } else {
            ylim  <-  c(0, 250)
        }
        plot(0, 0, type = 'n', xlab = '', xaxt = 'n', yaxt = 'n', cex = 1.5, ylim = ylim, xlim = c(1, lenUn(data$region)), ylab = '')
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col = 'grey90', border = NA)
        LoLinR::whiteGrid()
        polygon(usr[c(1, 2, 2, 1, 1)], c(usr[c(3, 3)], 0, 0, usr[3]), col = 'white', border = NA)
        lines(c(2.5, 2.5), usr[3:4], lty = 2, col = 'grey30', lwd = 0.8)
        lines(c(4.5, 4.5), usr[3:4], lty = 2, col = 'grey30', lwd = 0.8)
        lines(c(8.5, 8.5), usr[3:4], lty = 2, col = 'grey30', lwd = 0.8)
        lines(c(10.5, 10.5), usr[3:4], lty = 2, col = 'grey30', lwd = 0.8)
        box()

        orderedRegions  <-  c('nw_indian', 'western_indian', 'central_iwp', 'sw_pacific', 'central_pacific', 'polynesian', 'hawaiian', 'easter', 'offshore_tep', 'continental_tep', 'caribbean', 'sw_atlantic', 'offshore_sw_atlantic', 'eastern_atlantic')

        set.seed(1)
        for (f in seq_along(orderedRegions)) {
            tab  <-  dat[dat$regions == orderedRegions[f], ]
            if (nrow(tab) == 0) {
                next
            }
            cp   <-  cpointT[match(orderedRegions[f], regions)]
            pp   <-  ppoint[match(orderedRegions[f], regions)]
            
            if (length(tab$Asym) > 1) {
                den    <-  density(tab$Asym, from = min(tab$Asym), to = max(tab$Asym))
                y1     <-  linearRescale(den$y, c(f, f - 0.25))
                y2     <-  linearRescale(den$y, c(f, f + 0.25))
                polygon(c(y1, rev(y2)), c(den$x, rev(den$x)), col = cp, border = 'black')
                points(f, mean(c(0, usr[3])), pch = pp, col = 'black', bg = cp)
            } else {
                points(jitter(rep(f, nrow(tab)), amount = 0.1), tab$Asym, pch = pp, col = 'black', bg = cp, cex = 1.2)
                points(f, mean(c(0, usr[3])), pch = pp, col = 'black', bg = cp)
            }
            if (scale != 'region') {
                if (!isTryError) {
                    regionStats  <-  ano[['dat$regions']][grep(orderedRegions[f], rownames(ano[['dat$regions']])), ]
                    tukeyReg     <-  regionStats[, 'p adj'] <= 0.05
                    if (all(tukeyReg)) {
                        text(f, max(tab$Asym) + 0.05 * (usr[4] - usr[3]), '*', cex = 1.4, adj = c(0.5, 0.5), font = 2)
                    }
                }
            }
        }
        LoLinR::proportionalLabel(0.95, 0.9, sizes[k], adj = c(1, 0.5), font = 3, cex = 1)
        LoLinR::proportionalLabel(0.95, 0.8, substitute(a %+-% b, list(a = LoLinR::rounded(mean(dat$Coverage), 2), b = LoLinR::rounded(sd(dat$Coverage), 2))), adj = c(1, 0.5), font = 3, cex = 1)
        
        if (k %in% c(1, 3, 5)) {
            axis(2, las = 1)
        }
        if (k %in% 5:6) {
            text(1.5, -0.13 * (usr[4] - usr[3]), 'Western Indian', xpd = NA, cex = 0.8, adj = c(1, 1), font = 3, srt = 45)
            text(3.5, -0.13 * (usr[4] - usr[3]), 'Central Indo-Pacific', xpd = NA, cex = 0.8, adj = c(1, 1), font = 3, srt = 45)
            text(6.5, -0.13 * (usr[4] - usr[3]), 'Eastern Indo-Pacific', xpd = NA, cex = 0.8, adj = c(1, 1), font = 3, srt = 45)
            text(9.5, -0.13 * (usr[4] - usr[3]), 'Eastern Pacific', xpd = NA, cex = 0.8, adj = c(1, 1), font = 3, srt = 45)
            text(12.5, -0.13 * (usr[4] - usr[3]), 'Atlantic', xpd = NA, cex = 0.8, adj = c(1, 1), font = 3, srt = 45)
        }
        LoLinR::proportionalLabel(0.07, 0.93, paste0('(', letters[k], ')'), font = 3, adj = c(0.5,0.5))
        
    }
    mtext('Species richness', side = 2, line = 2, outer = TRUE, cex = 1.4)
}

makeFigure4  <-  function (dest, ...) {
    toPdf(fig4(...), dest, width = 9, height = 7.5)
    extrafont::embed_fonts(dest)
}

fig4  <-  function (data, regionSampleBasedRarefaction, figClist, posteriorPacking, regionsFitAll, colorScheme) {
    coord    <-  unique(data[, c('sites', 'locality', 'region', 'lon', 'lat')])
    coord[, c('color', 'shape')]  <-  colorScheme[match(coord$region, colorScheme$region), c('cpoint', 'ppoint')]

    cpoint        <-  tapply(coord$color, coord$region, unique)
    ppoint        <-  tapply(coord$shape, coord$region, unique)
    names(figClist)[names(figClist) %in% c('Gaspar_code', 'SpSize')]  <-  c('species', 'max_size')
    maxSizeClist  <-  plyr::ddply(figClist, .(Realm), function (x) {
        tabSizes  <-  plyr::ddply(x, .(Site_code), function (k) {
            modalSizeAllMaxSize(k)
        })
        data.frame(avMaxSize = mean(tabSizes$avModeSize), varMaxSize = var(tabSizes$avModeSize))
    })
    maxSizeStats  <-  plyr::ddply(data, .(region), function (data, nOfTrs) {
        plyr::ldply(nOfTrs, function (k, data) {
            avBodySize(data, k)
        }, data = data)
    }, nOfTrs = c(1, 5, 10, 15))
    maxSizeStats  <-  rbind(maxSizeStats, plyr::ddply(data, .(region), avMaxSizeAllTrs))
       
    #############
    # RAREFACTION
    #############
    layout(matrix(c(rep(c(1, 1, 2, 2), 2), rep(c(1, 1, 3, 3), 2)), 4, 4))
    par(omi = rep(0.5, 4), mar = c(3, 4.5, 0.5, 0), cex = 1, cex.axis = 0.8, cex.lab = 1.1)
    plot(0, 0, type = 'n', las = 1, xlab = '', ylab = '', xlim = c(0, 3100), yaxs = 'i', ylim = c(0, 1200), axes = FALSE)
    LoLinR::proportionalLabel(-0.2, 0.5, 'Species richness', adj = c(0.5, 0.5), xpd = NA, srt = 90, cex = 1.1)
    LoLinR::proportionalLabel(0.5, -0.23, 'Number of bootstrapped transects', adj = c(0.5, 0.5), xpd = NA, cex = 1.1)
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col = 'grey90', border = NA)
    LoLinR::whiteGrid()
    box()    
    axis(1, mgp = c(3, 0.5, 0))
    axis(2, las = 1, mgp = c(3, 0.7, 0))
    LoLinR::proportionalLabel(0.05, 0.91, '(a)', font = 3, adj = c(0.5, 0.5))
    plyr::l_ply(seq_along(regionSampleBasedRarefaction), function (x, rfct, cols) {
        plot(rfct[[x]], ci.type = 'polygon', ci.col = cols[x], add = TRUE, lwd = 0.3)
    }, rfct = regionSampleBasedRarefaction, cols = cpoint)
    
    #################
    # BODY SIZE STATS
    #################
    par(mar = c(3, 4.5, 1.5, 0))
    plot(0, 0, type = 'n', xlim = c(0.5, 7.5), ylim = c(10, 30), las = 1, xlab = '', ylab = '', xpd = NA, axes = FALSE)
    LoLinR::proportionalLabel(-0.2, 0.5, 'Mean [max. body size mode (cm)]', adj = c(0.5, 0.5), xpd = NA, srt = 90, cex = 1.1)
    LoLinR::proportionalLabel(0.5, -0.23, 'Number of bootstrapped transects', adj = c(0.5, 0.5), xpd = NA, cex = 1.1)
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col = 'grey90', border = NA)
    LoLinR::whiteGrid()
    axis(1, at = c(1, 2, 3, 4, 5.5, 7), c(1, 5, 10, 15, 'All', 'Checklist'), mgp = c(3, 0.5, 0))
    axis(2, las = 1, mgp = c(3, 0.7, 0))
    plyr::l_ply(seq_len(lenUn(maxSizeStats$region)), function (x, data, cols, sha) {
        points(c(1:4, 5.5), data$avMaxSize[data$region == unique(data$region)[x]], pch = sha[x], bg = cols[x], cex = 1, lwd = 0.7, type = 'b', lty = 1)
    }, data = maxSizeStats, cols = cpoint, sha = ppoint)
    segments(seq(4.25, 5.25, length = 7), usr[3], seq(4.25, 5.25, length = 7), usr[4], lwd = 3.5, col = 'grey90')
    plyr::l_ply(seq_len(lenUn(maxSizeClist$Realm)), function (x, data, cols, sha) {
        points(7, data$avMaxSize[data$Realm == unique(data$Realm)[x]], pch = sha[x], bg = cols[x], cex = 1)
    }, data = maxSizeClist, cols = cpoint, sha = ppoint)
    LoLinR::proportionalLabel(0.05, 0.91, '(b)', font = 3, adj = c(0.5, 0.5))
    box()
    
    plot(0, 0, type = 'n', xlim = c(0.5, 7.5), ylim = c(0, 15), las = 1, xlab = '', ylab = '', xpd = NA, axes = FALSE)
    LoLinR::proportionalLabel(0.5, -0.23, 'Transects', adj = c(0.5, 0.5), xpd = NA, cex = 1.1)
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col = 'grey90', border = NA)
    LoLinR::whiteGrid()
    axis(1, at = c(1, 2, 3, 4, 5.5, 7), c(1, 5, 10, 15, 'All', 'Checklist'), mgp = c(3, 0.5, 0))
    axis(2, las = 1, mgp = c(3, 0.7, 0))
    LoLinR::proportionalLabel(-0.17, 0.5, substitute(paste(sigma~'[max. body size mode (cm)]')), xpd = NA, adj = c(0.5, 0.5), xpd = NA, srt = 90, cex = 1.1)
    plyr::l_ply(seq_len(lenUn(maxSizeStats$region)), function (x, data, cols, sha) {
        points(c(1:4,5.5), sqrt(data$varMaxSize[data$region == unique(data$region)[x]]), pch = sha[x], bg = cols[x], cex = 1, lwd = 0.7, type = 'b', lty = 1)
    }, data = maxSizeStats, cols = cpoint, sha = ppoint)
    segments(seq(4.2, 5.3, length = 9), usr[3], seq(4.2, 5.3, length = 9), usr[4], lwd = 3.5, col = 'grey90')
    plyr::l_ply(seq_len(lenUn(maxSizeClist$Realm)), function (x, data, cols, sha) {
        points(7, sqrt(data$varMaxSize[data$Realm == unique(data$Realm)[x]]), pch = sha[x], bg = cols[x], cex = 1)
    }, data = maxSizeClist, cols = cpoint, sha = ppoint)
    LoLinR::proportionalLabel(0.05, 0.91, '(c)', font = 3, adj = c(0.5, 0.5))
    box()
}

makeFigure5  <-  function (dest, ...) {
    toPdf(fig5(...), dest, width = 8, height = 4.8)
    extrafont::embed_fonts(dest)
}

fig5  <-  function (data, brmsFit, colorScheme) {
    coord    <-  unique(data[, c('sites', 'locality', 'region', 'lon', 'lat')])
    coord[, c('color', 'shape')]  <-  colorScheme[match(coord$region, colorScheme$region), c('cpoint', 'ppoint')]
    bg            <-  tapply(coord$color, coord$region, unique)
    pch           <-  tapply(coord$shape, coord$region, unique)
    
    layout(matrix(c(rep(1, 4), rep(2:3, 2)), 2, 4))
    # Plot 1
    par(omi = rep(0.6, 4))
    par(mai = c(0.6732, 0.7412, 0.4, 0.0772))
    richness   <-  brmsFit$data$Asym
    maxSizeTo  <-  exp(brmsFit$data$lnAvModeSize)
    plot(NA, type = 'n', xlab = 'Species richness', ylab = 'Maximum body size mode (cm)', cex.lab = 1.5, ylim = c(10, 26), xlim = c(40, 900), axes = FALSE, log = 'xy', xpd = NA)
    usr  <-  par('usr')
    rect(10^usr[1], 10^usr[3], 10^usr[2], 10^usr[4], col = 'grey90', border = NA)
    LoLinR::whiteGrid(log = 'xy')
    axis(1, cex.axis = 1.2)
    axis(2, cex.axis = 1.2, las = 1)
    box()
    model  <-  summary(lm(log(maxSizeTo) ~ log(richness)))
    LoLinR::proportionalLabel(0.02, 0.98, '(a)', font = 3, cex = 1.2, log = 'xy')
    drawCircleLoglog(x = 65, y = 20, radius = 1.7, border = 'dodgerblue2', col = LoLinR::transparentColor('dodgerblue2', 0.3), lty = 1)
    drawCircleLoglog(x = 565, y = 12, radius = 1.7, border = 'tomato', col = LoLinR::transparentColor('tomato', 0.3), lty = 1)
    points(richness, maxSizeTo, pch = pch, col = 'black', bg = bg, cex = 2)
    lines(range(richness), exp(coef(model)[1] + coef(model)[2] * range(log(richness))), lty = 2)

    # Plots 4 and 5
    set.seed(3)
    par(mai = c(0.3, 0.3412, 0.4, 0.0772))
    x   <- seq(0, 70, length = 100)
    y   <- dnorm(x, mean = 35, sd = 15)
    plot(x, y*1.1, type = 'n', xlim = c(0, 100), ylab = '', xlab = '', axes = FALSE)
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col = 'grey90', border = NA)
    LoLinR::whiteGrid()
    box(bty = 'l')
    
    for (i in seq(0, 30, length = 5)) {
        polygon(c(x + i + rnorm(1, 0, 1)), c(y * rnorm(1, 0.4, 0.25)), col = LoLinR::transparentColor('dodgerblue2', 0.3), border = 'dodgerblue4')
    }
    LoLinR::proportionalLabel(0.02, 0.78, 'Low richness', adj = c(0, 0.5), cex = 1.2, font = 2)
    LoLinR::proportionalLabel(0.02, 0.63, 'Province = 5 spp.', adj = c(0, 0.5), cex = 1)
    LoLinR::proportionalLabel(0.02, 0.98, '(b)', font = 3, cex = 1.2, xpd = NA)
    LoLinR::proportionalLabel(0.98, 0.9, 'Larger size', adj = c(1, 0.5), cex = 1)
    LoLinR::proportionalLabel(0.98, 0.8, 'Larger range', adj = c(1, 0.5), cex = 1)
    LoLinR::proportionalLabel(c(0.65, 0.65), c(0, 0.5), text = FALSE, type = 'l', lty = 2)
    LoLinR::proportionalLabel(0.65, 0.52, 'Local = 4 spp.', cex = 1, adj = c(0.5, 0))

    par(mai = c(0.6732, 0.3412, 0.0268, 0.0772))
    x   <- seq(0, 20, length = 100)
    y   <- dnorm(x, mean = 10, sd = 4)
    plot(x, y * 0.8, type = 'n', xlim = c(0, 100), ylim = c(0.02, 0.12), axes = FALSE, xlab = 'Geographic range', cex.lab = 1.5, ylab = '')
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col = 'grey90', border = NA)
    LoLinR::whiteGrid()
    box(bty = 'l')
    
    for (i in seq(0, 80, length = 20)) {
        polygon(c(x + i + rnorm(1, 0, 1)), c(y * rnorm(1, 0.5, 0.25)), col = LoLinR::transparentColor('tomato', 0.3), border = 'tomato3')
    }
    LoLinR::proportionalLabel(0.02, 0.88, 'High richness', adj = c(0, 0.5), cex = 1.2, font = 2)
    LoLinR::proportionalLabel(0.02, 0.73, 'Province = 20 spp.', adj = c(0, 0.5), cex = 1)
    LoLinR::proportionalLabel(0.98, 0.9, 'Smaller size', adj = c(1, 0.5), cex = 1)
    LoLinR::proportionalLabel(0.98, 0.8, 'Smaller range', adj = c(1, 0.5), cex = 1)
    LoLinR::proportionalLabel(c(0.58, 0.58), c(0, 0.7), text = FALSE, type = 'l', lty = 2)
    LoLinR::proportionalLabel(0.58, 0.72, 'Local = 4 spp.', cex = 1, adj = c(0.5, 0))
    mtext('Abundance', 2, line = 0.9, at = 0.046, adj = -0.8, cex = 1)
    arrows(x0 = -80, y0 = 0.19, x1 = -7.27, y1 = 0.23, length = 0.08, angle = 30, xpd = NA, lwd = 1.2)
    arrows(x0 = -19.7, y0 = 0.065, x1 = -5.98, y1 = 0.067, length = 0.08, angle = 30, xpd = NA, lwd = 1.2)
}

makeFigureS1  <-  function (dest, ...) {
    toPdf(figS1(...), dest, width = 14, height = 5)
    extrafont::embed_fonts(dest)
}

figS1  <-  function (sitesFitAll, localitiesFitAll, regionsFitAll, ...) {
    par(mfrow = c(1, 3), mai = c(1.02, 0.62, 0.32, 0.3),  omi = c(0, 0.5, 0.5, 0), cex = 1, cex.lab = 1.4, cex.axis = 1.2, mgp = c(3, 0.5, 0), tck = -0.03)
    individualPPcheckPlot(brmsFitList = sitesFitAll, ylab = 'ln Observed richness', main = 'Site', matchCol = 'region', ...)
    individualPPcheckPlot(brmsFitList = localitiesFitAll, xlab = 'ln Predicted richness (Mean posterior)', main = 'Sub-province', matchCol = 'region', ...)
    individualPPcheckPlot(brmsFitList = regionsFitAll, main = 'Province', matchCol = '.id', ...)
}

makeFigureS2  <-  function (dest, ...) {
    toPdf(figS2(...), dest, width = 7, height = 7)
    extrafont::embed_fonts(dest)
}

figS2  <-  function (data, posteriorPacking, colorScheme) {
    coord           <-  unique(data[, c('sites', 'locality', 'region', 'lon', 'lat')])
    lon             <-  mapproj::mapproject(coord$lon, coord$lat, projection = 'gilbert', orientation = c(90, 0, 200))$x
    lonReg          <-  tapply(lon, coord$region, firstElement)
    orderedRegions  <-  names(lonReg)[order(lonReg)]
    packingBrms     <-  posteriorPacking$packingBrms
    coord[, c('color', 'shape')]  <-  colorScheme[match(coord$region, colorScheme$region), c('cpoint', 'ppoint')]

    allMcMc    <-  brms::posterior_samples(packingBrms)
    columns    <-  paste0('r_region[', orderedRegions, ',Intercept]')
    transRich  <-  data.frame(packing = exp(unlist(allMcMc[, 'b_Intercept'] + allMcMc[, columns])), region = rep(orderedRegions, each = nrow(allMcMc)), stringsAsFactors = FALSE)
    transRich       <-  transRich[order(match(transRich$region, orderedRegions)), ]
    transRich$num   <-  as.integer(factor(transRich$region, levels = orderedRegions))
    transRich[, c('color', 'shape')]  <-  colorScheme[match(transRich$region, colorScheme$region), c('cpoint', 'ppoint')]
    
    par(omi = c(0.8, 0.2, 0, 0), xpd = NA, cex = 1, cex.lab = 1.3)
    mai  <-  par('mai')
    par(mai = c(mai[1:2], 0.2, mai[4]))
    plot(packing ~ num, data = transRich, xlab = '', ylab = '', xlim = c(0.5, length(orderedRegions) + 0.2), ylim = c(0, 12), axes = FALSE, type = 'n', xpd = NA)
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col = 'grey90', border = NA)
    LoLinR::whiteGrid()
    axis(2, las = 1, mgp = c(3, 0.7, 0))
    box()
    LoLinR::proportionalLabel(-0.13, 0.5, substitute('Posterior species packing (species m'^{-2} ~ ')'), adj = c(0.5, 0.5), xpd = NA, srt = 90, cex = 1.3)
    points(packing ~ num, data = transRich, pch = transRich$shape, col = LoLinR::transparentColor(transRich$color, 0.01), bg = LoLinR::transparentColor(transRich$color, 0.01), cex = 0.7)
    for (z in seq_along(orderedRegions)) {
        dat   <-  transRich[transRich$region == orderedRegions[z], ]
        resp  <-  dat$packing
        bd    <-  quantile(resp, probs = c(0.005, 0.995))
        resp  <-  resp[resp >= bd[1] & resp <= bd[2]]
        dens  <-  density(resp, from = min(resp), to = max(resp))
        dens$y  <-  linearRescale(dens$y, c(z - 0.05, z - 0.5))
        polygon(c(dens$y, z - 0.05), c(dens$x, min(dens$x)), col = LoLinR::transparentColor(unique(dat$color), 0.8), border = NA)
        lines(dens$y[dens$x >= 0], dens$x[dens$x >= 0], col = 'grey30', lwd = 0.5)
        text(z, resp[which.min(abs(resp - dens$x[which.min(dens$y)]))], '*', cex = 1.5, col = 'black')
        regionLab  <-  strsplit(orderedRegions[z], '_')[[1]]
        if (length(regionLab) > 1) {
            for (j in seq_along(regionLab)) {
                chars  <-  nchar(regionLab[j])
                regionLab[j]  <-  ifelse(chars <= 3, toupper(regionLab[j]), Hmisc::capitalize(regionLab[j]))
            }
            regionLab  <-  paste0(regionLab, collapse = ' ')
        } else {
            regionLab  <-  Hmisc::capitalize(regionLab)
        }
        text(z, -0.8, regionLab, adj = c(1, 0), srt = 45, xpd = NA, cex = 1.3)
    }
}

makeFigureS4  <-  function (dest, ...) {
    toPdf(figS4(...), dest, width = 9, height = 8)
    extrafont::embed_fonts(dest)
}

figS4  <-  function (data, colorScheme) {
    
    pack.quant    <-  data[data$type == 'quant_class', ]
    coord         <-  unique(pack.quant[, c('sites', 'locality', 'region', 'lon', 'lat')])
    coord[, c('color', 'shape')]  <-  colorScheme[match(coord$region, colorScheme$region), c('cpoint', 'ppoint')]
    region        <-  sort(unique(data$region))
    bg            <-  tapply(coord$color, coord$region, unique)
    pch           <-  tapply(coord$shape, coord$region, unique)
    spTrans       <-  tapply(data$transect_id, data$species, lenUn)
    spSzCl        <-  tapply(data$size_class, data$species, unique)
    
    par(omi = rep(0.5, 4), mai = c(0.8, 0.5412, 0.5412, 0.4772), mfrow = c(2, 2), cex = 1)
    plot(NA, xlim = c(0.05, 0.95), ylim = c(-0.02, 1), axes = FALSE, xlab = '', ylab = '')
    LoLinR::proportionalLabel(-0.15, 0.5, 'Frequency of occurrence', xpd = NA, cex = 1.3, srt = 90, adj = c(0.5, 0.5))
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col = 'grey90', border = NA)
    LoLinR::whiteGrid()
    arrows(x0 = 0.05, y0 = -0.13, x1 = 0.95, length = 0.08, angle = 30, xpd = NA, lwd = 1.2)
    mtext('Body size', 1, line = 1.2, cex = 1.3)
    LoLinR::proportionalLabel(0.15, 0.7, 'Local\nabundance',  adj = c(0.5, 0.5), cex = 1)
    LoLinR::proportionalLabel(0.84, 0.7, 'Geographic\nrange', adj = c(0.5, 0.5), cex = 1)
    text(0.5, 0.5, 'Abundance x Range', cex = 1, adj = c(0.5, 0.5), col = 'dodgerblue2')
    size  <-  seq(0, 1, length.out = 20)
    polygon(c(0, size, 1), c(-1, 7 * ((size)^2) * ((1 - size)^2), -1), col = NA, border = 'dodgerblue2', lwd = 1.2, lty = 1)
    lines(size, (size * 0.8)^2, lty = 2, lwd = 0.8)
    lines(size, (0.8 - size * 0.8)^2, lty = 2, lwd = 0.8)
    LoLinR::proportionalLabel(0.07, 0.93, '(a)', font = 3, cex = 1.2, adj = c(0.5, 0.5))
    box()
    
    # Plot 2 - Number of transects in which each species was observed (sorted by size)
    par(mai = c(0.6732, 0.5412, 0, 0.4772), xpd = NA)
    plot(spSzCl + rnorm(length(spSzCl), 0, 0.05), spTrans, type = 'n', las = 1, xlab = '', ylab = 'Number of transects', log = 'y', xaxt = 'n', cex.lab = 1.3)
    usr  <-  par('usr')
    rect(usr[1], 10^usr[3], usr[2], 10^usr[4], col = 'grey90', border = NA)
    LoLinR::whiteGrid(log = 'y')
    box()
    points(spSzCl + rnorm(length(spSzCl), 0, 0.05), spTrans, pch = 21, col = 'transparent', bg = rgb(0, 0, 0, 0.1), cex = 1.1)
    for (i in 1:lenUn(region)) {
        subdata   <-  data[data$region == region[i], ]
        subTrans  <-  tapply(subdata$transect_id, subdata$species, lenUn)
        subSzCl   <-  tapply(subdata$size_class, subdata$species, unique)
        points(1:6 + rnorm(1, 0, 0.05), 10^tapply(log10(subTrans), subSzCl, mean), cex = 1.1, pch = pch[i], bg = bg[i], type = 'b', lwd = 0.8)
    }
    axis(1, at = 1:6, labels = NA, mgp = c(3, 0.5, 0))
    LoLinR::proportionalLabel((1:6 - usr[1]) / (usr[2] - usr[1]), rep(-0.05, 6), c('0-7', '7-15', '15-30', '30-50', '50-80', '> 80'), log = 'y', srt = 25, xpd = NA, adj = c(1, 1))
    LoLinR::proportionalLabel(0.02, 0.98, '(b)', font = 3, cex = 1.2, log = 'y')

    # Plot 3
    set.seed(1)
    plot(pack.quant$size_class + rnorm(dim(pack.quant)[1], 0, 0.05), pack.quant$abun, type = 'n', xlim = c(0.7, 6.4), xlab = 'Body size (cm)', ylab = 'Abundance within transects', axes = FALSE, log = 'y',  las = 1, cex.lab = 1.3)
    usr  <-  par('usr')
    rect(usr[1], 10^usr[3], usr[2], 10^usr[4], col = 'grey90', border = NA)
    LoLinR::whiteGrid(log = 'y')
    box()
    axis(1, at = 1:6, labels = NA, mgp = c(3, 0.5, 0))
    for (j in 0:5) {
        axis(2, at = 10^j, labels = substitute(10^a, list(a = j)), las = 1)
    }
    LoLinR::proportionalLabel((1:6 - usr[1]) / (usr[2] - usr[1]), rep(-0.05, 6), c('0-7', '7-15', '15-30', '30-50', '50-80', '> 80'), log = 'y', srt = 25, xpd = NA, adj = c(1, 1))
    
    points(pack.quant$size_class + rnorm(dim(pack.quant)[1], 0, 0.05), pack.quant$abun, lwd = 0.8, pch = 16, col = LoLinR::transparentColor('grey50', 0.05))
    for (i in seq_along(region)) {
        subset <- pack.quant[pack.quant$region == region[i], ]
        points(1:6 + rnorm(1, 0, 0.05), 10^tapply(log10(subset$abun), subset$size_class, mean), cex = 1.1, pch = pch[i], bg = bg[i], type = 'b', lwd = 0.8)
    }
    LoLinR::proportionalLabel(0.02, 0.98, '(c)', font = 3, cex = 1.2, xpd = NA, log = 'y')
    
    # Plot 4 - Distribution
    sizePack    <-  tapply(data$size_class, data$species, function (x) unique(x))
    regionPack  <-  tapply(data$region, data$species, lenUn)
    plot(sizePack + rnorm(length(sizePack), 0, 0.05), regionPack + rnorm(length(sizePack), 0, 0.05), xlim = c(0.7, 6.4), ylim = c(0.8, lenUn(data$region) + 5), type = 'n', axes = FALSE, xlab = 'Body size (cm)', ylab = 'Regions occupied', cex.lab = 1.3, log = 'y')
    usr  <-  par('usr')
    rect(usr[1], 10^usr[3], usr[2], 10^usr[4], col = 'grey90', border = NA)
    LoLinR::whiteGrid(log = 'y')
    box()
    points(sizePack + rnorm(length(sizePack), 0, 0.05), regionPack + rnorm(length(sizePack), 0, 0.05), pch = 16, col = LoLinR::transparentColor('black', 0.05), cex = 1.2)
    points(1:6, tapply(regionPack, sizePack, max),  cex = 1.5, pch = 22, bg = 'seagreen3', type = 'b', lwd = 0.8)
    points(1:6, tapply(regionPack, sizePack, mean), cex = 1.5, pch = 22, bg = 'dodgerblue', type = 'b', lwd = 0.8)
    LoLinR::proportionalLabel(0.05, 0.85, log = 'y', text = FALSE, cex = 1.2, pch = 22, bg = 'seagreen3', adj = c(0, 0.5))
    LoLinR::proportionalLabel(0.05, 0.78, log = 'y', text = FALSE, cex = 1.2, pch = 22, bg = 'dodgerblue', adj = c(0, 0.5))
    LoLinR::proportionalLabel(0.07, 0.85, log = 'y', 'Maximum', cex = 0.9, adj = c(0, 0.5))
    LoLinR::proportionalLabel(0.07, 0.78, log = 'y', 'Mean', cex = 0.9, adj = c(0, 0.5))
    axis(1, at = 1:6, labels = NA, mgp = c(3, 0.5, 0))
    LoLinR::proportionalLabel((1:6 - usr[1]) / (usr[2] - usr[1]), rep(-0.05, 6), c('0-7', '7-15', '15-30', '30-50', '50-80', '> 80'), log = 'y', srt = 25, xpd = NA, adj = c(1, 1))
    axis(2, at = c(1, 4, 7, 10, 13), las = 1)
    LoLinR::proportionalLabel(0.02, 0.98, '(d)', font = 3, cex = 1.2, xpd = NA, log = 'y')
}
