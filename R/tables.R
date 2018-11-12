tableRounding  <-  function (x, rounding = 2) {
    format(round(x, rounding), nsmall = rounding)
}

makeTableA1  <-  function (stanFitList) {
	data     <-  stanFitList$data
	tableA1  <-  Hmisc::rcorr(as.matrix(data[, c('lnAvModeSize', 'lnGravNpop', 'lnGravMark', 'lnReefArea', 'lnIsolationReef', 'lnIsolationCoast', 'avSst')]), type = 'pearson')
	significant  <-  tableA1$P < 0.001
	tableA1  <-  matrix(as.character(round(tableA1$r, 2)), nrow(tableA1$r), ncol(tableA1$r))
	tableA1  <-  ifelse(significant, paste0('\\textbf{', tableA1, '}'), tableA1)
	tableA1[is.na(tableA1)]  <-  '-'
	tableA1
}

makeTableA23  <-  function (sitesFitModel, localitiesFitModel, regionsFitModel) {
    # sites
	siteModel  <-  summary(sitesFitModel$modelFit)
	fixSite    <-  siteModel$fixed
	sdSite     <-  do.call(rbind.data.frame, siteModel$random)
    # locality
	locModel  <-  summary(localitiesFitModel$modelFit)
	fixLoc    <-  locModel$fixed
	sdLoc     <-  do.call(rbind.data.frame, locModel$random)
    # region
	regModel  <-  summary(regionsFitModel$modelFit)
	fixReg    <-  regModel$fixed
	cbind(
		rbind(as.matrix(tableRounding(sdSite[2:1, c(1, 3, 4)], 2)), matrix(rep('', 6), 2, 3), as.matrix(tableRounding(fixSite[, c(1, 3, 4)], 2))),
		rbind(matrix(rep('', 3), 1, 3), as.matrix(tableRounding(sdLoc[1, c(1, 3, 4)], 2)), matrix(rep('', 6), 2, 3), as.matrix(tableRounding(fixLoc[, c(1, 3, 4)], 2))),
		rbind(matrix(rep('', 12), 4, 3), as.matrix(tableRounding(fixReg[, c(1, 3, 4)], 2)), matrix(rep('', 18), 6, 3))
	)
}

makeTableA4  <-  function (sitesFitModel) {
    # sites
	siteModel  <-  summary(sitesFitModel$modelFit)
	fixSite    <-  siteModel$fixed
	sdSite     <-  do.call(rbind.data.frame, siteModel$random)
	rbind(as.matrix(tableRounding(sdSite[2:1, c(1, 3, 4)], 2)), matrix(rep('', 6), 2, 3), as.matrix(tableRounding(fixSite[, c(1, 3, 4)], 2)))
}

makeTableA5  <-  function (regionsFitModel) {
    # region
	regModel1  <-  summary(regionsFitModel$modelFitArea)
	fixReg1    <-  regModel1$fixed
	regModel2  <-  summary(regionsFitModel$modelFitSst)
	fixReg2    <-  regModel2$fixed
	cbind(
		rbind(as.matrix(tableRounding(fixReg1[, c(1, 3, 4)], 2)), matrix(rep('', 3), 1, 3)),
		rbind(as.matrix(tableRounding(fixReg2[1, c(1, 3, 4), drop = FALSE], 2)), matrix(rep('', 3), 1, 3), as.matrix(tableRounding(fixReg2[2, c(1, 3, 4), drop = FALSE], 2)))
	)
}

makeTableA6  <-  function (sitesFitAtl, sitesFitPac) {
    # sites
	atlModel  <-  summary(sitesFitAtl$modelFit)
	fixAtl    <-  atlModel$fixed
	sdAtl     <-  do.call(rbind.data.frame, atlModel$random)
    # locality
	pacModel  <-  summary(sitesFitPac$modelFit)
	fixPac    <-  pacModel$fixed
	sdPac     <-  do.call(rbind.data.frame, pacModel$random)

	cbind(
		rbind(as.matrix(tableRounding(sdAtl[2:1, c(1, 3, 4)], 2)), matrix(rep('', 6), 2, 3), as.matrix(tableRounding(fixAtl[, c(1, 3, 4)], 2))),
		rbind(matrix(rep('', 3), 1, 3), as.matrix(tableRounding(sdPac[1, c(1, 3, 4)], 2)), matrix(rep('', 6), 2, 3), as.matrix(tableRounding(fixPac[, c(1, 3, 4)], 2)))
	)
}

makeTableS1  <-  function (data) {
	effort    <-  tapply(data$transect_id, paste0(data$studyName, data$sites), lenUn)
	area      <-  tapply(data$area, paste0(data$studyName, data$sites), firstElement)
	richness  <-  tapply(data$species, paste0(data$studyName, data$sites), lenUn)
	tableS1               <-  unique(data[, c('region', 'locality', 'sites', 'lon', 'lat', 'studyName')])
	uniqueSites           <-  data.frame(sites = unique(tableS1$sites), stringsAsFactors = FALSE)
	uniqueSites$newSites  <-  sapply(uniqueSites$sites, function (x) {
		j    <-  strsplit(x, '_')[[1]]
		num  <-  nchar(j[length(j)])
		if (num == 1) {
			paste0(c(j[1:(length(j) - 1)], paste0('0', j[length(j)])), collapse = '_')
		} else {
			x
		}
	})
	tableS1$newSites  <-  uniqueSites$newSites[match(tableS1$sites, uniqueSites$sites)]
	tableS1           <-  tableS1[with(tableS1, order(region, locality, newSites, studyName)), ]
	tableS1$effort    <-  effort[match(paste0(tableS1$studyName, tableS1$sites), names(effort))]
	tableS1$area      <-  area[match(paste0(tableS1$studyName, tableS1$sites), names(area))]
	tableS1$richness  <-  richness[match(paste0(tableS1$studyName, tableS1$sites), names(richness))]
	names(tableS1)    <-  c('Province', 'Sub-province', 'Site', 'Longitude', 'Latitude', 'Main contact', 'New Site', 'Sampling effort', 'Transect area (m2)', 'Richness')
    for (z in 1:nrow(tableS1)) {
    	for (i in c('Province', 'Sub-province', 'Site')) {
		    scaleLabel  <-  strsplit(tableS1[[i]][z], '_')[[1]]
		    if (length(scaleLabel) > 1) {
		        for (j in seq_along(scaleLabel)) {
		            chars  <-  nchar(scaleLabel[j])
		            scaleLabel[j]  <-  ifelse(chars <= 3, toupper(scaleLabel[j]), Hmisc::capitalize(scaleLabel[j]))
		        }
		        scaleLabel  <-  paste0(scaleLabel, collapse = ' ')
		    } else {
		        scaleLabel  <-  Hmisc::capitalize(scaleLabel)
		    }
		    tableS1[[i]][z]  <-  scaleLabel
    	}
	}
	authors  <-  c('agreen_palau' = 'AL Green', 'agreen_raja_ampat_2012' = 'AL Green', 'agreen_samoa' = 'AL Green', 'agreen_solomon_2008' = 'AL Green', 'atlantic_tep' = 'SR Floeter', 'caribbean' = 'FA Rodriguez-Zaragoza', 'clipperton_alan_2017' = 'AM Friedlander', 'cocos' = 'AM Friedlander', 'curacao' = 'SR Floeter', 'easter_salaz' = 'AM Friedlander', 'fiji_biomass_2013a' = 'M Kulbicki', 'flores_1993_michel' = 'M Kulbicki', 'fpolynesia_biomass_2013' = 'M Kulbicki', 'hawaii' = 'AM Friedlander', 'loyalty' = 'AM Friedlander', 'mexico' = 'SR Floeter', 'mozambique_alan_2017' = 'AM Friedlander', 'new_caledonia_biomass_2013c1' = 'M Kulbicki', 'pitcairn' = 'AM Friedlander', 'revillagigedo_alan_2017' = 'AM Friedlander', 'rls_march_2017' = 'GJ Edgar, RD Stuart-Smith', 'seychelles_alan_2017' = 'AM Friedlander', 'tonga_biomass_2013' = 'M Kulbicki')
	tableS1[['Main contact']]  <-  authors[tableS1[['Main contact']]]
	tableS1[['Main contact']][tableS1[['Main contact']] == 'FA Rodriguez-Zaragoza' & tableS1[['Sub-province']] %in% c('Mexico Alacranes', 'Mexico Caribbean')]  <-  'JE Arias-Gonzalez'
	tableS1[['Main contact']][tableS1[['Main contact']] == 'FA Rodriguez-Zaragoza' & tableS1[['Sub-province']] %in% c('US Virgin Islands')]  <-  'AM Friedlander'
	tableS1
}

makeTableS2  <-  function (figClist) {
	tableS2  <-  tapply(figClist$Gaspar_code, figClist$Site_name, lenUn)
	tableS2  <-  data.frame(Province = figClist$Realm[match(names(tableS2), figClist$Site_name)], Site = names(tableS2), Longitude = figClist$Longitude[match(names(tableS2), figClist$Site_name)], Latitude = figClist$Latitude[match(names(tableS2), figClist$Site_name)], Richness = tableS2, row.names = NULL, stringsAsFactors = FALSE)
	tableS2  <-  tableS2[order(tableS2$Province), ]
    for (z in 1:nrow(tableS2)) {
	    regionLab  <-  strsplit(tableS2$Province[z], '_')[[1]]
	    if (length(regionLab) > 1) {
	        for (j in seq_along(regionLab)) {
	            chars  <-  nchar(regionLab[j])
	            regionLab[j]  <-  ifelse(chars <= 3, toupper(regionLab[j]), Hmisc::capitalize(regionLab[j]))
	        }
	        regionLab  <-  paste0(regionLab, collapse = ' ')
	    } else {
	        regionLab  <-  Hmisc::capitalize(regionLab)
	    }
	    tableS2$Province[z]  <-  regionLab
	}
	tableS2
}
