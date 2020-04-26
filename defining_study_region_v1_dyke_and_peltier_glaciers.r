### defining study region for green ash
### Adam B. Smith | Missouri Botanical Garden | 2019-12
### source('C:/ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/study_region/source('C:/ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/study_region/defining_study_region_v2_dalton_glaciers.r').r')
### source('D:/ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/study_region/source('C:/ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/study_region/defining_study_region_v2_dalton_glaciers.r').r')

### This script delineates a geographic domain to be used for simulating the biogeographic history of green ash.  As per meetings on 2019-11-12 and 2019-11-19, we will demarcate the region as such:
### * West: the North American continental divide
### * South: southern extent of the state of Tamaulipas, Mexico
### * East: Atlantic ocean (extended outward to account for seal level rise since 21 Kybp)
### * North: Northward to the extent of the Hudson Bay drainage

### Data sources:
### * Level I and Level II Watershed Boundaries from the USGS/Council on Environmental Cooperation (https://www.sciencebase.gov/catalog/item/4fb697b2e4b03ad19d64b47f)
### * GADM 3.6 political geography (https://gadm.org/)
### * Lorenz et al. climate layers to demarcate land versus sea since the LGM (Lorenz et al. 2016.  Downscaled and debiased climate simulations for North America from 21,000 years ago to 2100AD.  Scientific Data 3:160048.)

### CONTENTS ###
### setup ###
### copy/convert external data ###
### get pollen data for Fraxinus ###
### retrieve occurrence data ###
### create spatial polygon of study extent ###
### mask Dalton et al 2020 ice sheet layers onto land mass from Lorenz et al 2016 ###


### create interpolations of Dyke (vectorized) glacier layers ###
### maps of Dyke (vectorized) glaciers ###
### create interpolations of Peltier (rasterized) glacier layers ###
### maps of Peltier (rasterized) and Dyke (vectorized) glaciers ###


#############
### setup ###
#############

	memory.limit(memory.limit() * 2^30)
	rm(list=ls())
	gc()
	options(stringsAsFactors=FALSE)

	library(sp)
	library(raster)
	library(dismo)
	library(rgeos)
	library(neotoma)
	library(scales)
	library(ncdf4)
	# library(spatialEco)
	
	library(omnibus) # Adam's custom library on GitHub: adamlilith/omnibus
	library(enmSdm) # Adam's custom library on GitHub: adamlilith/enmSdm
	library(birdsEye) # Adam's custom library on GitHub: adamlilith/birdsEye

	rasterOptions(format='GTiff', overwrite=TRUE)
	
	library(compiler)
	enableJIT(1)
	setCompilerOptions(suppressUndefined=TRUE)

	drive <- 'C:'
	# drive <- 'D:'
	
	setwd(paste0(drive, '/ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/study_region'))

	# settings for glacial tweening... see birdsEye::interPolysByTween()
	method <- 'linear'
	# method <- 'cubic-in-out'
	delta <- 1000
	# delta <- 100

	# raster interpolation settings... see enmSdm::interpolateRasters()
	# rastInterpFx <- 'linear'
	rastInterpFx <- 'spline'
	
# ##################################
# ### copy/convert external data ###
# ##################################
	
	# # watershed boundaries from CEC
	# ws1 <- shapefile(paste0(drive, '/ecology/Watersheds (CEC)/NA_Watersheds_Level_I_Ocean_Drainage_Areas'))
	# ws23 <- shapefile(paste0(drive, '/ecology/Watersheds (CEC)/NA_Watersheds_Level_II_and_III'))

	# dirCreate('./watersheds_cec')
	# save(ws1, file='./watersheds_cec/watersheds_level_1.rda')
	# save(ws23, file='./watersheds_cec/watersheds_level_2_and_3.rda')

# ####################################
# ### get pollen data for Fraxinus ###
# ####################################

	# # Using this data in ArcMAP to visualize candidate study region.

	# availFraxData <- get_dataset(loc = c(-145, 10, -50, 71),
		# datasettype = 'pollen',
		# taxonname = 'Fraxinus%'
	# )

	# fraxData <- get_download(availFraxData, verbose = FALSE)
	
	# dirCreate('./fraxinus_pollen_neotoma')
	# save(fraxData, file='/fraxinus_pollen_neotoma/fraxinus_pollen_data_2020_01_03.rda')

	# compiledFrax <- compile_taxa(fraxData, 'P25')

	# fraxBySite <- data.frame()
	# for (i in seq_along(compiledFrax)) {
	
		# counts <- compiledFrax[[i]]$counts
		# whichCol <- grepl(pattern='Fraxinus', colnames(counts))
		
		# if (any(whichCol)) {
			
			# frax <- counts[ , whichCol]
			# sums <- rowSums(counts, na.rm=TRUE)
			# fraxPerc <- frax / sums
			# maxFraxProp <- max(fraxPerc, na.rm=TRUE)
			
		# } else {
			
			# maxFraxProp <- 0
			
		# }
		
		# fraxBySite <- rbind(
			# fraxBySite,
			# data.frame(
				# longitude = compiledFrax[[i]]$dataset$site$long,
				# latitude = compiledFrax[[i]]$dataset$site$lat,
				# maxFraxinusProp = maxFraxProp
			# )
		# )
		
	# }
	
	# fraxBySite <- SpatialPointsDataFrame(fraxBySite [ , c('longitude', 'latitude')], data=fraxBySite, proj4=getCRS('wgs84', TRUE))
	# shapefile(fraxBySite, './fraxinus_pollen_neotoma/fraxinus_pollen_data_2020_01_03_proportion_by_site')
	
# ################################
# ### retrieve occurrence data ###
# ################################

	# # Occurrence data is presumed to have been already downloaded (and possibly cleaned) using scripts in "enms" folder.

	# load(paste0(drive, '/ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/enms/species_records/00_Fraxinus_pennsylvanica_bien_all_occurrences.rda'))
	
	# occsRaw <- SpatialPoints(occsRaw[ , c('longitude', 'latitude')], getCRS('wgs84', TRUE))
	
	# dirCreate('./occurrences_bien')
	# shapefile(occsRaw, './occurrences_bien/00_Fraxinus_pennsylvanica_bien_all_occurrences')
	
# say('##############################################')
# say('### create spatial polygon of study extent ###')
# say('##############################################')

	# # watershed boundaries
	# load('./watersheds_cec/watersheds_level_1.rda')
	# load('./watersheds_cec/watersheds_level_2_and_3.rda')

	# # political boundaries from GADM
	# can <- getData('GADM', country='CAN', level=1, path='C:/ecology/!Scratch')
	# usa <- getData('GADM', country='USA', level=1, path='C:/ecology/!Scratch')
	# mex <- getData('GADM', country='MEX', level=1, path='C:/ecology/!Scratch')
	
	# nam <- rbind(can, usa, mex)
	# nam <- sp::spTransform(nam, CRS(projection(ws1)))
	
	# # clean watershed names
	# ws1$NAW1_EN <- trim(ws1$NAW1_EN)
	# ws1$NAW2_EN <- trim(ws1$NAW2_EN)
	# ws1$NAW3_EN <- trim(ws1$NAW3_EN)
	# ws1$NAW4_EN <- trim(ws1$NAW4_EN)

	# # major drainage areas
	# demesne <- ws1[ws1@data$NAW1_EN %in%  c('Atlantic Ocean', 'Hudson Bay', 'Gulf of Mexico'), ]
	
	# # remove Rio Grande drainage area *except* keep NAE Level 2 representing drainage into Gulf of Mexico to ensure spatial continuity
	# rioGrande <- demesne[demesne@data$NAW2_EN == 'Rio Grande', ]
	# lowerRioGrande <- rioGrande[rioGrande@data$NAW4_EN == 'Lower Rio Grande', ]
	
	# demesne <- demesne[demesne@data$NAW2_EN != 'Rio Grande', ]
	# demesne <- rbind(demesne, lowerRioGrande)
	
	# # add back in part of Rio Grande so study region is not "pinched" at southern end
	# keeps <- ws1[ws1@data$OBJECTID %in% c(284, 2228), ]
	# demesne <- rbind(demesne, keeps)
	
	# # remove El Salado basin
	# demesne <- demesne[demesne@data$NAW2_EN != 'El Salado', ]
	
	# # remove islands and areas mostly south of Tamaulipas
	# demesne <- demesne[!(demesne@data$NAW4_EN %in% c('Puerto Rico', 'Virgin Islands')), ]
	
	# demesne <- demesne[!(demesne@data$NAW4_EN %in% c('Pánuco', 'Tuxpan - Nautla', 'Papaloapan', 'Coatzacoalcos', 'Grijalva - Usumacinta', 'Yucatán Oeste (Campeche)', 'Yucatán Norte (Yucatán)')), ]
	# demesne <- demesne[!(demesne@data$OBJECTID %in% c(66, 79, 92, 95, 97, 98, 99, 130)), ] # 198?
	
	# # remove Sable Island
	# demesne <- demesne[!(demesne@data$OBJECTID %in% 2078), ]
	
	# # remove Baffin island and Southhampton Island
	# demesne <- demesne[!(demesne@data$NAW4_EN %in% c('Hudson Strait - Baffin and Southampton Islands', 'Hudson Bay - Southampton Island', 'Foxe Basin - Baffin Island', 'Foxe Basin - Southampton Island')), ]
	
	# # remove Melville Peninsula
	# demesne <- demesne[!(demesne@data$NAW4_EN %in% c('Foxe Basin - Melville Peninsula')), ]

	# # remove Akpotak Island
	# demesne <- demesne[!(demesne@data$OBJECTID %in% c(1772)), ]

	# # combine watershed-delineated study region with manually-drawn region that encompasses adjacent littorial area that was exposed during LGM
	# littoral <- shapefile('formerly_exposed_land_from_lorenz_et_al_2016/manually_drawn_study_region_encompassing_continental_shelf')
	# demesne <- gUnaryUnion(demesne)
	
	# demesne <- gUnion(demesne, littoral)
	# demesne <- spatialEco::remove.holes(demesne)

	# # load species data
	# occsRaw <- shapefile('./occurrences_bien/00_Fraxinus_pennsylvanica_bien_all_occurrences')
	# pollen <- shapefile('./fraxinus_pollen_neotoma/fraxinus_pollen_data_2020_01_03_proportion_by_site')
	# genetics <- shapefile(paste0(drive, '/Ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/genetic_data_locations_allan_strand_2019_11_12'))
	# dirCreate('./genetic_data')
	# save(genetics, file='./genetic_data/genetic_data_locations_2019_11_12')
	# shapefile(genetics, './genetic_data/genetic_data_locations_2019_11_12', overwrite=TRUE)
	
	# bien <- shapefile(paste0(drive, '/Ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/enms/regions/bien_range_map/Fraxinus_pennsylvanica'))
	# little <- shapefile(paste0(drive, '/Ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/enms/regions/littles_range_map/fraxpenn'))
	
	# occsRaw <- sp::spTransform(occsRaw, CRS(projection(demesne)))
	# pollen <- sp::spTransform(pollen, CRS(projection(demesne)))
	# genetics <- sp::spTransform(genetics, CRS(projection(demesne)))

	# bien <- sp::spTransform(bien, CRS(projection(demesne)))
	# little <- sp::spTransform(little, CRS(projection(demesne)))
	
	# png('./study_region_by_watershed_no_glaciers.png', width=1600, height=1200, res=300)
	
		# par(oma=c(0, 0, 0, 0), mar=rep(0.5, 4))
		# plot(demesne, col='darkolivegreen3', border='darkolivegreen')
		# plot(nam, border='black', add=TRUE)
		
		# plot(bien, lwd=2, border='darkblue', add=TRUE)
		# plot(little, lwd=2, border='darkgreen', add=TRUE)
		
		# points(occsRaw, pch=16, cex=0.25, col='darkorange')
		
		# points(pollen, pch=2, cex=0.6)
		# points(pollen[pollen$mxFrxnP >= 0.05, ], pch=24, cex=0.6, bg='yellow')
		
		# points(genetics, pch=22, bg='red', cex=0.8)

		# legend('bottomright', inset=-0.01, bty='n', cex=0.6,
			# legend=c(
				# 'study region',
				# 'Little\'s range',
				# 'BIEN range',
				# 'occurrences',
				# 'pollen with Fraxinus',
				# 'pollen with Fraxinus >= 5%',
				# 'genetic data'
			# ),
			# pch=c(
				# NA,
				# NA,
				# NA,
				# 16,
				# 2,
				# 24,
				# 22
			# ),
			# col=c(
				# NA,
				# 'darkgreen',
				# 'darkblue',
				# 'darkorange',
				# 'black',
				# 'black',
				# 'black'
			# ),
			# border=c(
				# 'darkolivegreen',
				# NA,
				# NA,
				# NA,
				# NA,
				# NA,
				# NA				
			# ),
			# fill=c(
				# 'darkolivegreen3',
				# NA,
				# NA,
				# NA,
				# NA,
				# NA,
				# NA
			# ),
			# pt.bg=c(
				# NA,
				# NA,
				# NA,
				# NA,
				# NA,
				# 'yellow',
				# 'red'
			# ),
			# lwd=c(
				# NA,
				# 2,
				# 2,
				# NA,
				# NA,
				# NA,
				# NA
			# )
		# )
		
		# title(sub=date(), line=0, cex.sub=0.3)
		
	# dev.off()
	
	# shapefile(demesne, './study_region_mask_without_glaciers', overwrite=TRUE)

say('#####################################################################################')
say('### mask Dalton et al 2020 ice sheet layers onto land mass from Lorenz et al 2016 ###')
say('#####################################################################################')	
	
	say('This step will create rasters of the land mass of North America through time using
	climate layers from Lorenz et al 2016 Sci Data. We will overlay the ice sheet layers from
	Dalton et al 2020 QSR onto the layers to yiel cell values ranging from 0 (no ice) to 1
	(complete ice cover). Andria will then try a spatiotemporal interpolation using these
	layers.')
	
	srcDir <- 'D:/Ecology/Climate/Lorenz et al 2016 North America 21Kybp to 2100 CE/CCSM/'
	
	lorenz <- raster(paste0(srcDir, '/22000BP/an_avg_ETR.tif')
	
	years <- seq(22000 - 500, 0, by=-500)
	
	for (year in years) {
		lorenz <- stack(
			lorenz,
			raster(paste0(srcDir, '/', year, 'BP/an_avg_ETR.tif'))
		)
	)
	
	
	
	
	
	
	
	
	
	
# say('#################################################################')
# say('### create interpolations of Dyke (vectorized) glacier layers ###')
# say('#################################################################')
	
	# dirCreate(paste0('./glaciers/glaciers_dyke/', method, '_delta', delta))
	
	# datesCrosswalk <- read.csv(paste0(drive, '/Ecology/Drive/Data/Paleoglaciers of North America - Dyke/Dates Crosswalk.csv'))
	
	# calYears <- seq(21000, 0, by=-30)
	
	# # for each time span
	# for (i in nrow(datesCrosswalk):2) {
		
		# say('   Tweening glaciers for ', datesCrosswalk$calendarYbp[i], ' to ', datesCrosswalk$calendarYbp[i - 1], '...')
		
		# # get years... missing ice for 8400 RCYBP
		# if (datesCrosswalk$iceLayers[i] == 'missing') {

			# startCalYear <- datesCrosswalk$calendarYbp[i + 1]
			# endCalYear <- datesCrosswalk$calendarYbp[i - 1]
			
			# startRadioYear <- datesCrosswalk$iceRadioYbp[i + 1]
			# endRadioYear <- datesCrosswalk$iceRadioYbp[i - 1]
		
		# # get years... missing ice for 8400 RCYBP		
		# } else if (datesCrosswalk$iceLayers[i - 1] == 'missing') {

			# # get years
			# startCalYear <- datesCrosswalk$calendarYbp[i]
			# endCalYear <- datesCrosswalk$calendarYbp[i - 2]
			
			# startRadioYear <- datesCrosswalk$iceRadioYbp[i]
			# endRadioYear <- datesCrosswalk$iceRadioYbp[i - 2]
		
		# # get years
		# } else {

			# # get years
			# startCalYear <- datesCrosswalk$calendarYbp[i]
			# endCalYear <- datesCrosswalk$calendarYbp[i - 1]
			
			# startRadioYear <- datesCrosswalk$iceRadioYbp[i]
			# endRadioYear <- datesCrosswalk$iceRadioYbp[i - 1]
			
		# }
			
		# say('interpolating from calendar year ', startCalYear, ' ybp to ', endCalYear, ' ybp', level=2)
		
		# startIce <- shapefile(paste0(drive, '/Ecology/Drive/Data/Paleoglaciers of North America - Dyke/IceMaps_1ka_All/', prefix(startRadioYear, 5), 'icell'))
		# endIce <- shapefile(paste0(drive, '/Ecology/Drive/Data/Paleoglaciers of North America - Dyke/IceMaps_1ka_All/', prefix(endRadioYear, 5), 'icell'))
		
		# theseCalYears <- calYears[calYears <= startCalYear & calYears >= endCalYear]
		# betweens <- (startCalYear - theseCalYears) / (startCalYear - endCalYear)

		# ices <- interpPolysByTween(
			# x1 = startIce,
			# x2 = endIce,
			# eaCrs = 'oae',
			# between = betweens,
			# delta = delta,
			# method = method,
			# verbose=FALSE
		# )
		
		# for (i in seq_along(ices)) {
		
			# ice <- ices[[i]]$poly
			# save(ice, file=paste0('./glaciers/glaciers_dyke/', method, '_delta', delta, '/ice_sheet_', prefix(theseCalYears[i], 5), '_ybp.rda'))
		
		# }
			
	# } # next time period

# say('##########################################')
# say('### maps of Dyke (vectorized) glaciers ###')
# say('##########################################')	

	# dirCreate('./glaciers/glaciers_dyke/', method, '_delta', delta, '/maps_by_year')
	# calYears <- seq(21000, 5010, by=-30)

	# # political boundaries from GADM
	# can <- getData('GADM', country='CAN', level=1, path='C:/ecology/!Scratch')
	# usa <- getData('GADM', country='USA', level=1, path='C:/ecology/!Scratch')
	# mex <- getData('GADM', country='MEX', level=1, path='C:/ecology/!Scratch')
	# nam <- rbind(can, usa, mex)

	# # study region
	# study <- shapefile('./study_region_mask_without_glaciers')
	
	# # starting ice
	# load(paste0('./glaciers/glaciers_dyke/', method, '_delta', delta, '/ice_sheet_21000_ybp.rda'))
	# ice21000 <- ice
	
	# # equal-area CRS
	# crs <- proj4string(ice21000)
	# ell <- ellipsoid(crs)
	# dat <- datum(crs)
	# cent <- geosphere::centroid(ice21000)
	# cent <- coordinates(cent)
	# eaCrs <- makeCRS('oae', long0=cent[1, 1], lat0=cent[1, 2], asCRS=TRUE)

	# # project
	# ice21000 <- sp::spTransform(ice21000, eaCrs)
	# nam <- sp::spTransform(nam, eaCrs)
	# study <- sp::spTransform(study, eaCrs)
	
	# # plot extent
	# ext <- extent(ice21000)
	# ext <- as(ext, 'SpatialPolygons')
	# projection(ext) <- as.character(eaCrs)
	
	# for (calYear in calYears) {
	
		# say('Making map for ', calYear, '...')
	
		# load(paste0('./glaciers/glaciers_dyke/', method, '_delta', delta, '/ice_sheet_', prefix(calYear, 5), '_ybp.rda'))
		# ice <- sp::spTransform(ice, eaCrs)
	
		# png(paste0('./glaciers/glaciers_dyke/', method, '_delta', delta, '/maps_by_year/', prefix(21000 - calYear, 5), ' Years After 21000 YBP.png'), width=1500, height=1000)
			# plot(ext, border=NA)
			# plot(nam, col='gray80', add=TRUE)
			# plot(study, col=alpha('orange', 0.3), add=TRUE)
			# plot(ice21000, border='turquoise4', col=alpha('turquoise4', alpha=0.2), lwd=2, add=TRUE)
			# plot(study, border='orange', lwd=2, add=TRUE)
			# plot(ice, col=alpha('cadetblue1', 0.8), add=TRUE)
			# title(main=paste(calYear, ' YBP'), sub=date(), cex.main=3.6, cex.sub=1, line=1)
		# dev.off()	
	
	# }

# say('####################################################################')
# say('### create interpolations of Peltier (rasterized) glacier layers ###')
# say('####################################################################')
	
	# dirCreate('./glaciers/glaciers_peltier/rastInterpFx_', rastInterpFx)
	
	# ### years
	
		# peltierYears <- c(seq(0, 20, by=0.5), 21) # base peltier glacier years
		# calYears <- seq(21000, 0, by=-30) # years for which we want plots

	# ### process peltier glaciers
	
		# if (exists('peltier')) rm(peltier)
		# for (peltierYear in peltierYears) {
		
			# inFileName <- paste0(drive, '/ecology/Drive/Data/Paleoglaciers of the World - Peltier/ICE-7G_NA/I7G_NA.VM7_1deg.', peltierYear, '.nc')
			
			# if (!exists('peltier')) {
				# peltier <- raster(inFileName, ncdf=TRUE)
			# } else {
				# thisRast <- raster(inFileName, ncdf=TRUE)
				# peltier <- stack(peltier, thisRast)
			# }
			
			# names(peltier)[nlayers(peltier)] <- paste0('peltier', prefix(1000 * peltierYear, 5), 'ybp')
		
		# }
		
		# # rotate to long = c(-180, 180)
		# peltier <- rotate(peltier)

		# # crop to Northern/Western Hemisphere
		# ext <- extent(c(-180, 0, 0, 90))
		# peltier <- crop(peltier, ext)
		
		# # interpolate
		# peltierInterp <- interpolateRasters(peltier, interpFrom=1000 * peltierYears, interpTo=calYears, type=rastInterpFx)
		# names(peltierInterp) <- paste0('peltier', prefix(calYears, 5), 'ybp')

		# for (i in 1:nlayers(peltierInterp)) {

			# rast <- peltierInterp[[i]]
			# name <- names(rast)
			# writeRaster(rast, paste0('./glaciers/glaciers_peltier/rastInterpFx_', rastInterpFx, '/', name), datatype='FLT4S')
		
		# }
	
# say('###################################################################')
# say('### maps of Peltier (rasterized) and Dyke (vectorized) glaciers ###')
# say('###################################################################')

	# outDir <- paste0('./glaciers/glaciers_peltier_and_dyke_with_method_', method, '_delta_', delta)
	# dirCreate(outDir)
	
	# ### years
	
		# calYears <- seq(0, 21000, by=30) # years for which we want plots

	# ## other data for plotting

		# # political boundaries from GADM
		# can <- getData('GADM', country='CAN', level=1, path='C:/ecology/!Scratch')
		# usa <- getData('GADM', country='USA', level=1, path='C:/ecology/!Scratch')
		# mex <- getData('GADM', country='MEX', level=1, path='C:/ecology/!Scratch')
		# nam <- rbind(can, usa, mex)

		# # study region
		# study <- shapefile('./study_region_mask_without_glaciers')

	# ### establish extent for plotting using Dyke glaciers from 21 Kybp
	
		# # starting ice
		# dyke21000 <- shapefile(paste0(drive, '/ecology/Drive/Data/Paleoglaciers of North America - Dyke/IceMaps_1ka_All/18000icell'))
		
		# # equal-area CRS
		# crs <- proj4string(dyke21000)
		# ell <- ellipsoid(crs)
		# dat <- datum(crs)
		# cent <- geosphere::centroid(dyke21000)
		# cent <- coordinates(cent)
		# eaCrs <- makeCRS('oae', long0=cent[1, 1], lat0=cent[1, 2], ell='GRS80', dat='NAD83', asCRS=TRUE)

		# # plot extent
		# ext <- extent(dyke21000)
		# ext <- as(ext, 'SpatialPolygons')
		# projection(ext) <- projection(dyke21000)

		# ext <- sp::spTransform(ext, eaCrs)
		# dyke21000 <- sp::spTransform(dyke21000, eaCrs)
		# nam <- sp::spTransform(nam, eaCrs)
		# study <- sp::spTransform(study, eaCrs)
		
		# ext <- gUnion(ext, study)
		# ext <- as(ext, 'SpatialPolygons')
		# projection(ext) <- projection(eaCrs)
		
	# ### process peltier glaciers
	
		# peltierLinear <- raster::stack(listFiles(paste0('./glaciers/glaciers_peltier/rastInterpFx_linear'), pattern='.tif'))
		# peltierSpline <- raster::stack(listFiles(paste0('./glaciers/glaciers_peltier/rastInterpFx_spline'), pattern='.tif'))
	
		# # crop to Dyke glacial extent
		# peltierLinear <- crop(peltierLinear, ext)
		# peltierSpline <- crop(peltierSpline, ext)
		
		# # project to equal-area CRS for plotting
		# beginCluster(4)
		# peltierLinear <- projectRaster(peltierLinear, crs=as.character(eaCrs))
		# peltierSpline <- projectRaster(peltierSpline, crs=as.character(eaCrs))
		# endCluster()
		
		# # force minimum value to 0
		# for (i in 1:nlayers(peltierLinear)) {
			# peltierLinear[[i]] <- calc(peltierLinear[[i]], fun=function(x) ifelse(x < 0, 0, x))
			# peltierSpline[[i]] <- calc(peltierSpline[[i]], fun=function(x) ifelse(x < 0, 0, x))
		# }
		
		# names(peltierLinear) <- paste0('peltier', prefix(calYears, 5), 'ybp')
		# names(peltierSpline) <- paste0('peltier', prefix(calYears, 5), 'ybp')
		
		# # settings for plotting
		# minPeltier <- min(minValue(peltierLinear), minValue(peltierSpline))
		# maxPeltier <- max(maxValue(peltierLinear), maxValue(peltierSpline))

		# breaks <- seq(0, maxPeltier, length.out=9)

	# ### plot

		# blues <- c('#deebf7', '#c6dbef', '#9ecae1', '#6baed6', '#4292c6', '#2171b5', '#08519c', '#08306b')
		
		# for (calYear in calYears) {
		
			# say('Making map for ', calYear, '...')
		
			# if (calYear >= 5000) {
			
				# load(paste0('./glaciers/glaciers_dyke/', method, '_delta', delta, '/ice_sheet_', prefix(calYear, 5), '_ybp.rda'))
				# dyke <- sp::spTransform(ice, eaCrs)
				
			# }
				
			# png(paste0(outDir, '/', prefix(21000 - calYear, 5), ' Years After 21000 YBP.png'), width=1600, height=800)
		
				# par(mfrow=c(1, 2), mar=rep(0, 4))
		
				# ## linear Peltier
				# plot(ext, border=NA)
				# plot(nam, col='gray80', add=TRUE)
				# plot(study, col=alpha('orange', 0.2), add=TRUE)

				# thisPeltier <- peltierLinear[[paste0('peltier', prefix(calYear, 5), 'ybp')]]
				# littleIce <- calc(thisPeltier, fun=function(x) ifelse(x > 0, 1, NA))
				# peltierInterpMasked <- peltierLinear[[paste0('peltier', prefix(calYear, 5), 'ybp')]] * littleIce
				# plot(peltierInterpMasked, breaks=breaks, col=blues, add=TRUE, ann=FALSE, legend=FALSE)
				
				# plot(study, border='orange', lwd=2, add=TRUE)
				# if (calYear >= 5000) {
					# plot(dyke, border='cadetblue', lwd=3, add=TRUE)
				# }

				# title(main='Peltier with Linear Interpolation\nDyke with Tweened Interpolation', cex.main=2, line=-6)
				
				# legend('bottomleft', inset=0.01, legend=c('Dyke glaciers', 'Peltier glaciers', 'Study region'), fill=c(NA, blues[5], alpha('orange', 0.2)), border=c(NA, NA, 'orange'), lwd=c(3, NA, NA), col=c('cadetblue4', NA, NA), bty='n', cex=2.2)

				# ## spline Peltier
				# plot(ext, border=NA)
				# plot(nam, col='gray80', add=TRUE)
				# plot(study, col=alpha('orange', 0.2), add=TRUE)

				# thisPeltier <- peltierSpline[[paste0('peltier', prefix(calYear, 5), 'ybp')]]
				# littleIce <- calc(thisPeltier, fun=function(x) ifelse(x > 0, 1, NA))
				# peltierInterpMasked <- peltierSpline[[paste0('peltier', prefix(calYear, 5), 'ybp')]] * littleIce
				# plot(peltierInterpMasked, breaks=breaks, col=blues, add=TRUE, ann=FALSE, legend=FALSE)
				
				# plot(study, border='orange', lwd=2, add=TRUE)
				# if (calYear >= 5000) {
					# plot(dyke, border='cadetblue', lwd=3, add=TRUE)
				# }
				
				# title(main='Peltier with Spline Interpolation\nDyke with Tweened Interpolation', cex.main=2, line=-6)
				
				# title(main=paste(calYear, ' YBP'), cex.main=3.6, cex.sub=1, line=-5.5, outer=TRUE)
				# title(sub=date(), cex.sub=1, outer=TRUE, line=-2)
				
			# dev.off()	

		# }
	
#################################	
say('DONE!!!', deco='%', level=1)
#################################
