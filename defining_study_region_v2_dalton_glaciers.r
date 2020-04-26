### defining study region for green ash
### Adam B. Smith | Missouri Botanical Garden | 2019-12
### source('C:/ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/study_region/defining_study_region_v2_dalton_glaciers.r')
### source('D:/ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/study_region/defining_study_region_v2_dalton_glaciers.r')

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

	# raster interpolation settings... see enmSdm::interpolateRasters()
	rastInterpFx <- 'linear'
	# rastInterpFx <- 'spline'
	
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

	# say('This subscript uses watershed boundaries to delineate a study region extent that covers the area of interest for projecting range dynamics of green ash from 21 Kybp to the present in eastern North America.')

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

# say('#####################################################################################')
# say('### mask Dalton et al 2020 ice sheet layers onto land mass from Lorenz et al 2016 ###')
# say('#####################################################################################')	
	
	# say('This step will create rasters of the land mass of North America through time using
	# climate layers from Lorenz et al 2016 Sci Data. We will overlay the ice sheet layers from
	# Dalton et al 2020 QSR onto the layers to yield cell values ranging from 0 (no ice) to 1
	# (complete ice cover).')

	# ### get rasters representing land from Lorenz et al

		# lorenzDir <- 'D:/Ecology/Climate/Lorenz et al 2016 North America 21Kybp to 2100 CE/CCSM/'
		
		# lorenz <- raster(paste0(lorenzDir, '/22000BP/an_avg_ETR.tif'))
		
		# years <- seq(22000 - 500, 0, by=-500)
		
		# for (year in years) {
			# lorenz <- stack(
				# lorenz,
				# raster(paste0(lorenzDir, '/', year, 'BP/an_avg_ETR.tif'))
			# )
		# }
		
		# lorenz <- lorenz * 0 + 1
		# lorenzYears <- seq(22000, 0, by=-500)
		# names(lorenz) <- paste0('yr', lorenzYears, 'bp')
	
	# ### get carbon and calendar years for each ice layer for Dalton et al ice sheet data
	# # note that the calendar years are not exactly the same as listed on Table 1 in Dalton et al 2020
		
		# daltonYears <- read.csv('C:/Ecology/Drive/Data/North American Ice Sheet Dalton et al 2020 Quaternary Science Reviews/Dalton et al 2020 QSR Dates from Shapefile Names.csv')

	# ### interpolate land rasters to time periods matching Durant et al ice sheet representations

		# interpTo <- -1000 * daltonYears$calKiloYear
		# lorenzInterp <- interpolateRasters(lorenz, interpFrom=-1 * lorenzYears, interpTo=interpTo, type=rastInterpFx)
		
		# names(lorenzInterp) <- paste0('yr', 1000 * daltonYears$calKiloYear, 'bp')
	
	# ### for each cell, assign a value from 0 to 1 indicating proportion covered by ice sheet
	
		# daltonDir <- 'C:/Ecology/Drive/Data/North American Ice Sheet Dalton et al 2020 Quaternary Science Reviews/RDA Files/'
	
		# for (countDalton in 1:nrow(daltonYears)) {

			# # ice shapefile
			# daltonCalYear <- daltonYears$calKiloYear[countDalton]
			# load(paste0(daltonDir, 'daltonEtAl2020_', sprintf('%02.2f', daltonCalYear), '_kiloCalYBP.rda'))
			
			# say('Extracting to ', daltonCalYear, ' kilo calendar years BP')
			
			# # extract, remember proportion of cell covered by ice
			# iceOnRast <- raster::extract(lorenzInterp[[countDalton]], daltonIce, cellnumbers=TRUE, weights=TRUE, normalizeWeights=FALSE)
		
			# # transfer values back to raster
			# vals <- values(lorenzInterp[[countDalton]])
			# vals <- vals * 0
			
			# for (countIce in seq_along(iceOnRast)) {
				
				# vals[iceOnRast[[countIce]][ , 'cell']] <- vals[iceOnRast[[countIce]][ , 'cell']] + iceOnRast[[countIce]][ , 'weight']
				
			# }
			
			# lorenzInterp[[countDalton]] <- setValues(lorenzInterp[[countDalton]], values=vals)
			
		# }
		
	# ### add "year 0" to Lorenz terrestrial
	
		# year0 <- lorenz[[nlayers(lorenz)]]
		# lorenzInterp <- stack(lorenzInterp, year0)

		# dirCreate('./glaciers/glaciers_dalton')
		# writeRaster(lorenzInterp, paste0('./glaciers/glaciers_dalton/daltonGlaciers_on_lorenzLand_', rastInterpFx, 'Interpolation'))
		# file.copy('C:/Ecology/Drive/Data/North American Ice Sheet Dalton et al 2020 Quaternary Science Reviews/Dalton et al 2020 QSR Dates from Shapefile Names.csv', to='C:/Ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/study_region/glaciers/glaciers_dalton/Dalton et al 2020 QSR Dates from Shapefile Names.csv')
		
		# sink(paste0('./glaciers/glaciers_dalton/README_', rastInterpFx, 'Interpolation.txt'), split=TRUE)
		
			# say('The raster stack in this folder represents land area used by Lorenz et al 2016 Sci Data overlaid with ice sheet coverage from Dalton et al 2020 Quat Sci Reviews. Values range from 0 (no portion of the cell covered by ice) to 1 (all of the cell covered by ice). Approximate CALENDAR year represented by each layer in the stack is provided by the layer name AND in the table "Dalton et al 2020 QSR Dates from Shapefile Names.csv". The ', rastInterpFx, ' interpolation method was used.', breaks=80)
			
		# sink()

say('#############################')
say('### generate study region ###')
say('#############################')

	### interpolate Dalton ice sheet layer to every 30 yr
	#####################################################

	dalton <- stack(paste0('./glaciers/glaciers_dalton/daltonGlaciers_on_lorenzLand_', rastInterpFx, 'Interpolation.tif'))
	daltonDates <- read.csv('./glaciers/glaciers_dalton/Dalton et al 2020 QSR Dates from Shapefile Names.csv')

	calYearsFrom <- c(-1 * daltonDates$calKiloYear * 1000, 0)
	calYearsTo <- seq(-21000, 0, by=30)
	
	names(dalton) <- paste0('yr', abs(calYearsFrom), 'bp')
	
	daltonInterp <- interpolateRasters(dalton, interpFrom=calYearsFrom, interpTo=calYearsTo, type=rastInterpFx)
	names(daltonInterp) <- paste0('yr', abs(calYearsTo), 'bp')

	### mask with study region polygon created above
	################################################
	
	studyRegionPoly <- shapefile('./study_region_mask_without_glaciers')
	studyRegionPoly <- sp::spTransform(studyRegionPoly, getCRS('wgs84'))
	
	studyRegionRast <- rasterize(studyRegionPoly, dalton[[1]])
	studyRegionRast <- studyRegionRast * 0 + 1
	
	daltonInterp <- daltonInterp * studyRegionRast
	daltonInterp <- trim(daltonInterp, padding=1)
	studyRegionRast <- trim(studyRegionRast, padding=1)

	### project raster
	##################
	
	daltonInterpEa <- projectRaster(daltonInterp, crs=getCRS('albersNA'))

	daltonInterpEa <- calc(daltonInterpEa, fun=function(x) ifelse(x > 1, 1, x))
	daltonInterpEa <- calc(daltonInterpEa, fun=function(x) ifelse(x < 0, 0, x))
	
	### mask out Great Lakes if they were not covered by ice
	########################################################

	usa <- getData('GADM', country='USA', level=2, path='C:/ecology/!Scratch')
	can <- getData('GADM', country='CAN', level=2, path='C:/ecology/!Scratch')

	nam2Sp <- rbind(can, usa)

	# get lakes for removal from other geographies... add small buffer because lake borders don't exactly align
	lakesSp <- nam2Sp[nam2Sp$ENGTYPE_2 == 'Water body', ]
	lakesSpEa <- sp::spTransform(lakesSp, getCRS('albersNA', TRUE))
	lakesSpEa <- gBuffer(lakesSpEa, width=10)
	lakesSpEa <- gUnaryUnion(lakesSpEa)
	
	daltonInterpEaLakesMasked <- daltonInterpEa
	
	# remove cell if ice is <1 and >50% or more of cell is covered by a lake
	for (countDate in seq_along(calYearsTo)) {
	
		calYearTo <- calYearsTo[countDate]

		this <- daltonInterpEaLakesMasked[[countDate]]

		lakesExtract <- extract(this, lakesSpEa, weights=TRUE, normalizeWeights=FALSE, cellnumber=TRUE)[[1]]
		lakesExtract <- as.data.frame(lakesExtract)

		# lake area --> NA when lake covers >50% of cell and ice < 1
		thisValues <- getValues(this)
		thisValues[lakesExtract$cell] <- ifelse(lakesExtract$value < 1 & lakesExtract$weight > 0.5, NA, lakesExtract$value)
	
		this <- setValues(this, thisValues)
	
		daltonInterpEaLakesMasked[[countDate]] <- this
	
	}
	
	### save!
	#########
	
	daltonInterpEa <- trim(daltonInterpEa, padding=1)
	daltonInterpEaLakesMasked <- trim(daltonInterpEaLakesMasked, padding=1)
	
	dirCreate('./!study_region_masks')
	names(daltonInterpEa) <- paste0('yr', abs(calYearsTo), 'bp')	
	names(daltonInterpEaLakesMasked) <- paste0('yr', abs(calYearsTo), 'bp')	

	writeRaster(daltonInterpEa, paste0('./!study_region_masks/study_region_dalton_ice_mask_no_lakes_', rastInterpFx, 'Interpolation'))
	writeRaster(daltonInterpEaLakesMasked, paste0('./!study_region_masks/study_region_dalton_ice_mask_lakes_masked_', rastInterpFx, 'Interpolation'))
	
	### plot
	########

	# political boundaries from GADM
	can <- getData('GADM', country='CAN', level=1, path='C:/ecology/!Scratch')
	usa <- getData('GADM', country='USA', level=1, path='C:/ecology/!Scratch')
	mex <- getData('GADM', country='MEX', level=1, path='C:/ecology/!Scratch')
	cub <- getData('GADM', country='CUB', level=1, path='C:/ecology/!Scratch')
	dom <- getData('GADM', country='DOM', level=1, path='C:/ecology/!Scratch')
	hti <- getData('GADM', country='HTI', level=1, path='C:/ecology/!Scratch')
	bah <- getData('GADM', country='BHS', level=1, path='C:/ecology/!Scratch')
	jam <- getData('GADM', country='JAM', level=1, path='C:/ecology/!Scratch')
	
	nam <- rbind(can, usa, mex, cub, dom, hti, bah)
	nam <- sp::spTransform(nam, getCRS('albersNA', TRUE))
	
	studyRegionRast <- projectRaster(studyRegionRast, crs=getCRS('albersNA'))
	studyRegionRast <- trim(studyRegionRast, padding=1)
	
	ext <- extent(studyRegionRast)
	ext <- as(ext, 'SpatialPolygons')
	projection(ext) <- getCRS('albersNA')
	
	breaks <- seq(0, 1, by=0.1)
	cols <- scales::alpha('cornflowerblue', seq(0, 1, length.out=length(breaks) - 1))

	# get dates for Dalton ice layers for overlap onto map
	daltonDates <- read.csv('C:/Ecology/Drive/Data/North American Ice Sheet Dalton et al 2020 Quaternary Science Reviews/Dalton et al 2020 QSR Dates from Shapefile Names.csv')
	
	daltonDates <- -1000 * daltonDates$calKiloYear
	
	dirCreate(paste0('./study_region_images/', rastInterpFx, 'Interpolation'))
	for (countYear in seq_along(calYearsTo)) {

		calYearTo <- calYearsTo[countYear]
		say('Plotting ', calYearTo, '...')

		png(paste0('./study_region_images/', rastInterpFx, 'Interpolation/study_region_', omnibus::prefix(21000 + calYearTo, 5), 'ybp_after_21000_ybp.png'), width=1600, height=1200, res=200)
		
			par(oma=c(0, 0, 0, 0), mar=rep(0.5, 4))
			plot(ext, ann=FALSE, border=NA, xpd=NA)
			plot(nam, border=NA, col='gray90', add=TRUE)
			thisStudyRegion <- studyRegionRast * daltonInterpEa[[countYear]] * 0 + 1
			plot(thisStudyRegion, col='gray70', ann=FALSE, legend=FALSE, add=TRUE)
			plot(nam, border='black', add=TRUE)
			plot(daltonInterpEa[[countYear]], legend=FALSE, col=cols, breaks=breaks, add=TRUE)
	
			entirelyIce <- calc(daltonInterpEa[[countYear]], fun=function(x) ifelse(x == 1, 1, NA))
			plot(entirelyIce, legend=FALSE, col='blue4', add=TRUE)
			
			# dalton shapefile overlay.. adding the two that are temporally closest to the given date
			daltonFrom <- daltonDates[calYearTo >= daltonDates]
			daltonFrom <- daltonFrom[length(daltonFrom)]
			daltonFrom <- sprintf('%02.2f', abs(daltonFrom / 1000))

			load(paste0('C:/Ecology/Drive/Data/North American Ice Sheet Dalton et al 2020 Quaternary Science Reviews/RDA Files/daltonEtAl2020_', daltonFrom, '_kiloCalYBP.rda'))
			daltonIceSpEa <- sp::spTransform(daltonIce, getCRS('albersNA', TRUE))
			plot(daltonIceSpEa, border='darkgoldenrod3', add=TRUE, lwd=2)

			daltonTo <- daltonDates[calYearTo <= daltonDates]
			
			if (length(daltonTo) == 0) {
				daltonTo <- NULL
			} else {
				daltonTo <- daltonTo[1]
				daltonTo <- sprintf('%02.2f', abs(daltonTo / 1000))
				load(paste0('C:/Ecology/Drive/Data/North American Ice Sheet Dalton et al 2020 Quaternary Science Reviews/RDA Files/daltonEtAl2020_', daltonTo, '_kiloCalYBP.rda'))
				daltonIceSpEa <- sp::spTransform(daltonIce, getCRS('albersNA', TRUE))
				plot(daltonIceSpEa, border='darkgoldenrod3', add=TRUE, lwd=2, lty='dotted')
			}
			
			legend('bottomright', inset=c(0, 0.22), bty='n', cex=1,
				title=paste(abs(calYearTo), 'YBP'),
				legend=c(
					'ice sheet = 100%',
					'ice sheet < 100%',
					'study region',
					'Period-starting ice sheet',
					'Period-ending ice sheet'
				),
				fill=c(
					'blue4',
					'cornflowerblue',
					'gray70',
					NA,
					NA
				),
				border=c(
					NA,
					NA,
					NA,
					NA,
					NA
				),
				col=c(
					NA,
					NA,
					NA,
					'darkgoldenrod3',
					'darkgoldenrod3'
				),
				lwd=c(
					NA,
					NA,
					NA,
					2,
					2
				),
				lty=c(
					NA,
					NA,
					NA,
					'solid',
					'dotted'
				)
			)
			
		dev.off()
		
	}
	
#################################	
say('DONE!!!', deco='%', level=1)
#################################
