### DEFINING STUDY REGION FOR GREEN ASH
### Adam B. Smith | Missouri Botanical Garden | 2019-12
###
### TO RUN:
### source('C:/Ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/study_region/defining_study_region.r')
### source('D:/Ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/study_region/defining_study_region.r')
### source('E:/Ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/study_region/defining_study_region.r')

### DESCRIPTION
### This script delineates a geographic domain to be used for simulating the biogeographic history of green ash.  As per meetings on 2019-11-12 and 2019-11-19, we will demarcate the region as such:
### * West: the North American continental divide
### * South: southern extent of the state of Tamaulipas, Mexico
### * East: Atlantic ocean (extended outward to account for seal level rise since 21 Kybp)
### * North: Northward to the extent of the Hudson Bay drainage

### DATA SOURCES
### * Level I and Level II Watershed Boundaries from the USGS/Council on Environmental Cooperation (https://www.sciencebase.gov/catalog/item/4fb697b2e4b03ad19d64b47f)
### * GADM 3.6 political geography (https://gadm.org/)
### * Version 2 of Lorenz et al. climate layers to demarcate land versus sea since the LGM (Lorenz et al. 2016.  Downscaled and debiased climate simulations for North America from 21,000 years ago to 2100AD. Scientific Data 3:160048.)
### * Pollen data from Neotoma
### * Occurrence data for Fraxinus americana from BIEN Version 4.1
### * Locations of genetic data from Allan Strand sent to Adam in late 2019
### * Shapefiles representing ice sheet cover from ~21 Kybp to the present from Dalton et al. Quaternary Science Reviews 234:106223.
### * Elevation and bathymetry data from ETOP (https://data.nodc.noaa.gov/cgi-bin/iso?id=gov.noaa.ngdc.mgg.dem:316). Using the "grid-registered" version because it's the authoritative version.

### CODE SECTION CONTENTS ###
### setup ###
### compile pollen and occurrence data for Fraxinus ###
### create spatial polygon encompassing study extent ###
### mask Dalton et al 2020 ice sheet layers onto land mass from Lorenz et al 2016 ###
### generate elevation raster for study region ###

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
	library(BIEN)
	# library(ncdf4)
	library(spatialEco)
	
	library(omnibus) # Adam's custom library on GitHub: adamlilith/omnibus
	library(enmSdm) # Adam's custom library on GitHub: adamlilith/enmSdm
	library(birdsEye) # Adam's custom library on GitHub: adamlilith/birdsEye

	rasterOptions(format='GTiff', overwrite=TRUE)
	
	library(compiler)
	enableJIT(1)
	setCompilerOptions(suppressUndefined=TRUE)

	# drive where study region is to be created
	# workDrive <- 'C:'
	# workDrive <- 'D:'
	workDrive <- 'E:'
	
	# drive where external (others') data is stored
	# extDrive <- 'C:'
	# extDrive <- 'D:'
	extDrive <- 'E:'
	
	setwd(paste0(workDrive, '/ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/study_region'))

	# raster interpolation settings... see enmSdm::interpolateRasters()
	rastInterpFx <- 'linear' # decided on using this because splines create apparent "reversals" of ice melt
	
	# rastInterpFx <- 'spline'
	# splineMethod <- 'hyman'
	
# say('#######################################################')
# say('### compile pollen and occurrence data for Fraxinus ###')
# say('#######################################################')

	# # Using this data in ArcMAP to visualize candidate study region.

	# say('get pollen data', level=2)
	# ###############################
	
		# availFraxData <- get_dataset(loc = c(-145, 10, -50, 71),
			# datasettype = 'pollen',
			# taxonname = 'Fraxinus%'
		# )

		# fraxData <- get_download(availFraxData, verbose = FALSE)
		
		# compiledFrax <- compile_taxa(fraxData, 'P25')

		# pollen <- data.frame()
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
			
			# pollen <- rbind(
				# pollen,
				# data.frame(
					# longitude = compiledFrax[[i]]$dataset$site$long,
					# latitude = compiledFrax[[i]]$dataset$site$lat,
					# maxFraxinusProp = maxFraxProp
				# )
			# )
			
		# }
		
		# pollen <- SpatialPointsDataFrame(pollen [ , c('longitude', 'latitude')], data=pollen, proj4=getCRS('wgs84', TRUE))
		# dirCreate('./fraxinus_pollen_neotoma')
		# save(pollen, file='./fraxinus_pollen_neotoma/fraxinus_pollen_data_proportion_by_site.rda')
		
	# say('get occurrence data', level=2)
	# ###################################
	
		# # Download from BIEN 4.1

		# occsRaw <- BIEN_occurrence_species(
			# species = 'Fraxinus pennsylvanica',
			# cultivated = FALSE,
			# only.new.world = TRUE,
			# all.taxonomy = FALSE,
			# native.status = FALSE,
			# natives.only = TRUE,
			# observation.type = TRUE,
			# political.boundaries = FALSE,
			# collection.info = FALSE
		# )
		
		# ### remove occurrences with missing coordinates and dates:
		# occs <- occsRaw[!(is.na(occsRaw$longitude) | is.na(occsRaw$latitude) | is.na(occsRaw$date_collected)), ]
		# dim(occs)

		# ### remove records <1950 (period covered by Lorenz et al. 2016 climate data is 1950-2005):
		# occs$date_collected <- as.Date(occs$date_collected)
		# occs <- occs[!is.na(occs$date_collected), ]
		# occs <- occs[which(occs$date_collected >= as.Date(paste0('1950-01-01'))), ]
		
		# dirCreate('./occurrences_bien')
		# save(occs, file='./occurrences_bien/fraxinus_pennsylvanica_bien_all_occurrences.rda')
	
# say('########################################################')
# say('### create spatial polygon encompassing study extent ###')
# say('########################################################')

	# say('This portion of the code uses watershed boundaries to delineate a study region extent that covers the area of interest for projecting range dynamics of green ash from 21 Kybp to the present in eastern North America. The watersheds mainly encompass the portions east of the Rocky Mountain continental divide, but sub-basins are included/excluded to create a contiguous region without choke-points that encompasses the present distribution, relevant pollen deposits, and potential past refugia.', breaks=80, pre=1)

	# say('We will use watershed basin boundaries to define the study region. Basins are taken from shapefiles from the Commission on Environmental Cooperation (http://www.cec.org/tools-and-resources/map-files/watersheds).', breaks=80, pre=1)
	
	# # watershed boundaries from CEC
	# load(paste0(extDrive, '/Ecology/Watersheds/CEC/watersheds_level_1.rda'))
	# load(paste0(extDrive, '/Ecology/Watersheds/CEC/watersheds_level_2_and_3.rda'))

	# # political boundaries from GADM
	# can <- getData('GADM', country='CAN', level=1, path=paste0(extDrive, '/ecology/!Scratch'))
	# usa <- getData('GADM', country='USA', level=1, path=paste0(extDrive, '/ecology/!Scratch'))
	# mex <- getData('GADM', country='MEX', level=1, path=paste0(extDrive, '/ecology/!Scratch'))
	
	# nam <- rbind(can, usa, mex)
	# nam <- sp::spTransform(nam, CRS(projection(ws1)))
	
	# # clean watershed names
	# ws1$NAW1_EN <- trim(ws1$NAW1_EN)
	# ws1$NAW2_EN <- trim(ws1$NAW2_EN)
	# ws1$NAW3_EN <- trim(ws1$NAW3_EN)
	# ws1$NAW4_EN <- trim(ws1$NAW4_EN)

	# # get major drainage areas
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

	# # combine watershed-delineated study region with manually-drawn region that encompasses adjacent littoral area that was exposed during LGM
	# littoral <- shapefile('./formerly_exposed_land_from_lorenz_et_al_2016/manually_drawn_study_region_encompassing_continental_shelf')
	# demesne <- gUnaryUnion(demesne)
	
	# demesne <- gUnion(demesne, littoral)
	# demesne <- spatialEco::remove.holes(demesne)

	# ### plot study region with spatial data on green ash
	# ####################################################
		
		# # load species data
		# load('./occurrences_bien/fraxinus_pennsylvanica_bien_all_occurrences.rda')
		# load('./fraxinus_pollen_neotoma/fraxinus_pollen_data_proportion_by_site.rda')
		# genetics <- shapefile(paste0(workDrive, '/Ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/genetic_data_locations_allan_strand_2019_11_12'))
		# dirCreate('./genetic_data')
		# save(genetics, file='./genetic_data/genetic_data_locations_2019_11_12')
		# shapefile(genetics, './genetic_data/genetic_data_locations_2019_11_12', overwrite=TRUE)
		
		# # convert data to spatial format
		# occs <- sp::SpatialPointsDataFrame(occs[ , c('longitude', 'latitude')], data=occs, proj4string=getCRS('wgs84', TRUE))
		
		# # range maps
		# bien <- shapefile(paste0(workDrive, '/Ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/enms/regions/bien_range_map/Fraxinus_pennsylvanica'))
		# little <- shapefile(paste0(workDrive, '/Ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/enms/regions/littles_range_map/fraxpenn'))
		
		# # project data
		# occs <- sp::spTransform(occs, CRS(projection(demesne)))
		# pollen <- sp::spTransform(pollen, CRS(projection(demesne)))
		# genetics <- sp::spTransform(genetics, CRS(projection(demesne)))

		# # project range maps
		# bien <- sp::spTransform(bien, CRS(projection(demesne)))
		# little <- sp::spTransform(little, CRS(projection(demesne)))
		
		# dirCreate('./study_region_images')
		# png('./study_region_images/initial_study_region_by_watershed_no_glaciers.png', width=1600, height=1200, res=300)
		
			# par(oma=c(0, 0, 0, 0), mar=rep(0.5, 4))
			# plot(demesne, col='darkolivegreen3', border='darkolivegreen')
			# plot(nam, border='black', add=TRUE)
			
			# plot(bien, lwd=2, border='darkblue', add=TRUE)
			# plot(little, lwd=2, border='darkgreen', add=TRUE)
			
			# points(occs, pch=16, cex=0.25, col='darkorange')
			
			# points(pollen, pch=2, cex=0.6)
			# points(pollen[pollen$maxFraxinusProp >= 0.05, ], pch=24, cex=0.6, bg='yellow')
			
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
		
	# dirCreate('./study_region_spatial_polygons')
	# demesneAlb <- sp::spTransform(demesne, getCRS('albersNA'))
	# shapefile(demesneAlb, './study_region_spatial_polygons/study_region_mask_without_glaciers', overwrite=TRUE)

# say('#####################################################################################')
# say('### mask Dalton et al 2020 ice sheet layers onto land mass from Lorenz et al 2016 ###')
# say('#####################################################################################')	
	
	# say('This step will create rasters of the land mass of North America through time using climate layers from Lorenz et al 2016 Sci Data. We will overlay the ice sheet layers from Dalton et al 2020 QSR onto the layers to yield cell values ranging from 0 (no ice) to 1 (complete ice cover), with values between 0 and 1 representing proportion of the cell covered in ice. We will create one raster per time period in the time series provided by the Dalton ice cover data product.', post=2, breaks=80)
	
	# say('Note, we have to use either the climate layers from the CCSM or ECBilt global circulation models, as these are the only two projected back through time.', breaks=80)

	# ### get rasters representing land from Lorenz et al

		# lorenzDir <- paste0(extDrive, '/Ecology/Climate/Lorenz et al 2016 North America 21Kybp to 2100 CE/Version 2017-06-16/ccsm3_22-0k_all_tifs/')
		
		# years <- seq(22000, 0, by=-500)

		# if (exists('lorenz')) rm(lorenz)
		# for (year in years) {
			# lorenz <- if (exists('lorenz')) {
				# stack(lorenz, raster(paste0(lorenzDir, '/', year, 'BP/an_avg_ETR.tif')))
			# } else {
				# raster(paste0(lorenzDir, '/', year, 'BP/an_avg_ETR.tif'))
			# }
		# }
		
		# lorenz <- lorenz * 0 + 1
		# lorenzYears <- seq(22000, 0, by=-500)
		# names(lorenz) <- paste0('yr', lorenzYears, 'bp')
	
	# ### get carbon and calendar years for each ice layer for Dalton et al ice sheet data
	# # note that the calendar years are not exactly the same as listed on Table 1 in Dalton et al 2020
		
		# daltonYears <- read.csv(paste0(workDrive, '/Ecology/Drive/Data/North American Ice Sheet Dalton et al 2020 Quaternary Science Reviews/Dalton et al 2020 QSR Dates from Shapefile Names.csv'))

	# ### interpolate land rasters to time periods matching Durant et al ice sheet representations

		# interpTo <- -1000 * daltonYears$calKiloYear
		# lorenzInterp <- interpolateRasters(lorenz, interpFrom=-1 * lorenzYears, interpTo=interpTo, type=rastInterpFx)
		
		# names(lorenzInterp) <- paste0('yr', 1000 * daltonYears$calKiloYear, 'bp')
	
	# ### for each cell, assign a value from 0 to 1 indicating proportion covered by ice sheet
	
		# daltonDir <- paste0(workDrive, '/Ecology/Drive/Data/North American Ice Sheet Dalton et al 2020 Quaternary Science Reviews/RDA Files/')
	
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
	
		# year0 <- lorenz[[nlayers(lorenz)]] * 0
		# lorenzInterp <- stack(lorenzInterp, year0)

		# dirCreate('./ice_sheet/glaciers_dalton')
		# writeRaster(lorenzInterp, paste0('./ice_sheet/glaciers_dalton/daltonGlaciers_on_lorenzLand'))
		# file.copy(paste0(workDrive, '/Ecology/Drive/Data/North American Ice Sheet Dalton et al 2020 Quaternary Science Reviews/Dalton et al 2020 QSR Dates from Shapefile Names.csv'), to=paste0(workDrive, '/Ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/study_region/ice_sheet/glaciers_dalton/Dalton et al 2020 QSR Dates from Shapefile Names.csv'))
		
		# sink(paste0('./ice_sheet/glaciers_dalton/README_', rastInterpFx, 'Interpolation.txt'), split=TRUE)
		
			# say('The raster stack in this folder represents land area used by Lorenz et al 2016 Sci Data overlaid with ice sheet coverage from Dalton et al 2020 Quat Sci Reviews. Values range from 0 (no portion of the cell covered by ice) to 1 (all of the cell covered by ice). Approximate CALENDAR year represented by each layer in the stack is provided by the layer name AND in the table "Dalton et al 2020 QSR Dates from Shapefile Names.csv". The ', rastInterpFx, ' interpolation method was used.', breaks=80)
			
		# sink()

# say('#######################################################', pre=2)
# say('### generate study region raster stack through time ###')
# say('#######################################################', post=2)

	# say('This section interpolates the land/ice rasters from the time periods provided by Dalton et al. to 30-yr time steps.', post=2, breaks=80)
	
	# say('In initial runs there was a problem with the overlaying of Dalton ice polygons onto rasters in cells adjacent to water bodies, especially for coastal cells along the northeast section of the study region. Since raster cells are rectangular, these cells occasionally had <100% ice cover even though the ice sheet polygon extended to the coast (and maybe slightly beyond). So, this script forces all cells that 1) are adjacent to an NA cell and have an ice cover of >0.5 to have an ice cover value of 1 (100% ice).', breaks=80, post=2)

	# ### interpolate Dalton ice sheet layer to every 30 yr
	# #####################################################

	# dalton <- stack(paste0('./ice_sheet/glaciers_dalton/daltonGlaciers_on_lorenzLand.tif'))
	# daltonDates <- read.csv('./ice_sheet/glaciers_dalton/Dalton et al 2020 QSR Dates from Shapefile Names.csv')

	# calYearsFrom <- c(-1 * daltonDates$calKiloYear * 1000, 0)
	# calYearsTo <- seq(-21000, 0, by=30)
	
	# names(dalton) <- paste0('yr', abs(calYearsFrom), 'bp')
	
	# daltonInterp <- interpolateRasters(dalton, interpFrom=calYearsFrom, interpTo=calYearsTo, type=rastInterpFx)
	# names(daltonInterp) <- paste0('yr', abs(calYearsTo), 'bp')

	# ### mask with study region polygon created above
	# ################################################
	
	# studyRegionPoly <- shapefile('./study_region_spatial_polygons/study_region_mask_without_glaciers')
	# studyRegionPoly <- sp::spTransform(studyRegionPoly, getCRS('wgs84'))
	
	# studyRegionRast <- rasterize(studyRegionPoly, dalton[[1]])
	# studyRegionRast <- studyRegionRast * 0 + 1
	
	# daltonInterp <- daltonInterp * studyRegionRast
	# daltonInterp <- trim(daltonInterp, padding=1)
	# studyRegionRast <- trim(studyRegionRast, padding=1)

	# ### project raster
	# ##################
	
	# daltonInterpEa <- projectRaster(daltonInterp, crs=getCRS('albersNA'))

	# daltonInterpEa <- calc(daltonInterpEa, fun=function(x) ifelse(x > 1, 1, x))
	# daltonInterpEa <- calc(daltonInterpEa, fun=function(x) ifelse(x < 0, 0, x))
	
	# ### mask out Great Lakes if they were not covered by ice
	# ########################################################

	# usa <- getData('GADM', country='USA', level=2, path='C:/ecology/!Scratch')
	# can <- getData('GADM', country='CAN', level=2, path='C:/ecology/!Scratch')

	# nam2Sp <- rbind(can, usa)

	# # get lakes for removal from other geographies... add small buffer because lake borders don't exactly align
	# lakesSp <- nam2Sp[nam2Sp$ENGTYPE_2 == 'Water body', ]
	# lakesSpEa <- sp::spTransform(lakesSp, getCRS('albersNA', TRUE))
	# lakesSpEa <- gBuffer(lakesSpEa, width=10)
	# lakesSpEa <- gUnaryUnion(lakesSpEa)
	
	# daltonInterpEaLakesMasked <- daltonInterpEa
	
	# # remove cell if ice is <1 and >50% or more of cell is covered by a lake
	# for (countDate in seq_along(calYearsTo)) {
	
		# this <- daltonInterpEaLakesMasked[[countDate]]

		# lakesExtract <- extract(this, lakesSpEa, weights=TRUE, normalizeWeights=FALSE, cellnumber=TRUE)[[1]]
		# lakesExtract <- as.data.frame(lakesExtract)

		# # lake area --> NA when lake covers >50% of cell and ice < 1
		# thisValues <- getValues(this)
		# thisValues[lakesExtract$cell] <- ifelse(lakesExtract$value < 1 & lakesExtract$weight > 0.5, NA, lakesExtract$value)
	
		# this <- setValues(this, thisValues)

		# daltonInterpEaLakesMasked[[countDate]] <- this
	
	# }
	
	# ### force cells adjacent to water with >x% of ice cover to have complete ice cover
	# ##################################################################################

	# forceIce <- function(x) if (any(is.na(x) & x[5] %>na% 0.5)) { 1 } else { x[5] }
	# w <- matrix(1, nc=3, nr=3)
	# for (i in 1:nlayers(daltonInterpEaLakesMasked)) {
	# for (i in 1) {
		# daltonInterpEa[[i]] <- focal(daltonInterpEa[[i]], w=w, fun=forceIce)
		# daltonInterpEaLakesMasked[[i]] <- focal(daltonInterpEaLakesMasked[[i]], w=w, fun=forceIce)
	# }

	# ### force cells complete surrounded by complete ice to be complete ice
	# ######################################################################

	# iceThreshold <- 0.99
	# for (i in 1:nlayers(daltonInterpEaLakesMasked)) {
		# daltonInterpEa[[i]] <- calc(daltonInterpEa[[i]], function(x) ifelse(x > iceThreshold, 1, x))
		# daltonInterpEaLakesMasked[[i]] <- calc(daltonInterpEaLakesMasked[[i]], function(x) ifelse(x > iceThreshold, 1, x))
	# }

	# ### save!
	# #########
	
	# daltonInterpEa <- trim(daltonInterpEa, padding=1)
	# daltonInterpEaLakesMasked <- trim(daltonInterpEaLakesMasked, padding=1)
	
	# dirCreate('./!study_region_raster_masks')
	# names(daltonInterpEa) <- paste0('yr', abs(calYearsTo), 'bp')	
	# names(daltonInterpEaLakesMasked) <- paste0('yr', abs(calYearsTo), 'bp')	

	# writeRaster(daltonInterpEa, paste0('./!study_region_raster_masks/study_region_daltonIceMask_noLakes_', rastInterpFx, 'IceSheetInterpolation'))
	# writeRaster(daltonInterpEaLakesMasked, paste0('./!study_region_raster_masks/study_region_daltonIceMask_lakesMasked_', rastInterpFx, 'IceSheetInterpolation'))
	
	# sink('./!study_region_raster_masks/README.txt', split=TRUE)
	
		# say('study_region_raster_masks')
		# say(date())
		# say('')
		# say('The files in this folder contain raster stacks that represent the exposed land and ice sheet cover of the study region for the green ash study. Values in the rasters are:', breaks=80)
		# say('* 0: Land.')
		# say('* 1: Completely ice.')
		# say('* between 0 and 1: Proportion of cell covered by ice.')
		# say('* NA: Either land that is outside the study region or water.')
		# say('')
		# say('The "topmost" layer (layer #1) in the stack represents the study region 21 Kybp, and the "bottommost" layer (layer #701) represents the study region at 0 ybp.')
		# say('')
		# say('Raster stacks represent either: a scenario with Great Lakes represented by NA cells if a cell was covered by more than 50% lake; or a scenario assuming the Lakes are entirely land. A few other inland cells are NA.  These are carryover from the original rasters from Lorenz et al. 2016 Scientific data.', breaks=80)
		# say('')
		# say('')
		# say('Cells that are adjacent to an NA cell and that have >50% ice cover have been forced to have 100% ice cover. Cells that are not NA and completely surrounded by ice are forced to be 100% ice.', breaks=80)
		# say('')
		# say('')
		# say('The interpolation method refers to the manner in which the cover of ice in cells was "smoothed" from the intervals at which ice cover was provided by Dalton et al. 2020 QSR to 30-yr intervals.', breaks=80)
		# say('')
		# say('The rasters are in Albers equal-area projection for North America.')
	
	# sink()

	# ### plot
	# ########

	# # political boundaries from GADM
	# can <- getData('GADM', country='CAN', level=1, path='C:/ecology/!Scratch')
	# usa <- getData('GADM', country='USA', level=1, path='C:/ecology/!Scratch')
	# mex <- getData('GADM', country='MEX', level=1, path='C:/ecology/!Scratch')
	# cub <- getData('GADM', country='CUB', level=1, path='C:/ecology/!Scratch')
	# dom <- getData('GADM', country='DOM', level=1, path='C:/ecology/!Scratch')
	# hti <- getData('GADM', country='HTI', level=1, path='C:/ecology/!Scratch')
	# bah <- getData('GADM', country='BHS', level=1, path='C:/ecology/!Scratch')
	# jam <- getData('GADM', country='JAM', level=1, path='C:/ecology/!Scratch')
	
	# nam <- rbind(can, usa, mex, cub, dom, hti, bah)
	# nam <- sp::spTransform(nam, getCRS('albersNA', TRUE))
	
	# studyRegionRast <- projectRaster(studyRegionRast, crs=getCRS('albersNA'))
	# studyRegionRast <- trim(studyRegionRast, padding=1)
	
	# ext <- extent(studyRegionRast)
	# ext <- as(ext, 'SpatialPolygons')
	# projection(ext) <- getCRS('albersNA')
	
	# breaks <- seq(0, 1, by=0.1)
	# cols <- scales::alpha('cornflowerblue', seq(0, 1, length.out=length(breaks) - 1))

	# ### plot one image per generation (30 yr)
	# #########################################
	
	# daltonDates <- -1000 * daltonDates$calKiloYear
	
	# dirCreate(paste0('./study_region_images/', rastInterpFx, '_ice_sheet_interpolation'))
	# for (countYear in seq_along(calYearsTo)) {

		# calYearTo <- calYearsTo[countYear]
		# say('Plotting ', calYearTo, '...')

		# png(paste0('./study_region_images/', rastInterpFx, '_ice_sheet_interpolation/study_region_', omnibus::prefix(21000 + calYearTo, 5), 'ybp_after_21000_ybp.png'), width=1000, height=750, res=100)
		
			# par(oma=c(0, 0, 0, 0), mar=rep(0.5, 4))
			# plot(ext, ann=FALSE, border=NA, xpd=NA)
			# plot(nam, border=NA, col='gray90', add=TRUE)
			# thisStudyRegion <- studyRegionRast * daltonInterpEa[[countYear]] * 0 + 1
			# plot(thisStudyRegion, col='gray70', ann=FALSE, legend=FALSE, add=TRUE)
			# plot(nam, border='black', add=TRUE)
			# plot(daltonInterpEa[[countYear]], legend=FALSE, col=cols, breaks=breaks, add=TRUE)
	
			# entirelyIce <- calc(daltonInterpEa[[countYear]], fun=function(x) ifelse(x == 1, 1, NA))
			# plot(entirelyIce, legend=FALSE, col='blue4', add=TRUE)
			
			# # dalton shapefile overlay... adding the two that are temporally closest to the given date
			# daltonFrom <- daltonDates[calYearTo >= daltonDates]
			# daltonFrom <- daltonFrom[length(daltonFrom)]
			# daltonFrom <- sprintf('%02.2f', abs(daltonFrom / 1000))

			# load(paste0(extDrive, '/Ecology/Drive/Data/North American Ice Sheet Dalton et al 2020 Quaternary Science Reviews/RDA Files/daltonEtAl2020_', daltonFrom, '_kiloCalYBP.rda'))
			# daltonIceSpEa <- sp::spTransform(daltonIce, getCRS('albersNA', TRUE))
			# # plot(daltonIceSpEa, border='darkgoldenrod3', add=TRUE, lwd=2)
			# plot(daltonIceSpEa, border='cyan', add=TRUE, lwd=2)

			# daltonTo <- daltonDates[calYearTo <= daltonDates]
			
			# if (length(daltonTo) == 0) {
				# daltonTo <- NULL
			# } else {
				# daltonTo <- daltonTo[1]
				# daltonTo <- sprintf('%02.2f', abs(daltonTo / 1000))
				# load(paste0(extDrive, '/Ecology/Drive/Data/North American Ice Sheet Dalton et al 2020 Quaternary Science Reviews/RDA Files/daltonEtAl2020_', daltonTo, '_kiloCalYBP.rda'))
				# daltonIceSpEa <- sp::spTransform(daltonIce, getCRS('albersNA', TRUE))
				# # plot(daltonIceSpEa, border='darkgoldenrod3', add=TRUE, lwd=2, lty='dotted')
				# plot(daltonIceSpEa, border='cyan', add=TRUE, lwd=2, lty='dotted')
			# }
			
			# legend('bottomright', inset=c(0, 0.28), bty='n', cex=1,
				# title=paste(abs(calYearTo), 'YBP'),
				# legend=c(
					# 'ice sheet = 100%',
					# 'ice sheet < 100%',
					# 'study region',
					# 'Period-starting ice sheet',
					# 'Period-ending ice sheet'
				# ),
				# fill=c(
					# 'blue4',
					# 'cornflowerblue',
					# 'gray70',
					# NA,
					# NA
				# ),
				# border=c(
					# NA,
					# NA,
					# NA,
					# NA,
					# NA
				# ),
				# col=c(
					# NA,
					# NA,
					# NA,
					# 'cyan',
					# 'cyan'
				# ),
				# lwd=c(
					# NA,
					# NA,
					# NA,
					# 2,
					# 2
				# ),
				# lty=c(
					# NA,
					# NA,
					# NA,
					# 'solid',
					# 'dotted'
				# )
			# )
			
		# dev.off()
		
	# }

# say('##################################################')
# say('### generate elevation raster for study region ###')
# say('##################################################')

	# # elevation data from ETOP: https://data.nodc.noaa.gov/cgi-bin/iso?id=gov.noaa.ngdc.mgg.dem:316
	# etop <- raster('E:/Ecology/Topography/ETOP/grid_registered/ETOPO1_Bed_g_geotiff.tif')
	# projection(etop) <- getCRS('wgs84')
	
	# # study region polygon
	# demesneAlb <- shapefile('./study_region_spatial_polygons/study_region_mask_without_glaciers')
	# demesneWgs84 <- sp::spTransform(demesneAlb, getCRS('wgs84', TRUE))
	
	# studyRegionRasts <- brick('./!study_region_raster_masks/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif')
	
	# # resample elevation to study region resolution and projection
	# etop <- crop(etop, demesneWgs84)
	# etopAlb <- projectRaster(etop, studyRegionRasts)
	
	# writeRaster(etopAlb, './!study_region_raster_masks/study_region_elevationInMeters_fromEtop')
	
#################################	
say('DONE!!!', deco='%', level=1)
#################################
