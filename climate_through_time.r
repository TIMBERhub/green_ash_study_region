### DEFINING STUDY REGION FOR GREEN ASH
### Adam B. Smith | Missouri Botanical Garden | 2019-12
###
### TO RUN:
### source('C:/Ecology/Drive/Research Active/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/study_region/climate_through_time.r')
### source('D:/Ecology/Drive/Research Active/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/study_region/climate_through_time.r')
### source('E:/Ecology/Drive/Research Active/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/study_region/climate_through_time.r')

### DESCRIPTION
### This script creates a plot of mean temperature through time using data from Lorenz, D.J., Nieto-Lugilde, D., Blois, J.L., Fitzpatrick, M.C., and Williams, J.W.  2016.  Downscaled and debiased climate simulations for North America from 21,000 years ago to 2100AD.  Scientific Data 3:160048.

### CODE SECTION CONTENTS ###
### setup ###
### compile table of mean annual temperatures through time ###

#############
### setup ###
#############

	memory.limit(memory.limit() * 2^30)
	rm(list=ls())
	gc()
	options(stringsAsFactors=FALSE)

	library(ggplot2)
	library(omnibus)
	library(raster)
	rasterOptions(format='GTiff', overwrite=TRUE)
	
	# drive where study region is to be created
	# workDrive <- 'C:'
	# workDrive <- 'D:'
	workDrive <- 'E:'
	
	# drive where external (others') data is stored
	# extDrive <- 'C:'
	# extDrive <- 'D:'
	extDrive <- 'E:'
	
	setwd(paste0(workDrive, '/ecology/Drive/Research Active/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/study_region'))

say('##############################################################')
say('### compile table of mean annual temperatures through time ###')
say('##############################################################')

	masks <- stack('./!study_region_raster_masks/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif')
	mask <- masks[[1]] * 0 + 1

	temperatures <- data.frame()
	for (gcm in c('ecbilt', 'ccsm')) {
	
		path <- if (gcm == 'ecbilt') {
			'E:/Ecology/Climate/Lorenz et al 2016 North America 21Kybp to 2100 CE/Version 2017-06-16/ecbilt_21-0k_all_tifs'
		} else if (gcm == 'ccsm') {
			'E:/Ecology/Climate/Lorenz et al 2016 North America 21Kybp to 2100 CE/Version 2017-06-16/ccsm3_22-0k_all_tifs'
		}
		
		for (ybp in seq(21000, 0, by=-500)) {
		
			say(gcm, ' ', ybp)
		
			tmax <- raster(paste0(path, '/', ybp, 'BP/an_avg_TMAX.tif'))
			tmin <- raster(paste0(path, '/', ybp, 'BP/an_avg_TMIN.tif'))
			
			temps <- stack(tmin, tmax)
			tmean <- mean(temps)
			
			tmean <- projectRaster(tmean, mask)
			tmean <- tmean * mask
			
			meanTemp <- cellStats(tmean, 'mean')
		
			temperatures <- rbind(
				temperatures,
				data.frame(
					gcm = gcm,
					ybp = ybp,
					meanAnnualTemp_C = meanTemp
				)
			)
		
		}
	}

	write.csv(temperatures, './Mean Annual Temperature in Study Region across Time.csv', row.names=FALSE)

	p <- ggplot(temperatures, aes(x=-1 * ybp, y=meanAnnualTemp_C, color=gcm)) +
		geom_line() +
		xlab('YBP') + ylab('Mean annual temperature (C)')
		
	print(p)
	

#################################	
say('DONE!!!', deco='%', level=1)
#################################
