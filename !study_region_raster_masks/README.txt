 study_region_raster_masks
 Tue Apr 28 23:47:40 2020
 
 The files in this folder contain raster stacks that represent the exposed land a
 nd ice sheet cover of the study region for the green ash study. Values in the rasters are:
 * 0: Land.
 * 1: Completely ice.
 * between 0 and 1: Proportion of cell covered by ice.
 * NA: Either land that is outside the study region or water.
 
 The "topmost" layer (layer #1) in the stack represents the study region 21 Kybp, and the "bottommost" layer (layer #701) represents the study region at 0 ybp.
 
 Raster stacks represent either: a scenario with Great Lakes represented by NA ce
 lls if a cell was covered by more than 50% lake; or a scenario assuming the Lake
 s are entirely land. A few other inland cells are NA.  These are carryover from 
 the original rasters from Lorenz et al. 2016 Scientific data.
 
 The interpolation method refers to the manner in which the cover of ice in cells
  was "smoothed" from the intervals at which ice cover was provided by Dalton et 
 al. 2020 QSR to 30-yr intervals.
 
 The rasters are in Albers equal-area projection for North America.
