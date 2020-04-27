# green_ash_study_region
The code and files in this repository are used to generate a study region for analysis of green ash biogeographic history. The code delineates an initial study region based on watershed basin boundaries, then is delineated using the process described in the code, but also requires a manual delineation to account for changing shoreline through time.

The code:
* Delineates a study region (polygon) based on watershed boundaries. The boundaries are chosen to create a contiguous area that encompasses the current and potential past distribution of green ash, relevant pollen deposits, and paleoshorelines, while representing a "reasonable" are (i.e., no geographic "choke points", no major islands, etc.).
* Expands this region using a manually-drawn polygon to account for shifting shorelines.
* Generates one raster representing exposed land mass for every 30 yr from 21 Kybp to 0 ybp.
* "Burns" values representing proportional ice cover in cells so they have either 0 (no ice), 1 (100% ice), some value between 0 and 1 (proportion of cell covered in ice), or `NA` (not land or outside study region).

The study region is represented by a raster stack in the folder `!study_region_masks`. This folder contains stack with/out the Great Lakes masked out (as `NA`s). There is one raster in each stack per 30 yr time period, with the first ("uppermost") raster representing land cover and ice 21 Kybp, and the last (number 701) present-day land cover.

Adam B. Smith
Last updated 2020-04-27
