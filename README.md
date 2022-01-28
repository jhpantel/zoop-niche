# zoop-niche
Data and analysis of niche occupancy and niche shifts in rural and urban environments for freshwater zooplankton.

# Citations
Pantel, J.H., J.M.T. Engelen, and L. De Meester. 2022. Niche use and co-occurrence patterns of zooplankton along a strong urbanization gradient. Ecography. In press.

# Table of Contents
The primary code file is main.R. This R file will upload all of the data, conduct all analyses included in the manuscript, and produce raw versions of all of the figures.

*data_files* contains all of the data files associated with this project. In some instances, analyses take a very long time (hours) or conduct randomizatio tests that may lead to p-values for significance tests that slightly differ from what appears in the manuscript. The code main.R thus will call an external .RData files (large_vars.RData) that contains these variables. All of the code and data needed to create what appears in large_vars.RData is incuded as well, and the commented code indicates that code. A saved version of the R environment that results when main.R is completed is included as 28_01_2022.RData.

- |-main.R
- |data_files
- | |-geographical.csv
- | |-spp-10-7-2015.csv
- | |-ZP-pa.csv
- | |-SPEEDY_R_continu2.csv
- | |-GBG.csv
- | |-BVA.csv
- | |-zoop_info.csv
- | |-I_D_test.R
- | |-large_vars.RData
- | |-Cscores.txt
- | |-Geomean.R
- | |-spat_interpol.csv
- | |-split_test_train.R
- | |-spat_interpolate.R
- | |-rural_urban_test.R
- | |-BUA.csv
- | |-28_01_2022.RData

- | |gadm36_BEL_shp
- | | |-gadm36_BEL_2.shp

- |raw_output
- | |Fig_3

There are 2 folders contained in *data_files*
*gadm36_BEL_shp* is a folder with shape / GIS files needed to create the map that is Figure S1. This was downloaded from https://www.gadm.org.

*raw_output* is an empty folder set that will hold the figures created in main.R.