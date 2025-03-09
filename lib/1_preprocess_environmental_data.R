##################################################################################
#
# Preprocessing of environmental & disturbance data layers for matching & analysis
#
#################################################################################


#---------- libraries

library(terra)
library(dplyr)

#---------- set working directory
setwd("~data/") 


#--------- load data

# note: we cannot provide the set aside forest area layer, as those partially include vector files provided under a data usage agreement
# from the regional authorities - for reproducibility, please request the respective layers at the regional authorities, while you can 
# openly access some here: 
# Bavaria: https://geoportal.bayern.de/geoportalbayern/seiten/dienste (search for "Naturwaldreservate")
# Saarland: https://geoportal.saarland.de/article/Download_Schutzgebietskataster/
# Sachsen-Anhalt, Hessen, Niedersachen, Schleswig-Holstein: contact this Group within "Nordwestliche Versuchsanstalt"  (https://www.nw-fva.de/wir/abteilungen/waldnaturschutz)

# check that set-aside forest areas cover at least 20 ha of forest and have been designated latest 1985!

setaside<- vect( "data/processed/setaside_forest_sites/all_setaside_forest_sites/setaside_forest.gpkg")
NPs <- vect("data/processed/setaside_forest_sites/NPs_Germany/NPs_Germany_25832.gpkg") # all national parks in Germany

sitecondition <- vect("data/raw/Wuchsbezirke/wb_2020.shp")
fcover<- rast("data/raw/disturbances/germany/forestcover_germany_25832.tiff")
dem <- rast("data/raw/dem/dem.tif")
slope <- rast("data/raw/dem/slope.tif")


#-------- preprocessing siteconditions

# merge all polygons of the same site condition together

sitecondition_merged <- aggregate(sitecondition, by ="wg_bu")
writeVector(sitecondition_merged, "data/processed/site_cond/siteconditions.gpkg")

# calculate average size per ecoregion (Wuchsgebiet)

sitecondition_merged$area_ha <- expanse(sitecondition_merged, unit = "ha")  # Calculate area in hectares
area_table <- as.data.frame(sitecondition_merged)
quantile(area_table$area_ha) # median 326,649 ha
mean(area_table$area_ha) # mean 435,876 ha

# create 30m res raster out of site conditions 

res <-  res(fcover)
extent <- ext(fcover) 
crs <- crs(fcover)
raster <- rast(extent=extent, resolution=res, crs=crs)
site.r <- rasterize(sitecondition_merged, raster, "wg_bu")

writeRaster(site.r, "data/processed/site_cond/siteconditions_raster.tif", overwrite=T)


#--------- preprocessing dem and slope

# resample to 30*30 m resolution

dem.30 <- resample(dem, fcover)
slope.30 <- resample(slope, fcover)

# reclassify values in DEM below 0 to 0

m <- c(-300, 0, 0)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
dem.30 <- classify(dem.30, rclmat, include.lowest=TRUE)


writeRaster(dem.30, "data/processed/dem/dem30.tif", overwrite=TRUE)
writeRaster(slope.30, "data/processed/dem/slope30.tif", overwrite=TRUE)



# -------- calculate Northerness (north-west exposition)

north.west <- cos((terra::terrain(dem.30, v = "aspect", unit="radians"))+45)
names(north.west) <- "north.west"

writeRaster(north.west, "data/processed/dem/north.west.tif", overwrite=TRUE)


#--------- preprocessing forest cover map

# reclassify forest cover map with all values >0.5 = forest 

m <- c(0, 0.5, NA, 0.5, 1, 1)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
fcover <- classify(fcover, rclmat, include.lowest=TRUE)

writeRaster(fcover, "data/processed/forestcover/forestcover.reclass.tif" )


# -- crop environmental variables to setaside and control forest cover area for matching in Germany

sitecondition <- rast("data/processed/site_cond/siteconditions_raster.tif")
fcover<- rast("data/processed/forestcover/forestcover.reclass.tif")
dem <- rast("data/processed/dem/dem30.tif")
slope <- rast("data/processed/dem/slope30.tif")
north.west <- rast("data/processed/dem/north.west.tif")

# --- mask environmental data to forest cover for comparability

dem.fcover <- crop(mask(dem, fcover), fcover)
slope.fcover <- crop(mask(slope, fcover), fcover)
ftype.fcover <- crop(mask(ftype, fcover), fcover)
sitecond.fcover <- crop(mask(sitecondition, fcover), fcover)
north.west.fcover <- crop(mask(north.west, fcover), fcover)


# -- crop to setaside forest sites

dem.setaside <- crop(mask(dem.fcover, setaside), setaside)
slope.setaside <- crop(mask(slope.fcover, setaside), setaside)
sitecond.setaside <- crop(mask(sitecond.fcover, setaside), setaside)
fcover.setaside <- crop(mask(fcover, setaside), setaside)
north.west.setaside <- crop(mask(north.west.fcover, setaside), setaside)

writeRaster(dem.setaside, "data/processed/dem/dem.setaside.tif", overwrite=T)
writeRaster(slope.setaside, "data/processed/dem/slope.setaside.tif", overwrite=T)
writeRaster(sitecond.setaside, "data/processed/site_cond/sitecond.setaside.tif", overwrite=T)
writeRaster(fcover.setaside, "data/processed/forestcover/fcover.setaside.tif", overwrite=T)
writeRaster(north.west.setaside, "data/processed/dem/north.west.setaside.tif", overwrite=T)

# -- crop to control areas (aka managed forests)

dem.control <- mask(dem.fcover, setaside, inverse=TRUE)
slope.control <- mask(slope.fcover, setaside, inverse=TRUE)
sitecond.control <- mask(sitecond.fcover, setaside, inverse=TRUE)
fcover.control <-  mask(fcover, setaside, inverse=TRUE)
north.west.control <-  mask(north.west.fcover, setaside, inverse=TRUE)

# mask out National Parks

dem.control <- mask(dem.control, NPs, inverse=TRUE)
slope.control <- mask(slope.control, NPs, inverse=TRUE)
sitecond.control <- mask(sitecond.control, NPs, inverse=TRUE)
fcover.control <-  mask(fcover.control, NPs, inverse=TRUE)
north.west.control <-  mask(north.west.control, NPs, inverse=TRUE)


writeRaster(dem.control, "data/processed/dem/dem.control.tif", overwrite=T)
writeRaster(slope.control, "data/processed/dem/slope.control.tif", overwrite=T)
writeRaster(sitecond.control, "data/processed/site_cond/sitecond.control.tif", overwrite=T)
writeRaster(fcover.control, "data/processed/forestcover/fcover.control.tif", overwrite=T)
writeRaster(north.west.control, "data/processed/dem/north.west.control.tif", overwrite=T)



# ---- disturbance information preprocessing ----

# --- dist agents

dist.agents <- rast("data/raw/disturbances/dist_agent/agent_classes_germany.tif")

# Define the projection string for CRS EPSG:3035
proj_string <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +datum=WGS84 +units=m +no_defs"

# Assign the projection to the raster
crs(dist.agents) <- proj_string

dist.agents <- project(dist.agents, severity)

dist.harvest <- dist.agents == 3 # harvest
dist.natural <- dist.agents == 1 # barkbeetle/Windthrow (2 = fire, only in Brandenburg, which is not included in our study)


#writeRaster(dist.harvest, "data/processed/disturbances/germany/disturbance_harvest.tif", overwrite=T)
writeRaster(dist.natural, "data/processed/disturbances/germany/disturbance_natural.tif", overwrite=T)

# - mask disturbance year map with agent map, to get dist year for natural disturbances only

# Load disturbance raster
disturbance.year.map <- terra::rast("data/raw/disturbances/germany/disturbance_year_germany_25832.tif")
disturbance.patches.natural <- terra::rast("data/processed/disturbances/germany/disturbance_natural.tif")

disturbance.year.map <- mask(disturbance.year.map, disturbance.patches.natural, maskvalues= 0)

writeRaster(disturbance.year.map, "data/processed/disturbances/germany/disturbance_natural_year.tif")


# --- severity

severity <- rast("data/raw/disturbances/germany/disturbance_severity_germany_25832.tif")

# Calculate high severity map

high.severity <- severity >= 0.8

# different threshold for sensitivity analysis - if you want to reproduce this analysis, you need to create the layers
#high.severity <- severity >= 0.7
#high.severity <- severity >= 0.9


names(high.severity) <- "high.severity"
writeRaster(high.severity, "data/processed/disturbances/germany/high.severity.tif")
#writeRaster(high.severity, "data/processed/disturbances/germany/high.severity_0.7.tif")
#writeRaster(high.severity, "data/processed/disturbances/germany/high.severity_0.9.tif")

high.severity.natural <- mask(high.severity, disturbance.patches.natural, maskvalues= 0)
writeRaster(high.severity.natural, "data/processed/disturbances/germany/high.severity.natural.tif")
#writeRaster(high.severity.natural, "data/processed/disturbances/germany/high.severity.natural_0.7.tif")
#writeRaster(high.severity.natural, "data/processed/disturbances/germany/high.severity.natural_0.9.tif")


