##################################################################################
#
# Get patch ID for disturbances in our sublandscape to extract patch size
#
#################################################################################


#---------- libraries

library(dplyr)
library(terra)
library(tidyr)
library(raster)

#---------- working directory
setwd("~/data/") 


# -----

folder <- "1o1" # change folder to respective matching setup, for which I want to extract disturbances

# 1o1_forest
# 1oMany_forest

# ----- 

# --- load data 


# !!! decide for the right matching setup !!!

# 1:1 Match
sublandscapes <- as.data.frame(readRDS("data/processed/matching/match.df.1o1_15perc.rds"))

# 1:Many Match
#sublandscapes <- as.data.frame(readRDS("data/processed/matching/match.df.1oMany_15perc_n30_filtered.rds"))


# load natural disturbances
dist <- rast("data/processed/disturbances/germany/disturbance_natural_year.tif")


sublandscapes <- vect(sublandscapes, geom=c("x", "y"), crs="EPSG:25832") # vectorize sublandscape centroid

# create buffer around sub landscapes
buffer_sublandscapes <- buffer(sublandscapes, width = 2000)

# mask disturbances with focused area around sublandscapes
dist.sub <- mask(dist, buffer_sublandscapes)
writeRaster(dist.sub, paste0("data/processed/disturbances/patch_size/", folder,"/nat_disturbances_year_in_buffer.tif"))

rm(dist, dist.sub, buffer_sublandscapes, sublandscapes) # clear RAM


# load subsetted natural disturbance map
dist.sub <- rast(paste0("data/processed/disturbances/patch_size/", folder,"/nat_disturbances_year_in_buffer.tif"))

# clump annual disturbance patches within buffer

years <- 1986:2020 


for (i in 1:length(years)){ # loop takes a while!
  
  year <- years[i]
  print(paste("----------------------", year, "---------------------"))
  temp.rast <- dist.sub 
  temp.rast[temp.rast != year] <- NA
  
  
  # natural Disturbance patches
  natural.patch.map <- terra::patches(terra::mask(dist.sub, temp.rast), directions=4, zeroAsNA=TRUE) # rook-contiguity
  
  # write to file
  writeRaster(natural.patch.map, paste0("data/processed/disturbances/patch_size/", folder,"/years_n/patch_n_", year,".tiff"))
  
  # calculate patch size 
  natural.patch.matrix <- freq(natural.patch.map)
  
  # Create a reclassification matrix
  reclass_matrix <- natural.patch.matrix %>%
    mutate(count = count * 0.09) %>%  # Convert count to hectares
    dplyr::select(value, count) %>%   
    as.matrix()
  
  
  # Reclassify the raster; assign patch size as cell value
  natural.patch.map_reclass <- classify(natural.patch.map, rcl = reclass_matrix) 
  
  
  # write to file
  writeRaster(natural.patch.map_reclass, paste0("data/processed/disturbances/patch_size/", folder,"/years_patch/patch_size_", year,".tiff"))
  
  rm(temp.rast, natural.patch.map, natural.patch.map_reclass) # free RAM at the end of each iteration
  
}

rm(list = ls()) # empty RAM

#% !!!!!!!!!!!!!!
# Use the patch ID information in the extraction, to calculate the size of the patch within the sublandscape
# You can use the patch_size layers in the extraction process to extract the complete patch size, also when it exceeds the sublandscape

# We decided to stick to the patch size WITHIN the landscape (so using the patch_n_ layers)
#% !!!!!!!!!!!!!!

# -----

folder <- "1o1_forest" # change folder to respective matching setup

# 1o1_forest
# 1oMany_forest

# ----- 


# load annual layers and mask them to managed and set aside forest area

# # Get a list of all .rds files in the folder
# raster_files <- list.files(paste0("data/processed/disturbances/patch_size/", folder,"/years_patch/"), all.files =TRUE, full.names = TRUE)
# # Remove the first two paths
# raster_files <- raster_files[-c(1, 2)]
# 
# # Read all rasters and stack them
# stack <- lapply(raster_files, terra::rast)
# 
# # Create a stack from raster list
# patch_size_stack <- rast(stack)
# 

# Get a list of all .rds files in the folder
raster_files <- list.files(paste0("data/processed/disturbances/patch_size/", folder,"/years_n/"), all.files =TRUE, full.names = TRUE)
# Remove the first two paths
raster_files <- raster_files[-c(1, 2)]

# Read all rasters and stack them
stack <- lapply(raster_files, terra::rast)

# Create a stack from raster list
patch_n_stack <- rast(stack)



# Rename layers in the stack for clarity

years <- 1986:2020
names(patch_size_stack) <- as.character(years)
#writeRaster(patch_size_stack,"processed/disturbances/patch_size/1oMany_forest/patch_size_stack.tif")

names(patch_n_stack) <- as.character(years)
#writeRaster(patch_n_stack,"processed/disturbances/patch_size/1oMany_forest/patch_n_stack.tif")


# separate patch_size_stack layer into set aside/reserve and managed

# setaside <- vect( "data/processed/setaside_forest_sites/all_setaside_forest_sites/setaside_forest.gpkg")
# 
# Research_patch_size_stack <- crop(mask(patch_size_stack, setaside), setaside)
# Control_patch_size_stack <- terra::mask(patch_size_stack, setaside, inverse = TRUE)
# 
# writeRaster(Research_patch_size_stack, paste0("data/processed/disturbances/patch_size/", folder,"/patch_size_stack_reserach.tif"))
# writeRaster(Control_patch_size_stack, paste0("data/processed/disturbances/patch_size/", folder,"/patch_size_stack_control.tif"))

Research_patch_n_stack <- crop(mask(patch_n_stack, setaside), setaside)
Control_patch_n_stack <- terra::mask(patch_n_stack, setaside, inverse = TRUE)

writeRaster(Research_patch_n_stack, paste0("data/processed/disturbances/patch_size/", folder,"/patch_n_stack_reserach.tif"))
writeRaster(Control_patch_n_stack, paste0("data/processed/disturbances/patch_size/", folder,"/patch_n_stack_control.tif"))


rm(list = ls()) # empty RAM

