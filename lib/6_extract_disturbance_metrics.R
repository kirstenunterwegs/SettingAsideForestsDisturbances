##################################################################################
#
# Extract disturbance information for sublandscapes
#
#################################################################################


#---------- libraries

library(dplyr)
library(terra)
library(tidyr)

#---------- working directory

setwd("~/data/")



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# --- create setaside disturbances stack for disturbance extraction ---


# Load disturbance rasters

disturbance.year.map <- terra::rast("data/processed/disturbances/germany/disturbance_natural_year.tif") # natural disturbances only
high.severity.map <- terra::rast("data/processed/disturbances/germany/high.severity.natural.tif") # natural disturbances only
forest.cover <- terra::rast("data/processed/forestcover/forestcover.reclass.tif")


# import research shapes
setaside <- vect( "data/processed/setaside_forest_sites/all_setaside_forest_sites/setaside_forest.gpkg")

# rasterize reserves and NPs (set aside forests)

res <-  res(forest.cover)
extent <- ext(forest.cover)
crs <- crs(forest.cover)
raster <- rast(extent=extent, resolution=res, crs=crs)

site_raster <- rasterize(setaside,raster ,field="newID")

# stack disturbance information and crop to control (managed) and research (setaside) area

disturbance_stack <- c(forest.cover, disturbance.year.map, high.severity.map, site_raster)
names(disturbance_stack) <- c("forest.cover", "disturbance.year", "high.severity","ID")

# Create Research Stack

Research_disturbance_stack <- crop(mask(disturbance_stack, setaside), setaside)

# Create Control Stack

Control_disturbance_stack <- terra::mask(disturbance_stack, setaside, inverse = TRUE)

rm(disturbance_stack, forest.cover, disturbance.year.map, high.severity.map, site_raster) # clear RAM

writeRaster(Research_disturbance_stack, "data/processed/disturbances/dist_stack_managed_setaside/research_stack_nat.agent.tif")
writeRaster(Control_disturbance_stack, "data/processed/disturbances/dist_stack_managed_setaside/control_stack_nat.agent.tif")



# ----- assign matching design for extraction:

folder <- "1o1_forest" # change folder to respective matching setup, for which I want to extract disturbances

# 1o1_forest
# 1oMany_forest

# ----- 

# load disturbance information stacks

#Research_disturbance_stack <- rast("data/processed/disturbances/dist_stack_managed_unmanaged/research_stack_nat.agent.tif")
#Control_disturbance_stack <- rast("data/processed/disturbances/dist_stack_managed_unmanaged/control_stack_nat.agent.tif")


# load patch id information for disturbance information extraction


stack_id_reserach <- rast(paste0("data/processed/disturbances/patch_size/",folder,"/patch_n_stack_reserach.tif"))
names(stack_id_reserach) <- paste0("patch.id.",1986:2020)
stack_id_control <- rast(paste0("data/processed/disturbances/patch_size/",folder,"/patch_n_stack_control.tif"))
names(stack_id_control) <- paste0("patch.id.",1986:2020)


# stack all disturbance information for extraction

Research_disturbance_stack <- c(Research_disturbance_stack, stack_id_reserach)
Control_disturbance_stack <- c(Control_disturbance_stack, stack_id_control)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### --- function to extract disturbance_information for sublandscape
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


get_disturbance_metrics <- function(x, y, Research, window){
 
  # define focal window/ extent to extract disturbance information depending on set aside (Research =1) or managed (Reserach = 0) forests
  
   if (Research == 1){
    temp_stack <- Research_disturbance_stack
    # Define Neighborhood matrix
    if (window %in% c(15,17,19,21,23,25,27,29,31)){
      neighborhood_matrix <- matrix(1,nrow=window,ncol=window)
    } 
    else{
      neighborhood_matrix <- matrix(1, nrow = window+1, ncol = window+1)
      neighborhood_matrix[window+1, ] <- neighborhood_matrix[, window+1] <- 0
    }
  } 
  
   if (Research == 0){
    temp_stack <- Control_disturbance_stack
    # Define Neighborhood matrix
    neighborhood_matrix <- matrix(1, nrow = window, ncol = window)
    
    # Modify the neighborhood_matrix if the window size is not odd sized
    if (window == 16) {
      neighborhood_matrix <- matrix(1, nrow = window+1, ncol = window+1)
      neighborhood_matrix[window+1, ] <- neighborhood_matrix[, window+1] <- 0 
    }
  }
  
  # identify cell IDs to extract disturbance information from
  
  cell  <- terra::cellFromXY(temp_stack, cbind(x, y)) #Testcell: x <- 4604600; y <- 2855400
  
  cellnumbers <- as.vector(terra::adjacent(temp_stack, cell, directions = neighborhood_matrix))
  disturbance_information <- as.data.frame(temp_stack[cellnumbers, drop = FALSE], xy=TRUE)
  
  cells_in_sublandscape <- nrow(dplyr::filter(disturbance_information, !is.na(forest.cover))) #ID / forest.cover
  
  
  #get numbers of disturbed cells per year
  
  cells.per.year <- disturbance_information %>%
    dplyr::filter(complete.cases(disturbance.year))%>%
    dplyr::group_by(disturbance.year) %>% 
    dplyr::summarise(ncells.per.year = n()) %>%
    dplyr::right_join(data.frame(disturbance.year = 1986:2020,
                                 default.cells.per.year = 0), by = c("disturbance.year")) %>% 
    dplyr::mutate(ncells.per.year = pmax(ncells.per.year, default.cells.per.year, na.rm=TRUE)) %>%
    dplyr::arrange(disturbance.year) %>%
    dplyr::select(-c(default.cells.per.year)) %>% 
    tidyr::spread(key=disturbance.year, value = ncells.per.year) %>%
    as.matrix()
  

  #get number of severe disturbed cells per year
  
  severity.cells.per.year <- disturbance_information %>%
    dplyr::filter(high.severity == 1) %>%
    dplyr::group_by(disturbance.year) %>%
    dplyr::summarise(ncells.per.year = n()) %>%
    dplyr::right_join(data.frame(disturbance.year = 1986:2020,
                                 default.cells.per.year = 0), by = c("disturbance.year")) %>%
    dplyr::mutate(ncells.per.year = pmax(ncells.per.year, default.cells.per.year, na.rm=TRUE)) %>%
    dplyr::arrange(disturbance.year) %>%
    dplyr::select(-c(default.cells.per.year)) %>%
    tidyr::spread(key=disturbance.year, value = ncells.per.year) %>%
    as.matrix()

  
  # get disturbance patch sizes (patches extend only within landscape! - crop disturbance patches to landscapes)
  
  patch.size <- disturbance_information %>%
    dplyr::filter(complete.cases(disturbance.year))%>%
    pivot_longer(cols = starts_with("patch.id."),
                 names_to = "year",
                 values_to = "patch_id") %>%
    filter(!is.na(patch_id)) %>% # Remove NA patch IDs
    count(year, patch_id) %>% # Count occurrences of each patch ID per year
    mutate(area_ha = n * 0.09) %>% # Convert counts to hectares
    select(-n)
  
  # maximum patch size in landscape
  max.patch <- patch.size %>%
    summarise(max = if (all(is.na(area_ha))) { 0 } 
              else { max(area_ha, na.rm = TRUE)}) %>%
    as.matrix()
  
  # (area-) weighted mean patch size (following Li & Archer 1997)
  wm_patch_size <- patch.size %>%
    count(area_ha) %>%
    mutate(n_s = n * area_ha) %>%
    mutate(n_s2 = n * area_ha^2) %>%
    summarise( sum_n_s2 = sum(n_s2),
               sum_n_s = sum(n_s)) %>%
    mutate(wmps = if (sum_n_s == 0) 0 else sum_n_s2 / sum_n_s) %>%
    select(wmps) %>% 
    as.matrix()
  
  # get density (n of disturbances patches per landscape)
    # Counting the number of unique values in patch.id. columns for each year
  
  density <- disturbance_information %>%
    select(starts_with("patch.id.")) %>%
    summarise(across(everything(), ~length(unique(na.omit(.))))) %>%
    sum(.)%>%
    as.matrix()
  
  # get frequency (n of disturbance events per landscape)
    # Counting the number of unique values in freq columns for each year
  
  frequency <- disturbance_information %>%
    select(disturbance.year) %>%
    filter(disturbance.year > 0) %>%
    summarise(n_events = length(unique(disturbance.year))) %>%
    as.matrix()
  


   return(list("ncells" = cells_in_sublandscape,
               "cells.per.year" = cells.per.year,
               "severity.cells.per.year" = severity.cells.per.year,
               "max.patch" = max.patch,
               "mean.patch" = wm_patch_size,
               "density" = density,
               "frequency" = frequency))
}

options(dplyr.summarise.inform = FALSE)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### --------------------- disturbance extraction -------------------- #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# load the respective matching data frame

if(folder == "1o1_forest") {
  
  export_df <- readRDS("data/processed/matching/match.df.1o1_15perc.rds")
  
} else if (folder == "1oMany_forest") {
  
  export_df <- readRDS("data/processed/matching/match.df.1oMany_15perc_n30_filtered.rds")
  
}


# --- assign the window size and forest cover of the research site to attain window for extracting control site disturbance info


export_df <- export_df %>%
  group_by(subclass) %>%
  mutate(fcover = fcover[which(Research == 1)[1]]) # to deal with pseudo replication in 1:Many sampling design of set aside sublandscapes

# decide on neighboorhood-matrix size depending on forestcover in Reserach site, to guarantee max. 2 ha difference
export_df$window <- ifelse(export_df$Research == 0 & export_df$fcover > 242 & export_df$fcover < 270, 16,
                                                ifelse(export_df$Research == 0 & export_df$fcover > 270, 17,
                                                  ifelse(export_df$Research == 0 & export_df$fcover <= 242, 15, export_df$window)))
                           


# extraction:

pb <- dplyr::progress_estimated(nrow(export_df))
for (i in 1:nrow(export_df)){
  sublandscape <- export_df[i,]
 disturbance_informations <- get_disturbance_metrics(x = sublandscape$x, y = sublandscape$y, Research = sublandscape$Research, window = sublandscape$window)
  if ((disturbance_informations$ncells - sublandscape$fcover) > 22 | (disturbance_informations$ncells - sublandscape$fcover) < -22) {print(paste("more than 2 ha difference in forest cover in ", i))} # 22*0.09 = 2.07 ha
  
  export_df[i, paste0("disturbed.cells.", 1986:2020)] <- disturbance_informations$cells.per.year 
  export_df[i, paste0("severe.disturbed.cells.", 1986:2020)] <- disturbance_informations$severity.cells.per.year
  export_df[i, "density"] <- disturbance_informations$density
  export_df[i, "max.patch"] <- disturbance_informations$max.patch
  export_df[i, "mean.patch"] <- disturbance_informations$mean.patch
  export_df[i, "frequency"] <- disturbance_informations$frequency
  export_df[i, "n.fcover"] <- disturbance_informations$ncells

  pb$tick()$print()
  
}

# round weighted mean patch size
export_df$mean.patch <- round(export_df$mean.patch, 2)

head(export_df)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# add climatic data to check whether site conditions properly reflect temperature and precipitation similarities
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# ----- climate data extraction ----

# unzip climate data (downloaded multi-annual averages from https://opendata.dwd.de/climate_environment/)


# precipitation

# ras_path <- "raw/climate/grids_germany_multi_annual_precipitation_1991-2020_17.asc.gz"
# 
# R.utils::gunzip(ras_path, remove = FALSE)
# 
# # average temperature
# 
# ras_path <- "raw/climate/grids_germany_multi_annual_air_temp_mean_1991_2020_17.asc.gz"
# 
# R.utils::gunzip(ras_path, remove = FALSE)

# 
# prec <- rast("raw/climate/grids_germany_multi_annual_precipitation_1991-2020_17.asc")
# temp <- rast("raw/climate/grids_germany_multi_annual_air_temp_mean_1991_2020_17.asc")

# # assign crs from documentation of dwd data
# crs(prec) <- "EPSG:31467"
# crs(temp) <- "EPSG:31467"
# 
# # reproject raster data for climate data extraction
# 
# prec <- project(prec, Research_disturbance_stack)
# temp <- project(temp, Research_disturbance_stack)
# 
# writeRaster(prec,"raw/climate/multi_annual_precipitation_1991-2020_17_25832.tiff" )
# writeRaster(temp,"raw/climate/multi_annual_air_temp_mean_1991_2020_17_25832.tiff" )


prec <- rast("data/raw/climate/multi_annual_precipitation_1991-2020_17_25832.tiff")
temp <- rast("data/raw/climate/multi_annual_air_temp_mean_1991_2020_17_25832.tiff" )

# create points from sublandscape df for raster extraction

subland.points <- vect(as.data.frame(export_df), geom=c("x", "y"), crs ="EPSG:25832" ) 

# extract climate metrics

subland.points <- terra::extract(prec, subland.points, bind=TRUE, method= "simple")
subland.points <- terra::extract(temp, subland.points, bind=TRUE, method= "simple")

dist.climate.df <- as.data.frame(subland.points, geom="XY")

# rename climate metrics

dist.climate.df <- dist.climate.df %>%
  rename(mean_prec91_20 = "grids_germany_multi_annual_precipitation_1991-2020_17",
         mean_temp91_20 = "grids_germany_multi_annual_air_temp_mean_1991_2020_17")


saveRDS(dist.climate.df, paste0("data/processed/dist.extraction/",folder,"/match.df_climate_dist.rds"))
