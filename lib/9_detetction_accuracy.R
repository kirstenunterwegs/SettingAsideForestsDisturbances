###################################################################################
#
# Identifying detection accuracy for disturbance map in managed & set-aside forests
# 
# (1) create stratification raster for sample point
# (2) create confusion matrix
#
##################################################################################

#---------- libraries

library(terra)
library(dplyr)
library(caret)
library(ggplot2)
library(reshape2)

#---------- working directory

setwd("~/data/")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# (1) create stratification raster for sample point
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# load data

landscapes <- vect("data/processed/landscapes/1o1_forest/landscapes.gpkg")
crs(landscapes) <- "EPSG:25832"

set_aside_forest <- vect( "data/processed/natural_forest_sites/all_natural_forest_sites/natural_forest.gpkg")

fcover<- rast("data/processed/forestcover/forestcover.reclass.tif")
dist.year<- rast("data/raw/disturbances/germany/disturbance_year_germany_25832.tif")
agent <- rast("data/raw/disturbances/dist_agent/agent_classes_germany_25832.tif")


# --- subset to those landscape pairs, where at least one was disturbed

dist.year[agent == 3] <- NA # subset disturbances to only natural disturbances

# Extract number of disturbed pixel within each landscape

n_dist_pixel <- extract(dist.year, landscapes, fun = function(x) sum(!is.na(x))) # does not work

landscapes$has_disturbance <- n_dist_pixel[, 2] 

landscapes_sub <- landscapes[landscapes$subclass %in% unique(landscapes$subclass[
  landscapes$has_disturbance >= 1]), ]

writeVector(landscapes_sub, "data/processed/detection_accuracy/landscapes_disturbed.gpkg", overwrite=T)


# --- create stratification raster

# create an empty raster with same crs, ext and res als dist.year
strata_raster <- rast(ext(dist.year), res = res(dist.year), crs = crs(dist.year))

# set all pixels within the landscapes_sub vector to 1
strata_raster <- rasterize(landscapes_sub, strata_raster, field = 1)

# Mask out harvest disturbances 
strata_raster[agent == 3] <- NA

# Mask out disturbances pre-year 2000 (as no images to verify)
strata_raster[dist.year <= 2000] <- NA

# mask to only forested areas 
strata_raster[is.na(fcover)] <- NA

# mask all landscapes with Research == 1 with the set-aside_forest layer 

  # subset managed landscapes for later merging
  managed_areas <- landscapes_sub[landscapes_sub$Research == 0, ]
  strata_raster_managed <- mask(strata_raster, managed_areas)

  # create masking layer for set-aside forests
  set_aside_areas <- landscapes_sub[landscapes_sub$Research == 1, ]
  set_aside_mask <- rasterize(set_aside_areas, strata_raster, field = 1)
  set_aside_mask <- mask(set_aside_mask, set_aside_forest) # mask to  forest in landscape extent
  
  strata_raster[is.na(set_aside_mask)] <- NA

  strata_raster<- cover(strata_raster, strata_raster_managed)


# set disturbed pixel post year 2000 to 2 (disturbed stratum)
  
  strata_raster[strata_raster == 1 & dist.year > 2000] <- 2

  strata_raster <- rast("data/processed/detection_accuracy/strata_raster.tif")
  
# now relabel the cell values to 3 and 4 if 1 in managed or 2 in set-aside; than create stratified point layer

landscapes_sub_raster <-  rasterize(landscapes_sub, strata_raster, field = "Research")
set_aside_areas <- landscapes_sub[landscapes_sub$Research == 1, ]
strata_raster_set_aside <- mask(strata_raster, set_aside_forest)


# Reclassify where Research == 0
reclass_matrix <- matrix(c(1, 3,   # 1 -> 3
                           2, 4),  # 2 -> 4
                         ncol = 2, byrow = TRUE)

strata_raster_final <- classify(strata_raster, reclass_matrix, others = NA, include.lowest = TRUE) # reclassify disturbed and forest strata


# use old strata label in set-aside forests and use new classified strata label for managed forests 
# 1 = intact forest in set-aside; 2 = disturbed forest in set-aside
# 3 = intact forest in managed; 4 = disturbed forest in managed

strata_raster_final <- ifel(!is.na(strata_raster_set_aside), strata_raster_set_aside, strata_raster_final)

writeRaster(strata_raster_final, "data/processed/detection_accuracy/strata_raster_final.tif" )
#strata_raster_final <- rast("data/processed/detection_accuracy/strata_raster_final.tif")


# -- convert stratification raster to df for sample draw

raster_points <- as.data.frame(strata_raster_final, xy = TRUE, na.rm = TRUE)
saveRDS(raster_points, "data/processed/detection_accuracy/strata_points.rds")
#raster_points <- readRDS("data/processed/detection_accuracy/strata_points.rds")

# draw samples 

sampled_points <- raster_points %>%
  group_by(layer) %>%  # 'layer' is the column with raster class values
  slice_sample(n = 110, replace = FALSE) %>%  # Sample 125 points per class, 500 points in total 
  ungroup()

sampled_points_vect <- vect(sampled_points, geom = c("x", "y"), crs = crs(strata_raster_final))

# create a random ID for labeling, so we do not know which point should represent which stratum:

n <- nrow(as.data.frame(sampled_points_vect))
random_ids <- sample(61:500, n, replace = FALSE)
sampled_points_vect$fid <- random_ids

# change id order and store new order, to mix up points

sorted_indices <- order(sampled_points_vect$fid)
sampled_points_vect <- sampled_points_vect[sorted_indices, ]

# create label layers and remove strata layer info 
sampled_points_vect$disturbed <- NA
sampled_points_vect$quality <- NA
sampled_points_vect$layer  <- NULL


writeVector(sampled_points_vect, "data/processed/detection_accuracy/sampled_points.shp", overwrite=T) 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# load points in QGis, delete the stratum column to not distort labeling and safe as ESRI 
# shapefile to load into Google Earth Pro
#
# label points as distrurbed (1 = yes, 0 = no) and add quality flag (1 = good, 2 = unclear, 3 = no clear images)
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# load labeled data back into work space:
sampled_points <- vect("processed/detection_accuracy/sampled_points.shp")

# extract stratum

stratum <- extract(strata_raster_final, sample_points) 

sample_points$stratum <- stratum[, 2] 

sample_points.df <- as.data.frame(sample_points)

# relabel the disturbed labels to compare with original strata

sample_points.df <- sample_points.df %>%
  mutate(manual_label = case_when(
    stratum %in% c(1, 2) & disturbed == 0 ~ 1,
    stratum %in% c(1, 2) & disturbed == 1 ~ 2,
    stratum %in% c(3, 4) & disturbed == 0 ~ 3,
    stratum %in% c(3, 4) & disturbed == 1 ~ 4,
    TRUE ~ NA_real_
  ))

sample_points.df_clean <- sample_points.df %>% filter(quality != 3) # 493 samples left
sample_points.df_clean <- sample_points.df %>% filter(quality == 1) # 433 samples left

sample_points.df_clean %>% group_by(stratum) %>% count()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# (2) create confusion matrix
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# ---- confusion matrix

# choose sample setup
sample_points <- sample_points.df
sample_points <- sample_points.df_clean

# calculate the confusion matrix

conf_matrix <- confusionMatrix(
  factor(sample_points$stratum, levels = unique(sample_points$stratum)),# predicted (disturbance map)
  factor(sample_points$manual_label, levels = unique(sample_points$stratum))# "true" results 
)

print(conf_matrix)


# Calculate standard error
standard_error <- sd(as.numeric(sample_points$stratum) - as.numeric(sample_points$manual_label)) /  sqrt(nrow(sample_points))


# Create a confusion matrix table
conf_table <- conf_matrix$table
conf_df <- as.data.frame(as.table(conf_table))


# Relabel using recode
conf_df <- conf_df %>%
  mutate(
    Reference = recode(Reference,
                       `1` = "Intact - set-aside",
                       `2` = "Disturbed - set-aside",
                       `3` = "Intact - managed",
                       `4` = "Disturbed - managed"),
    Prediction = recode(Prediction,
                        `1` = "Intact - set-aside",
                        `2` = "Disturbed - set-aside",
                        `3` = "Intact - managed",
                        `4` = "Disturbed - managed"))%>%
  mutate(Freq = ifelse(Freq == 0, NA, Freq))


# Plot confusion matrix as a heatmap

My_Theme = theme(
  axis.title.x = element_text(size = 30),
  axis.text.x = element_text(size = 30,angle = 45, hjust = 1), #
  axis.text.y = element_text(size = 25),
  axis.title.y = element_text(size = 30),
  legend.title = element_text(size=30),
  legend.text = element_text(size=25))


conf_matrix_plot <- ggplot(conf_df, aes(Reference , Prediction , fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), color = "black", size = 10) +  # Add text labels
  scale_fill_gradient(low = "white", high = "blue", na.value = "white") +
  labs(x = "Reference", y = "Predicted", fill = "Frequency") +
  theme_bw() +
  My_Theme
conf_matrix_plot

ggsave("data/results/plots/confusion_matrix.tiff", conf_matrix_plot, width = 15, height = 10)


# percentage mislabeled: 

conf_df <-  conf_df %>%
  mutate(Freq = ifelse(is.na(Freq), 0, Freq))

comission_error_set_aside <- round(sum(conf_df$Freq[conf_df$Reference == "Intact - set-aside" & conf_df$Prediction == "Disturbed - set-aside"]) /
  sum(conf_df$Freq[conf_df$Reference == "Intact - set-aside"]) * 100,2)

omission_error_set_aside <- round(sum(conf_df$Freq[conf_df$Reference == "Disturbed - set-aside" & conf_df$Prediction == "Intact - set-aside"]) /
  sum(conf_df$Freq[conf_df$Reference == "Disturbed - set-aside"]) * 100,2)

comission_error_managed <- round(sum(conf_df$Freq[conf_df$Reference == "Intact - managed" & conf_df$Prediction == "Disturbed - managed"]) /
  sum(conf_df$Freq[conf_df$Reference == "Intact - managed"]) * 100,2)

omission_error_managed <- round(sum(conf_df$Freq[conf_df$Reference == "Disturbed - managed" & conf_df$Prediction == "Intact - managed"]) /
  sum(conf_df$Freq[conf_df$Reference == "Disturbed - managed"]) * 100,2)

# Print results
cat("Comission error set-aside:", comission_error_set_aside, "%\n")
cat("Omission error set-aside:", omission_error_set_aside, "%\n")
cat("Comission error managed:", comission_error_managed, "%\n")
cat("oomission error managed:", omission_error_managed, "%\n")
