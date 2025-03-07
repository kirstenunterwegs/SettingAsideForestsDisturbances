#######################################################################
#
# Aggregating disturbances in Germany to define pulse disturbance years
# and background disturbance years
#
######################################################################

# --- libaries

library(terra)
library(ggplot2)

# --- set working dir
setwd("~/data/") 


### ----- load disturbances initiated by natural disturbance agent -----


# load all (natural) disturbances for Germany

dist <- rast("data/processed/disturbances/germany/disturbance_natural_year.tif")

# calculate area disturbed per year

year_counts <- freq(dist)

# calculate baseline disturbance rate for the whole time period

year_counts$disturbance_ha <- year_counts$count*0.09

# calculate average distrubances over whole period

year_counts$mean_ref <- mean(year_counts$disturbance_ha) 

# calculate the divergence of annual area disturbed from the reference baseline 

year_counts$anomaly <- year_counts$disturbance_ha / year_counts$mean_ref -1 

year_counts$class <- ifelse(year_counts$anomaly > 1, "pulse", "background")


write.csv(year_counts, "data/processed/disturbances/germany/pulse_background_classification_natural.dist.csv")

My_Theme = theme( # for Poster plot
  axis.title.x = element_text(size = 50),
  axis.text.x = element_text(size = 40),
  axis.text.y = element_text(size = 40),
  axis.title.y = element_text(size = 50))

anomaly_plot <- ggplot2::ggplot(year_counts, ggplot2::aes(x=value, y=anomaly)) +
  geom_point(size = 6) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_hline(yintercept = 1, linetype = "dashed", col="red")+
  geom_text(aes(label = value), vjust = -0.5, size = 8)+
  theme_light() +
  labs(x = "Year", y = "Anomaly")+
  My_Theme

ggsave("data/results/plots/anomaly_years.tiff", anomaly_plot, width = 20, height = 15,limitsize = FALSE) # Adjust dimensions as needed


# get vector of years as reference
reference_period <- 1986:2015
drought_period   <- 2018:2020


# Calculate anomalies
year_counts2 <-
  year_counts %>%
  group_by(gridindex) %>%
  filter(sum(disturbance_ha) > 35) %>% # Exclude areas with less than 1 ha/yr of disturbances on average
  filter(sum(disturbance_ha[year %in% reference_period]) > 30) %>% # Exclude areas with less than 1 ha/yr of disturbances on average
  mutate(mean_ref    = mean(disturbance_ha[year %in% reference_period], na.rm = T),
         sum_18_20   = sum( disturbance_ha[year %in% drought_period], na.rm = T)/3, 
         sum_19_20   = sum( disturbance_ha[year %in% c(2019,2020)], na.rm = T)/2, 
         anomaly     = disturbance_ha / mean(disturbance_ha[year %in% reference_period], na.rm = TRUE) - 1, 
         anomaly_18_20  = sum_18_20 / mean_ref - 1,
         anomaly_19_20  = sum_19_20 / mean_ref - 1) %>% 
  ungroup()
