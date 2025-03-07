##################################
#
# Compare natural potential vegetation 
# in Reserves & parks versus managed forests
#
##################################


#---------- libraries

library(terra)
library(dplyr)
library(ggplot2)
library(MetBrewer)

#---------- working directory

setwd("~/data/")

# ---- load the data:

npv <- vect("data/raw/NPV/Maps/npv_germany.gpkg")
npv_code <- read.csv("data/raw/NPV/npv_key_germany.csv", header=TRUE, sep=";")

fcover<- rast("data/processed/forestcover/forestcover.reclass.tif")
setaside_forest <- vect("data/processed/setaside_forest_sites/all_setaside_forest_sites/setaside_forest.gpkg")
nps <-  vect("data/processed/setaside_forest_sites/NPs_Germany/NPs_Germany_25832.gpkg")
wuchsgebiete <- vect("data/processed/site_cond/siteconditions_represented.gpkg")

forest_type <- rast("data/processed/forest_type/reclassification/forest.type.8.germany.30_covercrop.tif")
forest_type_key <- read.csv("data/processed/forest_type/reclassification/ftype_key.csv", header=TRUE, sep=";")
  
# rasterize setaside sites and National Parks for masking

res <-  res(fcover)
extent <- ext(fcover) 
crs <- crs(fcover)
raster <- rast(extent=extent, resolution=res, crs=crs)

setaside_forest.r <- rasterize(setaside_forest, raster, "ID")
nps.r <- rasterize(nps, raster, "WDPA_PID")

# rename npv code

npv <- merge(npv, npv_code[, c("code", "cover_aggregated")], by.x = "CODE", by.y = "code", all.x = TRUE)
npv.df <- as.data.frame(npv)


# mask NPV map to forest cover

npv.r <- rasterize(npv, raster, "CODE") # rasterize npv map
npv.forest <- mask(npv.r, fcover) # mask npv to forest cover

writeRaster(npv.forest, "data/processed/npv/npv.r_forest.tiff")
  

# extract NPV for setaside forests

npv.setaside <- mask(npv.forest, setaside_forest)
writeRaster(npv.setaside, "data/processed/npv/npv.setaside.tiff")
  

# extract NPV for managed forests 

# mask to those Wuchsgebiete represented in this study!

npv.managed <- mask(npv.forest, wuchsgebiete)
npv.managed <- mask(npv.managed, setaside_forest.r, inverse=TRUE) 
npv.managed <- mask(npv.managed, nps.r, inverse=TRUE) # mask out National Parks 

#writeRaster(npv.managed, "data/processed/npv/npv.managed.tiff")
writeRaster(npv.managed, "data/processed/npv/npv.managed_sub.tiff") # those masked to siteconditions represented

# compare NPV shares for managed and setaside forests

npv.setaside.df <- as.data.frame(npv.setaside)
npv.managed.df <- as.data.frame(npv.managed)

npv.setaside.summary <- npv.setaside.df %>% count(CODE)
npv.setaside.summary$type <- as.factor("setaside")
npv.setaside.summary <- npv.setaside.summary %>% mutate(n_pixel = sum(n)) # sum of pixels covered by setaside forests

npv.managed.summary <- npv.managed.df %>% count(CODE)
npv.managed.summary$type <- as.factor("Managed")
npv.managed.summary <- npv.managed.summary %>% mutate(n_pixel = sum(n))  # sum of pixels covered by managed forests

# recode vegetation type
npv.summary <- rbind(npv.setaside.summary, npv.managed.summary)
npv.summary <- merge(npv.summary, npv_code[, c("code", "cover_aggregated")], by.x = "CODE", by.y = "code", all.x = TRUE)

# shares according the vegetation type
npv.summary <- npv.summary %>% group_by(type, cover_aggregated) %>% 
                     mutate(n_cover = sum(n), # calculate amount of pixels in this NPV category per management
                            share = round(n_cover/n_pixel,2)) %>%  # calculate share per NPV category on area covered by management type
                            slice(1) # only one row per management - cover combination
  
  # Identify cover_aggregated groups with all share values equal to 0
  cover_to_remove <- npv.summary %>%
  group_by(cover_aggregated) %>%
  filter(all(share == 0)) %>%
  distinct(cover_aggregated) %>%
  ungroup()

  # Filter out rows from npv.summary that belong to these groups
  npv_summary_filtered <- npv.summary %>%
    anti_join(cover_to_remove, by = "cover_aggregated")

  
  
# set plotting theme  
My_Theme = theme(
  axis.title.x = element_text(size = 35),
  axis.text.x = element_text(size = 30, angle = 45,hjust=1),
  axis.text.y = element_text(size = 30),
  axis.title.y = element_text(size = 35),
  legend.key.height = unit(1, 'cm'),
  legend.key.width = unit(2, 'cm'),
  legend.title = element_text(size=35),
  legend.text = element_text(size=35),
  strip.text.x = element_text(size = 20),
  panel.spacing = unit(2, "lines"),
  legend.position = c(0.8, 0.8))

plot_npv <- ggplot(npv_summary_filtered, aes(x=cover_aggregated , y=share, fill=type)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  scale_fill_manual(values=met.brewer("Kandinsky"))+
  labs(x = "Potential setaside vegetation", y = "Share on area covered per management", fill = "") +
  theme_minimal()+
  My_Theme
plot_npv

ggsave("data/results/plots/npv_managed_setaside_sub.tiff", plot_npv, width = 20, height = 15)  #create theme first for better vis


# ----- compare PNV with realized forest type vegetation -----

ftype_setaside <- mask(forest_type, setaside_forest)

ftype_stack <- c(ftype_setaside, npv.setaside) 
names(ftype_stack) <- c("realized_ftype", "npv")
ftype_comparison_df <- as.data.frame(ftype_stack)

# recode forest types

ftype_comparison_df <- ftype_comparison_df %>% merge(npv_code[, c("code", "cover_aggregated")], by.x = "npv", by.y = "code", all.x = TRUE)
ftype_comparison_df <- ftype_comparison_df %>% merge(forest_type_key, by.x = "realized_ftype", by.y = "ID", all.x = TRUE)

write.csv(ftype_comparison_df,"data/processed/npv/npv_forest_type_comparison.csv")

sum(is.na(ftype_comparison_df)) # 62
ftype_comparison_df <- na.omit(ftype_comparison_df)

My_Theme = theme(
  axis.title.x = element_text(size = 35),
  axis.text.x = element_text(size = 30, angle = 45,hjust=1),
  axis.text.y = element_text(size = 30),
  axis.title.y = element_text(size = 35),
  legend.key.height = unit(1, 'cm'),
  legend.key.width = unit(2, 'cm'),
  legend.title = element_text(size=35),
  legend.text = element_text(size=30),
  strip.text.x = element_text(size = 20),
  panel.spacing = unit(2, "lines"))


# ggplot(ftype_comparison_df, aes(x = forest_type, fill = cover_aggregated)) +
#   geom_bar(position = "stack") +
#   labs(x = "Realized forest Type",
#        y = "Count",
#        fill = "PNV") +
#   scale_fill_manual(values=met.brewer("Renoir"))+
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))

pnv_ftype_count <- ggplot(ftype_comparison_df, aes(x = cover_aggregated, fill = forest_type)) +
  geom_bar(position = "stack") +
  labs(x = "PNV",
       y = "Count",
       fill = "Realized forest Type") +
  scale_fill_manual(values=met.brewer("Tiepolo"))+
  theme_minimal() +
  My_Theme


ftype_comparison_df_perc <- ftype_comparison_df %>%
  group_by(cover_aggregated, forest_type) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(cover_aggregated) %>%
  mutate(percentage = count / sum(count) * 100)


# Create the percentage stacked bar plot
pnv_ftype_perc <- ggplot(ftype_comparison_df_perc, aes(x = cover_aggregated, y = percentage, fill = forest_type)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "PNV",
       y = "Percentage",
       fill = "Realized forest Type") +
  scale_fill_manual(values=met.brewer("Tiepolo"))+
  theme_minimal()+
 My_Theme

ggsave("data/results/plots/pnv_forest_type_comparison_perc.tiff", pnv_ftype_perc, width = 20, height = 15) 
ggsave("data/results/plots/pnv_forest_type_comparison_count.tiff", pnv_ftype_count, width = 20, height = 15) 
