##################################################################################
#
# Matching landscapes
# Finding for the selected set aside landscape the most similar managed forest landscape
#
#(1) Create the search data frame which goes into the matchit function
#(2) Find matching landscapes (1:1 and 1:Many)
#(3) Filter in 1:Many matches the distance between multiple landscape matches to   
#    ensure minimal overlap
#(4) check on similarity in environmental features between matches &
#     calculate distance between matched sublandscapes
#
#################################################################################


#---------- libraries
library(terra)
library(dplyr)
library(data.table)
library(MatchIt)
library(ggplot2)
library(tidyr)
library(spatialEco)

#---------- working directory

setwd("~/data/") 



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# (1) Creating search data frame for matching
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# load preprocessed matching dataframes

set.aside.df <- readRDS("data/processed/matching/setaside.samples.df.rds")
control.df <- readRDS("data/processed/matching/control.df.rds")


# --- create searching data frame

# delete distance column in set aside sample data frame
set.aside.df <- set.aside.df[, -which(names(set.aside.df) %in% c("dist", "is.sample", "reserve", "forest.ha"))]

set.aside.df <- set.aside.df %>% 
  .[complete.cases(.[c("x", "y", "fcover", "elev.mean", "slope.mean", "aspect", "sitecond", 
                       "ftype1", "ftype2", "ftype3", "ftype4", "ftype5", "ftype6", "ftype7", "ftype8")]),] %>%
  dplyr::mutate(Research = 1) %>%
  as.data.frame(.)
head(set.aside.df)


control.df <- as.data.frame(control.df)

control.df <- control.df %>% 
  .[complete.cases(.[c("x", "y", "fcover", "elev.mean", "slope.mean", "aspect", "sitecond", 
                       "ftype1", "ftype2", "ftype3", "ftype4", "ftype5", "ftype6", "ftype7", "ftype8")]),] %>%
  dplyr::filter(fcover >= 220) %>%
  dplyr::mutate(Research = 0) %>%
  dplyr::mutate(window = 15) %>%
  as.data.frame(.)
head(control.df)

# combine set aside = Set-Aside site dataframe and managed site dataframe
search.df <- rbind(set.aside.df, control.df) 


saveRDS(search.df, "data/processed/matching/search.df.rds")
rm(list = ls())

search.df <- readRDS("data/processed/matching/search.df.rds")



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# (2) matching landscapes
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

threshold = 0.15 # % variation allowed for each forest type


# --- 1:1 

matchout <- MatchIt::matchit(Research ~ elev.mean + aspect + slope.mean, 
                             data = search.df, 
                             method = "nearest", 
                             distance = "mahalanobis",
                             exact = "sitecond",
                             caliper = c("ftype1" = threshold,
                                         "ftype2" = threshold,
                                         "ftype3" = threshold,
                                         "ftype4" = threshold,
                                         "ftype5" = threshold,
                                         "ftype6"= threshold,
                                         "ftype7"= threshold,
                                         "ftype8"= threshold),
                             std.caliper = FALSE,
                             replace = FALSE,
                             ratio=1, # looking for the one most similar managed landscape 
                             verbose = TRUE)


match_data <- MatchIt::match.data(matchout)
head(match_data)

saveRDS(match_data, "data/processed/matching/match.df.1o1_15perc.rds") # store for disturbance extraction


# --- 1:Many

matchout <- MatchIt::matchit(Research ~ elev.mean + aspect + slope.mean, 
                             data = search.df, 
                             method = "nearest", 
                             distance = "mahalanobis",
                             exact = "sitecond",
                             caliper = c("ftype1" = threshold,
                                         "ftype2" = threshold,
                                         "ftype3" = threshold,
                                         "ftype4" = threshold,
                                         "ftype5" = threshold,
                                         "ftype6"= threshold,
                                         "ftype7"= threshold,
                                         "ftype8"= threshold),
                             std.caliper = FALSE,
                             replace = FALSE,
                             ratio=30, # for the 1:Many match, 30 most similar matches are found
                             verbose = TRUE)


match_data <- MatchIt::match.data(matchout)
head(match_data)

saveRDS(match_data, "data/processed/matching/match.df.1oMany_15perc_n30.rds") # needs further filtering below!



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# (3) Distance filtering between 1:Many matches, to reduce spatial autocorrelation
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Create a terra SpatVector with the coordinates
match.points <- vect(as.data.frame(match_data), geom=c("x", "y"))

# subset to matched sites (Reserach=0) and get 10 smaples with dist <150

match.candidats <- subset(match.points, match.points$Research == 0) # managed (0)
matches.final <- subset(match.points, match.points$Research == 1) # set aside (1) 

n=10 # n samples

# i=99
# o <- as.data.frame(match.points) %>% group_by(subclass) %>% count()
# trial_sequence <- c(130:270)

for (i in unique(match.points$subclass) ) { 
  
  # condition on forest area per reserve
  
  print(i)
  
  # subset to reserves
  match.candid <- match.candidats[match.candidats$subclass == i,]
  
  if( length(match.candid$fcover)  <= n) {
   
     n <- length(match.candid$fcover) -1
     match.candid <- sf::st_as_sf(match.candid) #convert to sf object as input to subsample.distance function
     matches <- subsample.distance(match.candid, n, d=225, d.max=50000, replacement = FALSE)
     matches <- vect(matches)
     
     
     matches.final <- rbind(matches.final, matches)
     n <- 10 # to ensure we keep the goal match number 10
    
  }

    else{
    
  match.candid <- sf::st_as_sf(match.candid) #convert to sf object as input to subsample.distance function
  matches <- subsample.distance(match.candid, n, d=225, d.max=50000, replacement = FALSE)
  matches <- vect(matches)
  
  
  matches.final <- rbind(matches.final, matches)
    }
} 


crs(matches.final) <- "EPSG:25832"
matches.final.df <- as.data.frame(matches.final, geom="xy")


# calculate how often we found how many matches sublandscapes with distance filtering

matches.subclass <- matches.final.df %>% group_by(subclass) %>% count() 
matches.subclass$n <- matches.subclass$n -1
matches.subclass <- matches.subclass %>% group_by(n) %>% count()
names(matches.subclass) <- c("number.matches", "number of sublandscapes")

# 1-Many, 15% variation, inital n = 30 (with forest type matching)

# number.matches `number of sublandscapes`

            #  1                        14
            #  2                        54
            #  3                        54
            #  4                        50
            #  5                        61
            #  6                        26
            #  7                        26
            #  8                        14
            #  9                         5
            # 10                        10


saveRDS(matches.final.df, "data/processed/matching/match.df.1oMany_15perc_n30_filtered.rds") # store in match setup folder for disturbance extraction



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# (4) check on similarity in environmental features between matches &
#     calculate distance between matched sublandscapes
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# ---- check how environmental features align within matched sublandscapes ----

# load matching dataframes 

# with forest type
match_data <- readRDS("data/processed/matching/match.df.1o1_15perc.rds")
match_data <- readRDS("data/processed/matching/match.df.1oMany_15perc_n30_filtered.rds")

   
# --- compare similarity between set aside and control sublandscapes ---


# Remove columns 'x', 'y', 'subclass', and 'weights'
data_subset <- match_data %>%
  select(-x, -y, -subclass, -weights)

# Gather the data to a long format for plotting
gathered_data <- gather(data_subset, key = variable, value = value, -Research)
gathered_data <- subset(gathered_data, variable %in% c("ftype1", "ftype2", "ftype3", "ftype4", "ftype5", "ftype6", "ftype7", "ftype8", "elev.mean", "slope.mean", "aspect"))


# recode Reserach
gathered_data$Research <- recode(as.factor(gathered_data$Research),
                             "1" = "Set-Aside",
                             "0" = "Managed")
# Reverse factor levels
gathered_data$Research <- factor(gathered_data$Research, levels = rev(levels(gathered_data$Research)))


# Recode & Reorder the variables
gathered_data$variable <- factor(gathered_data$variable,
                                 levels = c("aspect", "elev.mean", "slope.mean", 
                                            "ftype1", "ftype2", "ftype3", 
                                            "ftype4", "ftype5", "ftype6", 
                                            "ftype7", "ftype8"),
                                 labels = c("Aspect", "Mean Elevation", "Mean Slope",
                                            "Spruce share", "Beech share", "Pine share",
                                            "Oak share", "Mixed-Conifer share", 
                                            "Mixed-Broadleaved share", "Mixed-Equal share", 
                                            "Other share"))


My_Theme = theme( # plotting theme
  title = element_text(size = 18),
  axis.title.x = element_text(size = 30),
  axis.text.x = element_text(size = 20,hjust=1),
  axis.text.y = element_text(size = 20),
  axis.title.y = element_text(size = 30),
  legend.key.height = unit(1, 'cm'),
  legend.title = element_text(size=30),
  legend.text = element_text(size=30),
  strip.text.x = element_text(size = 20),
  panel.spacing = unit(2, "lines"),
  legend.position = c(0.88, 0.15))
  #legend.position="top")


tiff("data/results/plots/distribution_matching_1o1.tiff", units="in", width=18, height=10, res=300)
# Plot density curves for each variable, faceted by Research (0,1)
ggplot(gathered_data, aes(x = value, fill = factor(Research))) +
  geom_density(alpha = 0.6) +
  facet_wrap(~ variable, scales = "free") +
  theme_minimal()+
  scale_fill_discrete(labels = c("Managed", "Set-Aside")) + 
  scale_fill_manual(values=met.brewer("Kandinsky"))+
  labs(fill = "", y="Density", x="Value") +
  My_Theme
dev.off()

tiff("data/results/plots/boxplot_matching_1o1.tiff", units="in", width=18, height=10, res=300)
ggplot(gathered_data, aes(x = value, y = Research, fill = factor(Research))) +
  geom_violin() +
  geom_boxplot(width = 0.2, color = "black", outlier.shape = NA) +  # Add boxplots with white outline
  facet_wrap(~ variable, scales = "free") +
  theme_minimal() +
  scale_fill_discrete(labels = c("Managed", "Set-Aside")) +
  scale_fill_manual(values=met.brewer("Kandinsky"))+
  labs(fill = "", y = "", x = "Value") +
  My_Theme +
  theme(axis.text.y = element_blank())
dev.off()



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# ----- calculate distance of matches to respective set aside sublandscapes ------


# load match data

match_data_1o1f <- readRDS("data/processed/matching/match.df.1o1_15perc.rds")
match_data_1o1f$match <- as.factor("1:1")
match_data_1oMf <- readRDS("data/processed/matching/match.df.1oMany_15perc_n30_filtered.rds")
match_data_1oMf$match <- as.factor("1:Many")


# Assuming match_data_1o1f has the correct column order
correct_order <- colnames(match_data_1o1f)

# Reorder columns of match_data_1oM
match_data_1oM <- match_data_1oM[correct_order]

# Reorder columns of match_data_1oMf
match_data_1oMf <- match_data_1oMf[correct_order]

#match_data <- rbind(as.data.frame(match_data_1o1f), as.data.frame(match_data_1o1), as.data.frame(match_data_1oM), as.data.frame(match_data_1oMf))
match_data <- rbind(as.data.frame(match_data_1o1f), as.data.frame(match_data_1oMf))
match_data$match <- as.character(match_data$match)


# Create SpatVector with the coordinates
match.points <- vect(as.data.frame(match_data), geom=c("x", "y"))

matches.final <- match.points
crs(matches.final) <- "EPSG:25832"

# Create an empty vector to store distances
distances_vector <- numeric()


# Loop over unique matches

for (match_val in unique(matches.final$match)) {
  
    # Subset points for the current 'match'
    match_subset <- matches.final[matches.final$match == match_val, ]
    
    # Loop over unique subclass values within each 'match' category
    for (subclass_val in unique(match_subset$subclass)) {
      
        # Subset points for the current subclass
        subset_points <- match_subset[match_subset$subclass == subclass_val, ]
        
        # Extract coordinates set aside landscape (Reserach = 1)
        research_coords <- subset_points[subset_points$Research == 1, c("x", "y")]
        
        # Extract coordinates of all managed landscapes (Research = 0)
        other_points_coords <- subset_points[subset_points$Research == 0, c("x", "y")]
        
        # Calculate the distance between sublandscapes
        distances <- terra::distance(research_coords, other_points_coords)
        
        distances_vector <- c(distances_vector, distances)
    }
}


# Add the distances column to the original df
matches.final$distance_to_research <- ifelse(matches.final$Research == 1, 0, distances_vector)
matches.final.df <- as.data.frame(matches.final, geom="xy")
matches.final.df$distance_to_research <- matches.final.df$distance_to_research/1000 # in km
distance_df <- subset(matches.final.df, Research == 0)


My_Theme = theme( # plotting theme
  axis.title.x = element_text(size = 25),
  axis.text.x = element_text(size = 25),
  axis.text.y = element_text(size = 25),
  axis.title.y = element_text(size = 25))


tiff("data/results/plots/matching_distance.tiff", units="in", width=18, height=6, res=300)  

ggplot(distance_df, aes(x=match,y=distance_to_research)) + 
  theme_minimal() +
  geom_boxplot() +
  stat_summary(aes(group = match),fun.y = mean, geom = "point", shape = 17, size = 8, color = "red", show.legend = FALSE)+
  stat_summary(aes(group = match, label = round(..y.., 2)), fun.y = mean, geom = "text", hjust = -0.4, show.legend = FALSE, size=7) +
  #scale_x_discrete(labels = NULL, breaks = NULL)+
  labs(x="Match type", y= "Distance between matches in km")+
  My_Theme

dev.off()

distance_df %>% group_by(match) %>% # get descriptive metrics on distances between matched landscapes
    summarize(avg = mean(distance_to_research),
              median = median(distance_to_research),
              q5 = quantile(distance_to_research, 0.05),
              q95 = quantile(distance_to_research, 0.95))


