###############################################################################
#
# Analysis script to check sensitivity of high severity disturbance threshold
# on high severity disturbance rate per mngt stratum
#
###############################################################################

#---------- libraries

library(dplyr)
library(tidyr)
library(ggplot2)
library(glmmTMB)
library(MetBrewer)
library(cowplot)
library(patchwork)
library(terra)
library(DHARMa)
library(stringr)
library(broom)

#---------- working directory

setwd("/data/") 



# --- load matching data frames and wrangle data for modelling input 


years <- 1986:2020

# preclassification of years to assign pulse and background years
dist_year_class <- read.csv("data/processed/disturbances/germany/pulse_background_classification_natural.dist.csv")

matches <- c("1o1_forest")


for (i in matches) {

  # load matching df
  #dist.df <-readRDS(paste0("data/processed/dist.extraction/",i,"/match.df_climate_dist_0.7_severity.rds"))
  dist.df <-readRDS(paste0("data/processed/dist.extraction/",i,"/match.df_climate_dist_0.9_severity.rds"))
  
  # process the data frame:

  filtered_data <- dist.df %>%

    # create column with number of not disturbed cells
    mutate(across(.cols = paste0("disturbed.cells.", years),
                  ~ n.fcover - .,  # substract disturbed cells from all forested cells in landscape to get not disturbed cells
                  .names = "no.{.col}")) %>%

    # calculate high severity rate
    rowwise() %>%
    mutate(
      total_disturbed_cells = sum(c_across(starts_with("disturbed.cells"))),
      total_high_severity = sum(c_across(starts_with("severe.disturbed.cells"))),
      #dist.rate = ((total_disturbed_cells / n.fcover) / 35), # annual disturbance rate
      high.severity.rate = total_high_severity / total_disturbed_cells,
    ) %>%
    ungroup()  # Important to remove the row wise grouping after operations


   # Filter the dataframe to include only the desired columns, depending on matching scheme
      filtered_data <- filtered_data %>% select(paste0("disturbed.cells.", years), paste0("no.disturbed.cells.", years),
           high.severity.rate, density, frequency, max.patch, mean.patch,
           Research, subclass, n.fcover, sitecond, ftype1:ftype8)


  # First, gather the disturbed.cells columns
  gathered_disturbed <- filtered_data %>%
    pivot_longer(
      cols = starts_with("disturbed.cells."),
      names_to = "year",
      names_prefix = "disturbed.cells.",
      values_to = "dist.cells"
    ) %>%
    select(-starts_with("no.disturbed.cells."))

  # Then, gather the no.disturbed.cells columns
  gathered_no_disturbed <- filtered_data %>%
    pivot_longer(
      cols = starts_with("no.disturbed.cells."),
      names_to = "year",
      names_prefix = "no.disturbed.cells.",
      values_to = "no.dist"
    ) %>%
    select(-starts_with("disturbed.cells."))

  # combine both
  gathered_data_long <- gathered_disturbed
  gathered_data_long$no.dist <-  gathered_no_disturbed$no.dist


  # add pulse year information
  gathered_data_long <- merge(gathered_data_long, dist_year_class[ , c("class", "value")],
                             by.x= "year", by.y="value", all.x=T)
  gathered_data_long$class <- as.factor(gathered_data_long$class)


  # write modelling df to file
  #saveRDS(gathered_data_long, paste0("data/processed/dist.extraction/",i,"/modelling.df_long_0.7_severity.rds"))
   saveRDS(gathered_data_long, paste0("data/processed/dist.extraction/",i,"/modelling.df_long_0.9_severity.rds"))
  
}


# --- load modelling dfs --- 

dist.df_1o1f_long_0.7 <- readRDS("data/processed/dist.extraction/1o1_forest/modelling.df_long_0.7_severity.rds")
dist.df_1o1f_long_0.9 <- readRDS("data/processed/dist.extraction/1o1_forest/modelling.df_long_0.9_severity.rds")
dist.df_1o1f_long_original <- readRDS("data/processed/dist.extraction/1o1_forest/modelling.df_long.rds")



# ------- model high severity rate metrics ---------- 


# --- 0.7 high severity threshold ---

# Summarize 'dist.cells' to create 'dist.sum' column
dist_cells_summary <- dist.df_1o1f_long_0.7 %>% 
  group_by(subclass, Research) %>% 
  summarise(dist.sum = sum(dist.cells, na.rm = TRUE)) %>% 
  ungroup()

# subset df to one row per landscape (excluding year specific disturbance information)
dist.df_1o1f_long.sub <- dist.df_1o1f_long_0.7 %>%
  left_join(dist_cells_summary, by = c("subclass", "Research")) %>%
  mutate(`high.severity.rate` = ifelse(is.na(`high.severity.rate`), 0, `high.severity.rate`)) %>%# Replace NaN with 0 in all columns before subsetting
  select(Research, subclass, frequency, density, `max.patch`, `mean.patch`, `high.severity.rate`,`dist.sum`, `sitecond`) %>%
  group_by(subclass, Research) %>%
  slice(1) %>%
  ungroup()


# subset df to only disturbed sites: 
dist.df_1o1f_long.sub.dist <- subset(dist.df_1o1f_long.sub, dist.sum > 0)



# calculate number severely disturbed and non-severly disturbed cells
dist.df_1o1f_long.sub.dist_0.7 <- dist.df_1o1f_long.sub.dist %>% 
                                      mutate(dist.severe = round(dist.sum * high.severity.rate,0),
                                      dist.no.severe = dist.sum - dist.severe)


# --- 0.8 high severity threshold --- USED in study!

# Summarize 'dist.cells' to create 'dist.sum' column
dist_cells_summary <- dist.df_1o1f_long_original %>% 
  group_by(subclass, Research) %>% 
  summarise(dist.sum = sum(dist.cells, na.rm = TRUE)) %>% 
  ungroup()

# subset df to one row per landscape (excluding year specific disturbance information)
dist.df_1o1f_long.sub <- dist.df_1o1f_long_original %>%
  left_join(dist_cells_summary, by = c("subclass", "Research")) %>%
  mutate(`high.severity.rate` = ifelse(is.na(`high.severity.rate`), 0, `high.severity.rate`)) %>%# Replace NaN with 0 in all columns before subsetting
  select(Research, subclass, frequency, density, `max.patch`, `mean.patch`, `high.severity.rate`,`dist.sum`, `sitecond`) %>%
  group_by(subclass, Research) %>%
  slice(1) %>%
  ungroup()


# subset df to only disturbed sites: 
dist.df_1o1f_long.sub.dist <- subset(dist.df_1o1f_long.sub, dist.sum > 0)


# calculate number severely disturbed and non-severly disturbed cells
dist.df_1o1f_long.sub.dist_0.8 <- dist.df_1o1f_long.sub.dist %>% 
  mutate(dist.severe = round(dist.sum * high.severity.rate,0),
         dist.no.severe = dist.sum - dist.severe)




# --- 0.9 high severity threshold ---

# Summarize 'dist.cells' to create 'dist.sum' column
dist_cells_summary <- dist.df_1o1f_long_0.9 %>% 
  group_by(subclass, Research) %>% 
  summarise(dist.sum = sum(dist.cells, na.rm = TRUE)) %>% 
  ungroup()

# subset df to one row per landscape (excluding year specific disturbance information)
dist.df_1o1f_long.sub <- dist.df_1o1f_long_0.9 %>%
  left_join(dist_cells_summary, by = c("subclass", "Research")) %>%
  mutate(`high.severity.rate` = ifelse(is.na(`high.severity.rate`), 0, `high.severity.rate`)) %>%# Replace NaN with 0 in all columns before subsetting
  select(Research, subclass, frequency, density, `max.patch`, `mean.patch`, `high.severity.rate`,`dist.sum`, `sitecond`) %>%
  group_by(subclass, Research) %>%
  slice(1) %>%
  ungroup()


# subset df to only disturbed sites: 
dist.df_1o1f_long.sub.dist <- subset(dist.df_1o1f_long.sub, dist.sum > 0)


# calculate number severely disturbed and non-severly disturbed cells
dist.df_1o1f_long.sub.dist_0.9 <- dist.df_1o1f_long.sub.dist %>% 
  mutate(dist.severe = round(dist.sum * high.severity.rate,0),
         dist.no.severe = dist.sum - dist.severe)




# --- high severity rate ---


# lower bound 
fit_severity_0.7 <- glmmTMB(cbind(dist.severe, dist.no.severe) ~ 1+ Research + (1|subclass), #+ (1|sitecond) makes effect more pronounced
                          data = dist.df_1o1f_long.sub.dist_0.7, 
                          family = betabinomial(link = "logit"))

summary(fit_severity_0.7)

simulationOutput <- simulateResiduals(fittedModel = fit_severity_0.7, plot = F)
plot(simulationOutput)

# 0.8 - used in paper
fit_severity_0.8 <- glmmTMB(cbind(dist.severe, dist.no.severe) ~ 1+ Research + (1|subclass), #+ (1|sitecond) makes effect more pronounced
                            data = dist.df_1o1f_long.sub.dist_0.8, 
                            family = betabinomial(link = "logit"))

summary(fit_severity_0.8)

simulationOutput <- simulateResiduals(fittedModel = fit_severity_0.8, plot = F)
plot(simulationOutput)


# upper bound
fit_severity_0.9 <- glmmTMB(cbind(dist.severe, dist.no.severe) ~ 1+ Research + (1|subclass), #+ (1|sitecond) makes effect more pronounced
                            data = dist.df_1o1f_long.sub.dist_0.9, 
                            family = betabinomial(link = "logit"))

summary(fit_severity_0.9)

simulationOutput <- simulateResiduals(fittedModel = fit_severity_0.9, plot = F)
plot(simulationOutput)


# --- table of model fits ---

# Create a list of your models
models_list <- list(
  severity_0.7 = list(model = fit_severity_0.7, type = "betabinomial"),
  severity_0.8 = list(model = fit_severity_0.8, type = "betabinomial"),
  severity_0.9 = list(model = fit_severity_0.9, type = "betabinomial")
)

# Function to extract effect size, CI lower, and CI upper
get_effect_size <- function(model_info) {
  # Extract the model estimate and confidence interval
  ci <- confint(model_info$model, method = "profile", parm = "beta_") %>% as.data.frame()
  ci <- ci[!grepl("^zi~", row.names(ci)), ]
  
  est <- fixef(model_info$model) 
  
  ci$estimate <- NA  # Initialize the estimate column with NA
  ci$estimate[1] <- as.numeric(est$cond["(Intercept)"])
  ci$estimate[2] <- as.numeric(est$cond["Research"])
  
  ci[2, ] <- ci[1, ] + ci[2, ] # get model estimate for unmanaged forests (Intercept (= managed ) + Research1 (effect of non-management) )
  
  
  # Apply the inverse logit transformation if it's a betabinomial model
  if (model_info$type == "betabinomial") {
    ci <- apply(ci, 2, boot::inv.logit)
  }
  
  # Apply exponential transformation for log link models!
  if (model_info$type %in% c("lognormal", "truncated_poisson")) {
    ci <- apply(ci, 2, exp)
  }
  
  return(ci)
}

# Apply the function to each model and create a data frame
results <- lapply(models_list, get_effect_size)


# Create an empty dataframe to store the combined results
result_df <- data.frame(Model = character(), Condition = character(), `2.5 %` = numeric(), `97.5 %` = numeric(), estimate = numeric(), stringsAsFactors = FALSE)

# Iterate through the list and bind rows
for(model_name in names(results)) {
  df <- as.data.frame(results[[model_name]])  # Making sure it's a dataframe
  df$Condition <- row.names(df)  # Adding Condition based on row names
  df$Model <- model_name  # Assigning the current model name
  
  # Adjusting Condition names
  df$Condition <- ifelse(df$Condition == "(Intercept)", "Managed", "Unmanaged")
  
  # Combining with the result dataframe
  result_df <- rbind(result_df, df[, c("Model", "Condition", "2.5 %", "97.5 %", "estimate")])
}


# Select and reorder the columns
result_df <- result_df %>%
  select(Model, Condition, estimate, `2.5 %`, `97.5 %`)

# Remove row names by converting them to NULL
rownames(result_df) <- NULL

# round model estimates
result_df[,c("estimate","2.5 %" ,"97.5 %") ] <- round(result_df[,c("estimate","2.5 %" ,"97.5 %") ],2)


write.table(result_df, "data/results/effects/sensitivity_high_severity.csv", 
            sep = ",", eol = ";\n", row.names = FALSE, col.names = TRUE, 
            quote = FALSE)



# > summary(fit_severity_0.7)
# Family: betabinomial  ( logit )
# Formula:          cbind(dist.severe, dist.no.severe) ~ 1 + Research + (1 | subclass)
# Data: dist.df_1o1f_long.sub.dist_0.7
# 
# AIC      BIC   logLik deviance df.resid 
# 1529.2   1543.7   -760.6   1521.2      270 
# 
# Random effects:
#   
#   Conditional model:
#   Groups   Name        Variance Std.Dev.
# subclass (Intercept) 0.05211  0.2283  
# Number of obs: 274, groups:  subclass, 195
# 
# Dispersion parameter for betabinomial family (): 4.26 
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  0.57240    0.09063   6.316 2.69e-10 ***
#   Research    -0.32405    0.12749  -2.542    0.011 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# > summary(fit_severity_0.8)
# Family: betabinomial  ( logit )
# Formula:          cbind(dist.severe, dist.no.severe) ~ 1 + Research + (1 | subclass)
# Data: dist.df_1o1f_long.sub.dist_0.8
# 
# AIC      BIC   logLik deviance df.resid 
# 1557.6   1572.0   -774.8   1549.6      270 
# 
# Random effects:
#   
#   Conditional model:
#   Groups   Name        Variance Std.Dev.
# subclass (Intercept) 0.05141  0.2267  
# Number of obs: 274, groups:  subclass, 195
# 
# Dispersion parameter for betabinomial family (): 3.41 
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)  0.18578    0.09267   2.005  0.04499 * 
#   Research    -0.39013    0.13601  -2.868  0.00413 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



# > summary(fit_severity_0.9)
# Family: betabinomial  ( logit )
# Formula:          cbind(dist.severe, dist.no.severe) ~ 1 + Research + (1 | subclass)
# Data: dist.df_1o1f_long.sub.dist_0.9
# 
# AIC      BIC   logLik deviance df.resid 
# 1524.6   1539.0   -758.3   1516.6      270 
# 
# Random effects:
#   
#   Conditional model:
#   Groups   Name        Variance  Std.Dev. 
# subclass (Intercept) 1.184e-07 0.0003441
# Number of obs: 274, groups:  subclass, 195
# 
# Dispersion parameter for betabinomial family ():  2.8 
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -0.32439    0.09581  -3.386  0.00071 ***
#   Research    -0.37893    0.14685  -2.580  0.00987 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1