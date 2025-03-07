##########################################################################################
#
# Analysis script to model disturbance metrics for managed and Set-aside forest landscapes
#
# (1) model disturbance metrics
#     (1.1) 1:1 match
#     (1.2) 1:Many match
# (2)  Extracting estimate and standard error for Research from the model summary
# (3)  Model effect direction of forest type shares on disturbance rates 
# (4)  Effect of pulse and background year 
# (5) Additional analysis & checks
#
##########################################################################################

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

setwd("~/data/") 


###############################################################################
# --- load matching dfs and wrangle data for modelling input 
###############################################################################

# run this code ones to prepare modelling data frames:

years <- 1986:2020

# preclassification of years to assign pulse and background years
dist_year_class <- read.csv("data/processed/disturbances/germany/pulse_background_classification_natural.dist.csv")

matches <- c("1o1_forest", "1oMany_forest")


for (i in matches) {

  # load matching df
  dist.df <-readRDS(paste0("data/processed/dist.extraction/",i,"/match.df_climate_dist.rds"))

  # process the dataframe:

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
  saveRDS(gathered_data_long, paste0("data/processed/dist.extraction/",i,"/modelling.df_long.rds"))

}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# summarize hectar of forest per treatment

dist.df <-readRDS(paste0("data/processed/dist.extraction/1o1_forest/match.df_climate_dist.rds"))
dist.df %>% group_by(Research) %>% summarize(sum_forest = sum(n.fcover*0.09))

# Research sum_forest
# <dbl>      <dbl>
#    0      6329.
#    1      6281.

dist.df %>% group_by(Research) %>% summarize(low_forest = min(n.fcover*0.09),
                                             high_forest = max(n.fcover*0.09))
# Research low_forest high_forest
# <dbl>      <dbl>       <dbl>
#    0       19.9        24.9
#    1       19.2        24.6
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# --- load modelling dfs --- 

dist.df_1o1f_long <- readRDS("data/processed/dist.extraction/1o1_forest/modelling.df_long.rds")
dist.df_1oManyf_long <- readRDS("data/processed/dist.extraction/1oMany_forest/modelling.df_long.rds")


# variability in disturbed area among landscapes

dist.df_1o1f_long %>% group_by(Research, subclass) %>% 
  summarise(area = sum(dist.cells*0.09)) %>% ungroup() %>%
  group_by(Research) %>% summarise(mean = mean(area),
                                   q2.5 = quantile(area, 0.025),
                                   q97.5 = quantile(area, 0.975),
                                   max = max(area))



################################################################################
# --- (1) ------------------------- model disturbance metrics -----------------
################################################################################



################################################################################
# --- (1.1) 1:1 Match
################################################################################


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% disturbance rate

dist.df_1o1f_long$Research <- as.factor(dist.df_1o1f_long$Research)

# disturbance rate

betabinom_1o1f <- glmmTMB(cbind(dist.cells, no.dist) ~ 1+ Research + (1|year) + (1|subclass), 
                  #dispformula = ~ Research, -> do I really want that?
                  data = dist.df_1o1f_long, 
                  family = betabinomial(link = "logit"))


summary(betabinom_1o1f)

model_summary <- summary(betabinom_1o1f)

fixed_effect <- fixef(betabinom_1o1f)$cond["Research1"]
se <-  model_summary$coefficients$cond["Research1", "Std. Error"]

1-exp(fixed_effect) 

# 95% confidence interval
1-exp(fixed_effect -  1.96 * se) 
1-exp(fixed_effect +  1.96 * se)

# annual dist rate unmanaged
boot::inv.logit(-7.11749 -0.25005)*100 # 0.06310211

# annual dst rate manages
boot::inv.logit(-7.11749)*100 # 0.08101424


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% density frequency, max patch size & severity 

# --- 1o1  ---

# Summarize 'dist.cells' to create 'dist.sum' column
dist_cells_summary <- dist.df_1o1f_long %>% 
  group_by(subclass, Research) %>% 
  summarise(dist.sum = sum(dist.cells, na.rm = TRUE)) %>% 
  ungroup()

# subset df to one row per landscape (excluding year specific disturbance information, as these metrics are for the whole observation period)
dist.df_1o1f_long.sub <- dist.df_1o1f_long %>%
  left_join(dist_cells_summary, by = c("subclass", "Research")) %>%
  mutate(`high.severity.rate` = ifelse(is.na(`high.severity.rate`), 0, `high.severity.rate`)) %>%# Replace NaN with 0 in all columns before subsetting
  select(Research, subclass, frequency, density, `max.patch`, `mean.patch`, `high.severity.rate`,`dist.sum`, `sitecond`) %>%
  group_by(subclass, Research) %>%
  slice(1) %>%
  ungroup()

#dist.df_1o1f_long.sub$Research <- as.factor(dist.df_1o1f_long.sub$Research)


# subset df to only disturbed sites: 
dist.df_1o1f_long.sub.dist <- subset(dist.df_1o1f_long.sub, dist.sum > 0)


# calculate number severely disturbed and non-severely disturbed cells
dist.df_1o1f_long.sub.dist <- dist.df_1o1f_long.sub.dist %>% 
                                      mutate(dist.severe = round(dist.sum * high.severity.rate,0),
                                      dist.no.severe = dist.sum - dist.severe)


# --- set model variables ---

model_df <- dist.df_1o1f_long.sub

# ---


# --- frequency (all landscapes)

ggplot(model_df, aes(x=frequency)) +
  geom_density()


fit_freq <-  glmmTMB(frequency ~ 1+ Research +  (1|subclass) , 
                     ziformula = ~Research,
                    data = model_df,
                    family = truncated_poisson(link = "log")) 


summary(fit_freq)


# --- patch density (n patches) (all landscapes)


ggplot(model_df, aes(x=density)) +
  geom_density()

fit_den <-  glmmTMB(density ~ 1+ Research + (1|subclass), #+ (1|sitecond)
                    ziformula = ~Research,
                     data = model_df,
                     family = truncated_poisson(link = "log")) 


summary(fit_den)



# --- set model variables ---

model_df <- dist.df_1o1f_long.sub.dist

# ---


# ---- patch size 

ggplot(model_df) +
  geom_density(aes(x=max.patch, group=Research,col=Research))

ggplot(model_df) +
  geom_density(aes(x=log(max.patch), group=Research,col=Research))



fit_size <-  glmmTMB(max.patch ~ 1+ Research + (1|subclass) ,  
                     data = model_df,  
                     family = lognormal(link = "log"))


summary(fit_size)



# --- high severity rate 


ggplot(model_df) +
  geom_density(aes(x=high.severity.rate, group=Research,col=Research))



fit_severity <- glmmTMB(cbind(dist.severe, dist.no.severe) ~ 1+ Research + (1|subclass), #+ (1|sitecond) makes effect more pronounced
                          data = dist.df_1o1f_long.sub.dist, 
                          family = betabinomial(link = "logit"))

summary(fit_severity)



##############################################################
# find simulation based model checks at the end of the script!
##############################################################



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# --------Table of model effects  1:1 Match ---------------#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Create a list of your models
models_list_1o1f <- list(
  betabinom_1o1f = list(model = betabinom_1o1f, type = "betabinomial"),
  fit_severity = list(model = fit_severity, type = "betabinomial"),
  fit_size = list(model = fit_size, type = "lognormal"),
  fit_den = list(model = fit_den, type = "truncated_poisson"),
  fit_freq = list(model = fit_freq, type = "truncated_poisson")
)



# Function to extract effect size, CI lower, and CI upper
get_effect_size <- function(model_info) {
  # Extract the model estimate and confidence interval
  ci <- confint(model_info$model, method = "profile", parm = "beta_") %>% as.data.frame()
  ci <- ci[!grepl("^zi~", row.names(ci)), ]
  
  est <- fixef(model_info$model) 
  
  ci$estimate <- NA  # Initialize the estimate column with NA
  ci$estimate[1] <- as.numeric(est$cond["(Intercept)"])
  ci$estimate[2] <- as.numeric(est$cond["Research1"])
  
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
results <- lapply(models_list_1o1f, get_effect_size)
#results <- sapply(models_list, get_effect_size)

# Create an empty dataframe to store the combined results
result_df <- data.frame(Model = character(), Condition = character(), `2.5 %` = numeric(), `97.5 %` = numeric(), estimate = numeric(), stringsAsFactors = FALSE)

# Iterate through the list and bind rows
for(model_name in names(results)) {
  df <- as.data.frame(results[[model_name]])  # Making sure it's a dataframe
  df$Condition <- row.names(df)  # Adding Condition based on row names
  df$Model <- model_name  # Assigning the current model name
  
  # Adjusting Condition names
  df$Condition <- ifelse(df$Condition == "(Intercept)", "Managed", "Set-aside")
  
  # Combining with the result dataframe
  result_df <- rbind(result_df, df[, c("Model", "Condition", "2.5 %", "97.5 %", "estimate")])
}


# Select and reorder the columns
result_df <- result_df %>%
  select(Model, Condition, estimate, `2.5 %`, `97.5 %`)

# Remove row names by converting them to NULL
rownames(result_df) <- NULL

# multiply estimates * 100 for disturbance rate to get %
result_df[c(1:2),c("estimate","2.5 %" ,"97.5 %") ] <- result_df[c(1:2),c("estimate","2.5 %" ,"97.5 %") ] * 100
# round model estimates
result_df[,c("estimate","2.5 %" ,"97.5 %") ] <- round(result_df[,c("estimate","2.5 %" ,"97.5 %") ],2)


write.table(result_df, "data/results/effects/model_effects_1o1f.csv", 
            sep = ",", eol = ";\n", row.names = FALSE, col.names = TRUE, 
            quote = FALSE)




################################################################################
# --- (1.2) 1:Many Match
################################################################################


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% disturbance rate

dist.df_1oManyf_long$Research <- as.factor(dist.df_1oManyf_long$Research)

betabinom_1oManyf <- glmmTMB(cbind(dist.cells, no.dist) ~ 1+ Research + (1|year) + (1|subclass), 
                          #dispformula = ~ Research, -> do I really want that?
                          data = dist.df_1oManyf_long, 
                          family = betabinomial(link = "logit"))

summary(betabinom_1oManyf)
fixed_effects <- fixef(betabinom_1oManyf)
1-exp(fixed_effects$cond[2]) 


#%%%%%%%%%%%%%%%%%%%%%%%%%%  density frequency, max patch size & severity


dist.df_1oManyf<- readRDS("data/processed/dist.extraction/1oMany_forest/match.df_climate_dist.rds")

  years <- 1986:2020
  
  dist.df_1oManyf <- dist.df_1oManyf %>%

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


   # Filter the dataframe to include only the desired columns
  
  dist.df_1oManyf <- dist.df_1oManyf %>%
    select(total_disturbed_cells, high.severity.rate, 
           density, frequency, max.patch, mean.patch,
           Research, subclass)


# subset df to one row per landscape (excluding year specific disturbance information)
dist.df_1oManyf <- dist.df_1oManyf %>%
  mutate(dist.sum =total_disturbed_cells) %>%
  mutate(`high.severity.rate` = ifelse(is.na(`high.severity.rate`), 0, `high.severity.rate`)) %>%# Replace NaN with 0 in all columns before subsetting
  select(Research, subclass, frequency, density, `max.patch`, `mean.patch`, `high.severity.rate`,`dist.sum`)

dist.df_1oManyf$Research <- as.factor(dist.df_1oManyf$Research)


# subset df to only disturbed sites: 
dist.df_1oManyf.dist <- subset(dist.df_1oManyf, dist.sum > 0)



# calculate number severely disturbed and non-severely disturbed cells
dist.df_1oManyf.dist <- dist.df_1oManyf.dist %>% 
  mutate(dist.severe = round(dist.sum * high.severity.rate,0),
         dist.no.severe = dist.sum - dist.severe)



# --- set model variables ---

model_df <- dist.df_1oManyf

# ---


# --- frequency 

ggplot(model_df, aes(x=frequency)) +
  geom_density()


fit_freq <-  glmmTMB(frequency ~ 1+ Research +  (1|subclass), # more realistic value range in max frequency than normal poission
                     ziformula = ~Research,
                     data = model_df,
                     family = truncated_poisson(link = "log")) 


summary(fit_freq)


# --- patch density (n patches) 


ggplot(model_df, aes(x=density)) +
  geom_density()

fit_den <-  glmmTMB(density ~ 1+ Research + (1|subclass),
                    ziformula = ~Research,
                    data = model_df,
                    family = truncated_poisson(link = "log")) 


summary(fit_den)



# --- set model variables ---

model_df <- dist.df_1oManyf.dist

# ---


# ---- patch size 


fit_size <-  glmmTMB(max.patch ~ 1+ Research + (1|subclass),  # seems to be the perfect fit, but doesn't simulate any variability???
                     data = model_df,  
                     family = lognormal(link = "log"))


summary(fit_size)



# --- high severity rate (only disturbed landscapes)


fit_severity <- glmmTMB(cbind(dist.severe, dist.no.severe) ~ 1+ Research + (1|subclass), 
                        data = model_df, 
                        family = betabinomial(link = "logit"))

summary(fit_severity)




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ------ Table of model effects  1:Many Match ------------#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Create a list of your models
models_list_1oManyf <- list(
  betabinom_1oManyf = list(model = betabinom_1oManyf, type = "betabinomial"),
  fit_severity = list(model = fit_severity, type = "betabinomial"),
  fit_size = list(model = fit_size, type = "lognormal"),
  fit_den = list(model = fit_den, type = "truncated_poisson"),
  fit_freq = list(model = fit_freq, type = "truncated_poisson")
)


# Function to extract effect size, CI lower, and CI upper
get_effect_size <- function(model_info) {
  # Extract the model estimate and confidence interval
  ci <- confint(model_info$model, method = "profile", parm = "beta_") %>% as.data.frame()
  ci <- ci[!grepl("^zi~", row.names(ci)), ]
  
  est <- fixef(model_info$model) 
  
  ci$estimate <- NA  # Initialize the estimate column with NA
  ci$estimate[1] <- as.numeric(est$cond["(Intercept)"])
  ci$estimate[2] <- as.numeric(est$cond["Research1"])
  
  ci[2, ] <- ci[1, ] + ci[2, ] # get model estimate for Set-aside forests (Intercept (= managed ) + Research1 (effect of non-management) )
  
  
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
results <- lapply(models_list_1oManyf, get_effect_size)
#results <- sapply(models_list, get_effect_size)

# Create an empty dataframe to store the combined results
result_df <- data.frame(Model = character(), Condition = character(), `2.5 %` = numeric(), `97.5 %` = numeric(), estimate = numeric(), stringsAsFactors = FALSE)

# Iterate through the list and bind rows
for(model_name in names(results)) {
  df <- as.data.frame(results[[model_name]])  # Making sure it's a dataframe
  df$Condition <- row.names(df)  # Adding Condition based on row names
  df$Model <- model_name  # Assigning the current model name
  
  # Adjusting Condition names
  df$Condition <- ifelse(df$Condition == "(Intercept)", "Managed", "Set-aside")
  
  # Combining with the result dataframe
  result_df <- rbind(result_df, df[, c("Model", "Condition", "2.5 %", "97.5 %", "estimate")])
}


# Select and reorder the columns
result_df <- result_df %>%
  select(Model, Condition, estimate, `2.5 %`, `97.5 %`)

# Remove row names by converting them to NULL
rownames(result_df) <- NULL

# multiply estimates * 100 for disturbance rate to get %
result_df[c(1:2),c("estimate","2.5 %" ,"97.5 %") ] <- result_df[c(1:2),c("estimate","2.5 %" ,"97.5 %") ] * 100
# round model estimates
result_df[,c("estimate","2.5 %" ,"97.5 %") ] <- round(result_df[,c("estimate","2.5 %" ,"97.5 %") ],2)


write.table(result_df, "data/results/effects/model_effects_1oManyf.csv", 
            sep = ",", eol = ";\n", row.names = FALSE, col.names = TRUE, 
            quote = FALSE)


##################################################################################
# (2) ---- Extracting estimate and standard error for Research from the model summary
##################################################################################

# 1:1 Match

# --- disturbance rate ---

model <- betabinom_1o1f

estimate_research <- summary(model)$coefficients$cond["Research1", "Estimate"]
std_error_research <- summary(model)$coefficients$cond["Research1", "Std. Error"]
estimate_intercept <- summary(model)$coefficients$cond["(Intercept)", "Estimate"]
std_error_intercept <- summary(model)$coefficients$cond["(Intercept)", "Std. Error"]


# Converting logit estimates to probabilities (in percent)
prob_research <- boot::inv.logit(estimate_intercept + estimate_research) * 100
prob_intercept <- boot::inv.logit(estimate_intercept) * 100

# Converting logit confidence intervals to probabilities (in percent)
lower_research <- boot::inv.logit((estimate_intercept + estimate_research) - 1.96 * std_error_research) * 100
upper_research <- boot::inv.logit((estimate_intercept + estimate_research) + 1.96 * std_error_research) * 100
lower_intercept <- boot::inv.logit(estimate_intercept - 1.96 * std_error_intercept) * 100
upper_intercept <- boot::inv.logit(estimate_intercept + 1.96 * std_error_intercept) * 100

# Create a data frame for plotting
effect_size_df <- data.frame(
  Research = c("Managed", "Set-aside"),
  estimate = c(prob_intercept, prob_research),
  lower = c(lower_intercept, lower_research),
  upper = c(upper_intercept, upper_research)
)


My_Theme = theme(
  axis.title.x = element_text(size = 50),
  axis.text.x = element_text(size = 50),
  axis.text.y = element_text(size = 50),
  axis.title.y = element_text(size = 60),
  strip.text.x = element_text(size = 50))


dist.df_1o1f_long <- dist.df_1o1f_long %>%
  mutate(Research = recode(Research, '0' = "Managed", '1' = "Set-aside"))
dist.df_1o1f_long$Research <- as.factor(dist.df_1o1f_long$Research)
dist.df_1o1f_long$Research <- factor(dist.df_1o1f_long$Research, levels = rev(levels(dist.df_1o1f_long$Research)))
dist.df_1o1f_long$dist.rate <- dist.df_1o1f_long$dist.cells / (dist.df_1o1f_long$dist.cells + dist.df_1o1f_long$no.dist)

#exclude 0 disturbance rate, as this is captured by the zero-inflation model
dist.df_1o1f_long_nona <- dist.df_1o1f_long %>% 
  filter(dist.rate != 0) 

My_Theme = theme( 
  axis.title.x = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_text(size = 25),
  axis.title.y = element_text(size = 30),
  strip.text.x = element_text(size = 25))

# Plot
pdist <- ggplot(dist.df_1o1f_long_nona, aes(x=Research, y=dist.rate, fill=Research)) +
  geom_violin() + 
  scale_fill_manual(values=met.brewer("Kandinsky"), guide=FALSE)+
  geom_point(data= effect_size_df, aes(x = Research, y = estimate),size=4) +
  geom_errorbar(data= effect_size_df, aes(x = Research, y = estimate,ymin = lower, ymax = upper), width = 0.1, size=1.5) +
  labs(x = "", y = expression("Disturbance rate [" * "%" * " " * year^-1 * "]"))+
  theme_bw() +
  coord_cartesian(ylim = c(0, 0.3)) +
  My_Theme 
pdist

#ggsave("data/results/plots/effect_dist.rate_0excluded.tiff", pdist, width = 20, height = 15) 


# --- Frequency ---

model <- fit_freq

estimate_research <- summary(model)$coefficients$cond["Research1", "Estimate"]
std_error_research <- summary(model)$coefficients$cond["Research1", "Std. Error"]
estimate_intercept <- summary(model)$coefficients$cond["(Intercept)", "Estimate"]
std_error_intercept <- summary(model)$coefficients$cond["(Intercept)", "Std. Error"]


# Converting logit estimates to probabilities (in percent)
prob_research <- exp(estimate_intercept + estimate_research)
prob_intercept <- exp(estimate_intercept)

# Converting logit confidence intervals to probabilities (in percent)
lower_research <- exp((estimate_intercept + estimate_research) - 1.96 * std_error_research)
upper_research <- exp((estimate_intercept + estimate_research) + 1.96 * std_error_research)
lower_intercept <- exp(estimate_intercept - 1.96 * std_error_intercept) 
upper_intercept <- exp(estimate_intercept + 1.96 * std_error_intercept)

# Create a data frame for plotting
effect_size_df <- data.frame(
  Research = c("Managed", "Set-aside"),
  estimate = c(prob_intercept, prob_research),
  lower = c(lower_intercept, lower_research),
  upper = c(upper_intercept, upper_research)
)


My_Theme = theme(
  axis.title.x = element_text(size = 50),
  axis.text.x = element_text(size = 50),
  axis.text.y = element_text(size = 50),
  axis.title.y = element_text(size = 50),
  strip.text.x = element_text(size = 50))


dist.df_1o1f_long.sub <- dist.df_1o1f_long.sub %>%
  mutate(Research = recode(Research, '0' = "Managed", '1' = "Set-aside"))
dist.df_1o1f_long.sub$Research <- as.factor(dist.df_1o1f_long.sub$Research)
dist.df_1o1f_long.sub$Research <- factor(dist.df_1o1f_long.sub$Research, levels = rev(levels(dist.df_1o1f_long.sub$Research)))

My_Theme = theme( 
  axis.title.x = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_text(size = 25),
  axis.title.y = element_text(size = 30),
  strip.text.x = element_text(size = 25))

# Plot
pfreq <- ggplot(dist.df_1o1f_long.sub, aes(x=Research, y=frequency, fill=Research)) +
  geom_violin(adjust = 2) + 
  scale_fill_manual(values=met.brewer("Kandinsky"), guide=FALSE)+
  geom_point(data= effect_size_df, aes(x = Research, y = estimate),size=4) +
  geom_errorbar(data= effect_size_df, aes(x = Research, y = estimate,ymin = lower, ymax = upper), width = 0.1, size=1.5) +
  labs(x = "", y = expression("Frequency"))+
  theme_bw() +
  My_Theme 
pfreq


#ggsave("data/results/model_checks/effect_frequency_1o1_with_obs.tiff", pfreq, width = 20, height = 15) 


# --- Density ---

model <- fit_den

estimate_research <- summary(model)$coefficients$cond["Research1", "Estimate"]
std_error_research <- summary(model)$coefficients$cond["Research1", "Std. Error"]
estimate_intercept <- summary(model)$coefficients$cond["(Intercept)", "Estimate"]
std_error_intercept <- summary(model)$coefficients$cond["(Intercept)", "Std. Error"]


# Converting logit estimates to probabilities (in percent)
prob_research <- exp(estimate_intercept + estimate_research)
prob_intercept <- exp(estimate_intercept)

# Converting logit confidence intervals to probabilities (in percent)
lower_research <- exp((estimate_intercept + estimate_research) - 1.96 * std_error_research)
upper_research <- exp((estimate_intercept + estimate_research) + 1.96 * std_error_research)
lower_intercept <- exp(estimate_intercept - 1.96 * std_error_intercept) 
upper_intercept <- exp(estimate_intercept + 1.96 * std_error_intercept)

# Create a data frame for plotting
effect_size_df <- data.frame(
  Research = c("Managed", "Set-aside"),
  estimate = c(prob_intercept, prob_research),
  lower = c(lower_intercept, lower_research),
  upper = c(upper_intercept, upper_research)
)


My_Theme = theme( 
  axis.title.x = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_text(size = 25),
  axis.title.y = element_text(size = 30),
  strip.text.x = element_text(size = 25))

# Plot
pdens <- ggplot(dist.df_1o1f_long.sub, aes(x=Research, y=density, fill=Research)) +
  geom_violin() + 
  scale_fill_manual(values=met.brewer("Kandinsky"), guide=FALSE)+
  geom_point(data= effect_size_df, aes(x = Research, y = estimate),size=4) +
  geom_errorbar(data= effect_size_df, aes(x = Research, y = estimate,ymin = lower, ymax = upper), width = 0.1, size=1.5) +
  labs(x = "", y = expression("N disturbance patches"))+
  theme_bw() +
  scale_y_sqrt() +
  My_Theme 
pdens


#ggsave("data/results/model_checks/effect_density_1o1f_with_obs.tiff", pdens, width = 20, height = 15) 


# --- patch size ---

model <- fit_size

estimate_research <- summary(model)$coefficients$cond["Research1", "Estimate"]
std_error_research <- summary(model)$coefficients$cond["Research1", "Std. Error"]
estimate_intercept <- summary(model)$coefficients$cond["(Intercept)", "Estimate"]
std_error_intercept <- summary(model)$coefficients$cond["(Intercept)", "Std. Error"]


# Converting logit estimates to probabilities (in percent)
prob_research <- exp(estimate_intercept + estimate_research)
prob_intercept <- exp(estimate_intercept)

# Converting logit confidence intervals to probabilities (in percent)
lower_research <- exp((estimate_intercept + estimate_research) - 1.96 * std_error_research)
upper_research <- exp((estimate_intercept + estimate_research) + 1.96 * std_error_research)
lower_intercept <- exp(estimate_intercept - 1.96 * std_error_intercept) 
upper_intercept <- exp(estimate_intercept + 1.96 * std_error_intercept)

# Create a data frame for plotting
effect_size_df <- data.frame(
  Research = c("Managed", "Set-aside"),
  estimate = c(prob_intercept, prob_research),
  lower = c(lower_intercept, lower_research),
  upper = c(upper_intercept, upper_research)
)


dist.df_1o1f_long.sub.dist <- dist.df_1o1f_long.sub.dist %>%
  mutate(Research = recode(Research, '0' = "Managed", '1' = "Set-aside"))
dist.df_1o1f_long.sub.dist$Research <- as.factor(dist.df_1o1f_long.sub.dist$Research)
dist.df_1o1f_long.sub.dist$Research <- factor(dist.df_1o1f_long.sub.dist$Research, levels = rev(levels(dist.df_1o1f_long.sub.dist$Research)))


My_Theme = theme( # for Poster plot
  axis.title.x = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_text(size = 25),
  axis.title.y = element_text(size = 30),
  strip.text.x = element_text(size = 25))

# Plot
psize <- ggplot(dist.df_1o1f_long.sub.dist, aes(x=Research, y=mean.patch, fill=Research)) + #max.patch
  geom_violin() + 
  scale_fill_manual(values=met.brewer("Kandinsky"), guide=FALSE)+
  geom_point(data= effect_size_df, aes(x = Research, y = estimate),size=4) +
  geom_errorbar(data= effect_size_df, aes(x = Research, y = estimate,ymin = lower, ymax = upper), width = 0.1, size=1.5) +
  labs(x = "", y = expression("Max patch size (ha)"))+
  theme_bw() +
  scale_y_sqrt() +
  My_Theme 
psize



# --- High severity rate ---

model <- fit_severity

estimate_research <- summary(model)$coefficients$cond["Research1", "Estimate"]
std_error_research <- summary(model)$coefficients$cond["Research1", "Std. Error"]
estimate_intercept <- summary(model)$coefficients$cond["(Intercept)", "Estimate"]
std_error_intercept <- summary(model)$coefficients$cond["(Intercept)", "Std. Error"]


# Converting logit estimates to probabilities (in percent)
prob_research <- boot::inv.logit(estimate_intercept + estimate_research)*100
prob_intercept <- boot::inv.logit(estimate_intercept)*100

# Converting logit confidence intervals to probabilities (in percent)
lower_research <- boot::inv.logit((estimate_intercept + estimate_research) - 1.96 * std_error_research)*100
upper_research <- boot::inv.logit((estimate_intercept + estimate_research) + 1.96 * std_error_research)*100
lower_intercept <- boot::inv.logit(estimate_intercept - 1.96 * std_error_intercept) *100
upper_intercept <- boot::inv.logit(estimate_intercept + 1.96 * std_error_intercept)*100

# Create a data frame for plotting
effect_size_df <- data.frame(
  Research = c("Managed", "Set-aside"),
  estimate = c(prob_intercept, prob_research),
  lower = c(lower_intercept, lower_research),
  upper = c(upper_intercept, upper_research)
)

dist.df_1o1f_long.sub.dist <- dist.df_1o1f_long.sub.dist %>%
  mutate(Research = recode(Research, '0' = "Managed", '1' = "Set-aside"),
         high.severity.rate.perc = high.severity.rate * 100)


My_Theme = theme(
  axis.title.x = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_text(size = 25),
  axis.title.y = element_text(size = 30),
  strip.text.x = element_text(size = 25))

# Plot
psev <- ggplot(dist.df_1o1f_long.sub.dist, aes(x=Research, y=high.severity.rate.perc, fill=Research)) +
  geom_violin() + 
  scale_fill_manual(values=met.brewer("Kandinsky"), guide=FALSE)+
  geom_point(data= effect_size_df, aes(x = Research, y = estimate),size=4) +
  geom_errorbar(data= effect_size_df, aes(x = Research, y = estimate,ymin = lower, ymax = upper), width = 0.1, size=1.5) +
  labs(x = "", y = expression("High severity rate (%)"))+
  theme_bw() +
  My_Theme 
psev

#ggsave("data/results/model_checks/effect_severity_1o1f_with_obs.tiff", psev, width = 20, height = 15) 



# ---- Create a plot to extract the legend
psev_legend <- ggplot(dist.df_1o1f_long.sub.dist, aes(x = Research, y = high.severity.rate.perc, fill = Research)) +
  geom_violin() + 
  scale_fill_manual(values = met.brewer("Kandinsky")) +
  theme(legend.position = "top", legend.title = element_blank(),legend.text = element_text(size = 30), legend.key.width = unit(1.6, "cm"))

# Extract the legend
components <- get_plot_component(psev_legend, "guide-box", return_all = TRUE)
legend <- components[[4]]


#------ Add significance bars

pdist <- pdist + 
  geom_segment(aes(x = 1, xend = 2, y = 0.15, yend = 0.15), size = 1.5) + 
  geom_segment(aes(x = 1, xend = 1, y = 0.08, yend = 0.15), size = 1.5) + 
  geom_segment(aes(x = 2, xend = 2, y = 0.125, yend = 0.15), size = 1.5) + 
  annotate("text", x = 1.5, y = 0.16, label = "*", size = 12)

psev <- psev + 
  geom_segment(aes(x = 1, xend = 2, y = 70, yend = 70), size = 1.5) + 
  geom_segment(aes(x = 1, xend = 1, y = 55, yend = 70), size = 1.5) + 
  geom_segment(aes(x = 2, xend = 2, y = 65, yend = 70), size = 1.5) + 
  annotate("text", x = 1.5, y = 75, label = "*", size = 12)

# Combine plots as before
empty_plot <- ggplot() + theme_void() + theme(panel.background = element_blank(), panel.grid = element_blank())

combined_plot <- cowplot::plot_grid(
  pdist, pfreq, pdens, psize, psev, 
  ncol = 5
)

# Create a plot grid with the legend at the top
final_plot <- plot_grid(legend, combined_plot, ncol = 1, rel_heights = c(0.1, 1))


# Save the final plot
ggsave("data/results/plots/dist_metrics_observations_model.estimates_significance.tiff", final_plot, width = 20, height = 6,limitsize = FALSE,  dpi = 300) 
ggsave("data/results/plots/dist_metrics_observations_model.estimates_significance.png", final_plot, width = 20, height = 6,limitsize = FALSE,  dpi = 300) 


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# --- (3)--- model effect direction of forest type shares on disturbance rates ---
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



fit_all <- glmmTMB(cbind(dist.cells, no.dist) ~ 1+ 
                     Research*ftype1 + 
                     Research*ftype2 +
                     Research*ftype5 +
                     Research*ftype6 +
                     Research*ftype7 +
                     (1|year) + (1|subclass),
                   data = dist.df_1o1f_long, 
                   family = betabinomial(link = "logit"))
summary(fit_all)



# simulate with this model

################################################## create simulation data frame

#set.seed(12)

# Define discrete values for forest types
ftype_values <- seq(0, 1, length.out = 30)

# Number of repetitions for each combination
n_repeats <- 1000

# Create an empty data frame for predictions
pred_data <- data.frame()

# Loop over Research values (0 and 1)
for (research in c(0, 1)) {
  # Create strict tradeoff combinations (one ftype is 1, others are 0)
  tradeoff_combinations <- diag(1, nrow = 5, ncol = 5) # 5 ftypes (strict tradeoff)
  
  # Create random combinations with discrete values
  random_combinations <- expand.grid(
    ftype1 = ftype_values,
    ftype2 = ftype_values,
    ftype5 = ftype_values,
    ftype6 = ftype_values,
    ftype7 = ftype_values
  )
  
  # Filter out invalid rows where the sum of ftypes is not 1 (partial tradeoff)
  random_combinations <- random_combinations[rowSums(random_combinations) == 1, ]
  
  # Combine strict tradeoff and valid random combinations
  all_combinations <- rbind(
    data.frame(ftype1 = tradeoff_combinations[, 1],
               ftype2 = tradeoff_combinations[, 2],
               ftype5 = tradeoff_combinations[, 3],
               ftype6 = tradeoff_combinations[, 4],
               ftype7 = tradeoff_combinations[, 5]),
    random_combinations
  )
  
  # Add Research column to all combinations
  combined_data <- cbind(all_combinations, Research = research)
  
  # Bootstrapping for equal share of all ftypes in all shares
  
  bootstrapped_data_list <- list() # Initialize list to store bootstrapped data for each ftype
  
  # Bootstrapping for all ftypes
  set.seed(123)  # Set a seed for reproducibility
  for (ftype in c("ftype1", "ftype2", "ftype5", "ftype6", "ftype7")) {
    resampled_ftype <- combined_data %>%
      group_by_at(ftype) %>%
      summarise(
        bootstrapped_samples = list(sample_n(cur_data(), size = n_repeats, replace = TRUE))
      ) %>%
      unnest(bootstrapped_samples) %>%
      ungroup()
    
    # Store the resampled data for the current ftype
    bootstrapped_data_list[[ftype]] <- resampled_ftype
  }
  
  # Merge all bootstrapped dataframes
  resampled_data <- bind_rows(bootstrapped_data_list)
  
  # Append to the final prediction data frame
  pred_data <- rbind(pred_data, resampled_data)
}
##################################################


# Randomly assign 'year' and 'subclass' values for each replication
set.seed(66) # Set a seed for reproducibility
pred_data$year <- sample(1986:2020, size = nrow(pred_data), replace = TRUE)
pred_data$subclass <- sample(1:314, size = nrow(pred_data), replace = TRUE)


# Replace 'fit1.2' with your model object
pred_data$predicted <- predict(fit_all, newdata = pred_data, re.form = NULL, type = "response")


# aggregate the predictions for plotting:

agg_predictions <- pred_data %>%
  group_by(ftype1, ftype2,  ftype5, ftype6, ftype7, Research) %>%
  summarize(
    mean_predicted = mean(predicted),
    lower_quantile = quantile(predicted, probs = 0.05),
    upper_quantile = quantile(predicted, probs = 0.95),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = starts_with("ftype"),
    names_to = "group",
    values_to = "ftype"
  )%>%
  group_by(group, ftype, Research) %>%
  summarize(
    avg_mean_predicted = mean(mean_predicted),
    avg_lower_quantile = min(lower_quantile),
    avg_upper_quantile = mean(upper_quantile),
    .groups = "drop"
  )



agg_predictions$Research <- as.factor(agg_predictions$Research)
agg_predictions$group <- as.factor(agg_predictions$group)

# change labels for plotting
agg_predictions <- agg_predictions %>%
  mutate(
    group = case_when(
      group == "ftype1" ~ "Spruce",
      group == "ftype2" ~ "Beech",
      group == "ftype5" ~ "Conifer - dominated",
      group == "ftype6" ~ "Broadleaved - dominated",
      group == "ftype7" ~ "Mixed",
      TRUE ~ group # Keep other values unchanged
    ),
    Research = case_when(
      Research == 0 ~ "Managed",
      Research == 1 ~ "Set-aside",
      TRUE ~ as.character(Research) # Keep other values unchanged
    )
  )


# set plotting theme
My_Theme = theme(
  title = element_text(size = 18),
  axis.title.x = element_text(size = 30),
  axis.text.x = element_text(size = 23,hjust=1),
  axis.text.y = element_text(size = 23),
  axis.title.y = element_text(size = 30),
  legend.key.height = unit(1, 'cm'),
  legend.key.width = unit(2, 'cm'),
  legend.title = element_text(size=30),
  legend.text = element_text(size=30),
  strip.text.x = element_text(size = 28),
  strip.background = element_rect(fill = "transparent"), # Transparent background
  panel.spacing = unit(2, "lines"),
  # Legend position in the empty space (bottom left)
  legend.position = c(0.65, 0.05), 
  legend.justification = c("left", "bottom"),
  legend.box.margin = margin(t = 10, r = 10, b = 10, l = 10))


# prepare original data points for the rugs
dist.df_long_sub <- dist.df_long %>%
  mutate(
    group = case_when(
      group == "ftype1" ~ "Spruce",
      group == "ftype2" ~ "Beech",
      group == "ftype5" ~ "Conifer - dominated",
      group == "ftype6" ~ "Broadleaved - dominated",
      group == "ftype7" ~ "Mixed",
      TRUE ~ group # Keep other values unchanged
    ),
    Research = case_when(
      Research == 0 ~ "Managed",
      Research == 1 ~ "Set-aside",
      TRUE ~ as.character(Research) # Keep other values unchanged
    )
  ) %>%
  filter(!group %in% c("ftype3", "ftype4"))

# bring forest type labels into right order
agg_predictions$group <- factor(agg_predictions$group, 
                                levels = c("Spruce", 
                                           "Conifer - dominated", 
                                           "Beech", 
                                           "Broadleaved - dominated", 
                                           "Mixed"))
dist.df_long_sub$group <- factor(dist.df_long_sub$group, 
                                 levels = c("Spruce", 
                                            "Conifer - dominated", 
                                            "Beech", 
                                            "Broadleaved - dominated", 
                                            "Mixed"))

# plot
ftype_dist<-  ggplot(agg_predictions, aes(x = ftype, y = avg_mean_predicted, color = as.factor(Research))) +
  geom_ribbon(aes(ymin = avg_lower_quantile, ymax = avg_upper_quantile, fill = factor(Research)), alpha = 0.4) +
  geom_line(size = 2) +
  geom_rug(
    data = dist.df_long_sub,
    aes(x = value,  color = Research), 
    sides = "b",
    inherit.aes = FALSE
  ) +
  # Panel layout: 3 rows and 2 columns
  facet_wrap(~group, nrow = 3, ncol = 2, drop = FALSE) +
  scale_color_manual(values = c("#ce9642", "#3b7c70")) +
  scale_fill_manual(values = c("#ce9642", "#3b7c70")) +
  labs(
    x = "Forest type share [%]",
    y = expression("Predicted disturbance rate [" * "%" * " " * year^-1 * "]"),
    fill = "Treatment") +
  guides(colour = FALSE) +
  #theme_classic() +
  theme_test() +
  My_Theme 
ftype_dist


ggsave("results/plots/ftype_dist_simulations.tiff", ftype_dist, width = 15, height = 13, dpi = 200)
ggsave("results/plots/ftype_dist_simulations.png", ftype_dist, width = 15, height = 13, dpi = 150)


################################################################################
# (4)---- effect of pulse and background year -------
################################################################################

#dist.df_1o1f_long$Research <- as.factor(dist.df_1o1f_long$Research)

betabinom_1o1f.class <- glmmTMB(cbind(dist.cells, no.dist) ~ 1+ Research*class + (1|year) + (1|subclass), 
                          data = dist.df_1o1f_long, 
                          family = betabinomial(link = "logit"))

summary(betabinom_1o1f.class)


# --- bootstrap for model coefficient uncertainty estimation ---


# Run the bootstrap and modeling

data <- dist.df_1o1f_long


n_bootstrap <- 100
sample_size <- round(0.8 * nrow(data))
estimates_df <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(estimates_df) <- c("coeffs", "Estimate", "Std.Error", "z.value", "Pr(>|z|)", "sim")


for (i in 1:n_bootstrap) {
  
  # Bootstrap sample
  sample_data <- data[sample(nrow(data), sample_size, replace = TRUE), ]
  
  # Fit the model
  fit <- glmmTMB(cbind(dist.cells, no.dist) ~ 1 + Research*class + (1|year) + (1|subclass), 
                 data = sample_data, 
                 family = betabinomial(link = "logit"))
  
  # Store estimates
  estimates <- as.data.frame(summary(fit)$coeff$cond)
  estimates$sim <- as.numeric(i)
  estimates <- tibble::rownames_to_column(estimates, var = "coeffs")
  
  estimates_df <- rbind(estimates_df, estimates)
  
}

saveRDS(estimates_df, "data/processed/modelling/sim.estimates.pulse.back_df.rds")
#estimates_df <- readRDS("data/processed/modelling/sim.estimates.pulse.back_df.rds" )

# Plot the distributions of the estimates
ggplot(estimates_df, aes(x = Estimate, fill = coeffs)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ coeffs, scales = "free") +
  labs(title = "Distribution of Model Estimates", x = "Estimate", y = "Density")

# calculate model parameters
parameter_df <- estimates_df %>%
  group_by(sim) %>%
  summarize(
    managed_back = sum(Estimate[coeffs == "(Intercept)"]),
    unmanaged_back = sum(Estimate[coeffs == "(Intercept)"], Estimate[coeffs == "Research1"], na.rm = TRUE),
    managed_pulse = sum(Estimate[coeffs == "(Intercept)"], Estimate[coeffs == "classpulse"], na.rm = TRUE),
    unmanaged_pulse = sum(Estimate[coeffs == "(Intercept)"], Estimate[coeffs == "Research1"], 
                          Estimate[coeffs == "classpulse"], Estimate[coeffs == "Research1:classpulse"], na.rm = TRUE)
  )

# convert back from logit scale
parameter_df.inv <- parameter_df %>%
  mutate(across(c(managed_back, unmanaged_back, managed_pulse, unmanaged_pulse), boot::inv.logit))

summary(parameter_df.inv)

# create plotting data frame
parameter_df.inv.long <- parameter_df.inv %>%
  pivot_longer(
    cols = c(managed_back, unmanaged_back, managed_pulse, unmanaged_pulse),
    names_to = "parameter",
    values_to = "value"
  )%>%
  mutate(year.type = factor(ifelse(str_detect(parameter, "pulse"), "Pulse year", "Background year")))%>%
  mutate(treatment = factor(ifelse(str_detect(parameter, "unmanaged"), "Set-aside", "Managed"))) %>% 
  mutate(value = value*100) #%


parameter_df.inv.long %>% group_by(year.type, treatment) %>% 
  summarise(mean_rate = mean(value))

# year.type       treatment mean_rate

# Background year Managed      0.0426
# Background year Set-aside    0.0404
# Pulse year      Managed      0.228 
# Pulse year      Set-aside    0.135 

# ----- plot the effects: -----

My_Theme = theme(
  title = element_text(size = 18),
  axis.title.x = element_text(size = 30),
  axis.text.x = element_text(size = 25,hjust=0.5, margin = margin(t = 28)),
  axis.text.y = element_text(size = 20),
  axis.title.y = element_text(size = 30),
  legend.key.height = unit(1, 'cm'),
  legend.key.width = unit(2, 'cm'),
  legend.title = element_text(size=25),
  legend.text = element_text(size=25),
  legend.position = c(0.3,0.9))

parameter_df.inv.long$treatment <- factor(parameter_df.inv.long$treatment, 
                                          levels = c("Set-aside", "Managed"))

# Plot the distributions of the estimates
p1 <- ggplot(parameter_df.inv.long, aes(x = year.type, y= value, fill = treatment)) +
  geom_boxplot() +
  labs(x = "", y = expression("Disturbance rate [" * "%" * " " * year^-1 * "]"), fill="")+
  scale_fill_manual(values = c("#3b7c70","#ce9642" )) +
  theme_classic()+
  My_Theme
p1


#ggsave("data/results/plots/bootstrap_estimates_eff_pulse_back.tiff", p1, width = 20, height = 15) 


# timeline plot of disturbance rates 

dist.df_1o1f_long$dist.rate <- dist.df_1o1f_long$dist.cells / dist.df_1o1f_long$n.fcover

# Calculate mean and quantiles for each research and subclass
summary_dat <- dist.df_1o1f_long %>%
  group_by(Research, year) %>%
  summarize(mean_dist_rate = mean(dist.rate),
            q5_dist_rate = quantile(dist.rate, probs = 0.05),
            q95_dist_rate = quantile(dist.rate, probs = 0.95))


My_Theme = theme(
  axis.title.x = element_text(size = 30),
  axis.text.x = element_text(size = 22, angle = 45,hjust=1),
  axis.text.y = element_text(size = 25),
  axis.title.y = element_text(size = 30),
  legend.key.height = unit(1, 'cm'),
  legend.key.width = unit(2, 'cm'),
  legend.title = element_text(size=30),
  legend.text = element_text(size=25),
  strip.text.x = element_text(size = 20),
  panel.spacing = unit(2, "lines"),
  legend.position = c(0.75, 0.8))

pulse_years <- c(1990, 2000, 2007, 2018,2019,2020)
summary_dat$year <- as.numeric(summary_dat$year)


rect_data <- data.frame(
  xmin = pulse_years - 0.5, # Slightly before the year starts
  xmax = pulse_years + 0.5, # Slightly after the year ends
  ymin = -Inf, # Cover the whole y-axis
  ymax = Inf
)

# Plot the disturbance rate per year and Research with mean 
plot_time <- ggplot(summary_dat, aes(x = year, y = mean_dist_rate, color = Research)) +
  geom_rect(data = rect_data, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="lightgrey", alpha=0.5, inherit.aes = FALSE) + # Add this line
  geom_line(linewidth = 3) +
  geom_line(linewidth = 3) +
  labs(x = "Year", y = expression("Disturbance rate [" * "%" * " " * year^-1 * "]"), color = "") +
  scale_color_manual(values=met.brewer("Kandinsky"))+
  theme_classic()+
  theme(legend.position = "top") +
  scale_x_continuous(breaks = seq(min(summary_dat$year), max(summary_dat$year), by = 2)) +
  My_Theme
plot_time

combined_plot <- plot_grid(plot_time, p1, labels = c("(A)", "(B)"), ncol = 2)

combined_plot <- plot_grid(plot_time, p1, 
                           labels = c("(A)", "(B)"), 
                           ncol = 2, 
                           rel_widths = c(2.1, 1),  # Relative widths
                           label_size = 30,
                           label_colour = "gray45",
                           label_x = c(0, - 0.05)) #, hjust = -1.5
combined_plot

# Save the final plot
ggsave("data/results/plots/pulse_years_obs_prediction4.tiff", combined_plot, width = 20, height = 7.5) # Adjust dimensions as needed



################################################################################
# (5)---- additional analysis  -------
################################################################################

#################################################################################
# --- compare random effects ----
#################################################################################

ranef(betabinom_1o1f)

# random effects year

intercepts <- ranef(betabinom_1o1f) %>% 
  as_tibble() %>%
  subset(grpvar== "year")%>%
  mutate(group = paste("year", 1986:2020)) %>%
  rename(intercept_variation = condval)


intercepts


intercepts_year <- ggplot(intercepts, 
                          aes(x = group, 
                              y = intercept_variation)) +
  geom_point() +
  geom_segment(aes(xend = group, 
                   y = 0, 
                   yend = intercept_variation, 
                   group = group)) +
  coord_flip() +
  theme_classic() +
  labs(x = NULL,
       y = "Difference in intercept") +
  geom_hline(yintercept = 0, linetype = "dashed", col = "grey")

ggsave("data/processed/modelling/intercepts.year.tiff", intercepts_year , width = 6, height = 7)


ggplot(intercepts, 
       aes(x = grp, 
           y = intercept_variation)) +
  geom_point() +
  geom_segment(aes(xend = grp, 
                   y = 0, 
                   yend = intercept_variation, 
                   group = group)) +
  coord_flip() +
  theme_classic() +
  labs(x = NULL,
       y = "Difference in intercept") +
  geom_hline(yintercept = 0, linetype = "dashed", col = "grey")

# random effects ecozone

ranef(betabinom_1o1f.site)

intercepts <- ranef(betabinom_1o1f.site) %>% 
  as_tibble() %>%
  subset(grpvar== "sitecond")%>%
  mutate(group = paste("ecoregion", grp)) %>%
  rename(intercept_variation = condval)


intercepts_eco <- ggplot(intercepts, 
                         aes(x = group, 
                             y = intercept_variation)) +
  geom_point() +
  geom_segment(aes(xend = group, 
                   y = 0, 
                   yend = intercept_variation, 
                   group = group)) +
  coord_flip() +
  theme_classic() +
  labs(x = NULL,
       y = "Difference in intercept") +
  geom_hline(yintercept = 0, linetype = "dashed", col = "grey")

ggsave("data/processed/modelling/intercepts.ecoregion.tiff", intercepts_eco , width = 6, height = 7)



#################################################################################
# --- check disturbance rate model per ecoregion ---
# check weather a specific ecoregion with high representation is dominating the global model
#################################################################################

betabinom_1o1f.site <- glmmTMB(cbind(dist.cells, no.dist) ~ 1+ Research + (1|year) + (1|subclass) + (1|sitecond), 
                          data = dist.df_1o1f_long, 
                          family = betabinomial(link = "logit"))

summary(betabinom_1o1f.site)
fixed_effects <- fixef(betabinom_1o1f)
1-exp(fixed_effects$cond[2]) 


# check effect of different ecoregions:

#---- subset data to different sites and store model results

dist.df_1o1f_long$sitecond <- as.numeric(dist.df_1o1f_long$sitecond)
ecoreg <- unique(dist.df_1o1f_long$sitecond)
ecoreg_df <- data.frame(prob_research = numeric(),
                         lower_bound = numeric(),
                         upper_bound = numeric(),
                        ecoregion = numeric())


for (i in ecoreg) {
  
  print(i)
  
  dist.df.sub <- dist.df_1o1f_long[dist.df_1o1f_long$sitecond == i,] #dist.df_1oManyf_long dist.df_1o1f_long
  
  model <- glmmTMB(cbind(dist.cells, no.dist) ~ 1+ Research + (1|year) + (1|subclass), 
                               data = dist.df.sub, 
                               family = betabinomial(link = "logit"))

  
  # Extracting the estimate and standard error for Research1
  estimate_research <- summary(model)$coefficients$cond["Research1", "Estimate"]
  std_error_research <- summary(model)$coefficients$cond["Research1", "Std. Error"]
  
  # Calculate the 95% confidence interval
  critical_value <- 1.96  # For a 95% confidence interval
  lower_bound <- estimate_research - critical_value * std_error_research
  upper_bound <- estimate_research + critical_value * std_error_research
  
  # percent more or less disturbances in managed forests
  prob_research <- 1-exp(estimate_research)
  
  # Converting logit confidence intervals to probabilities (in percent)
  lower_research <- 1-exp(lower_bound)
  upper_research <- 1-exp(upper_bound)
  
  # Store in an array
  research_values <- c(prob_research, lower_research, upper_research ,i)
  
  # Append as a new row to the dataframe
  ecoreg_df[nrow(ecoreg_df) + 1, ] <- research_values
  
}

# set all estimates where model didn't converge(1) = 0

cols_to_check <- c("prob_research", "lower_bound", "upper_bound")

# Find rows where any of the specified columns contain NaN or -Inf
rows_to_modify <- apply(ecoreg_df[cols_to_check], 1, function(x) any(is.nan(x) | x == -Inf))

# Set values to NA in these rows and specified columns
ecoreg_df[rows_to_modify, cols_to_check] <- NA

saveRDS(ecoreg_df, "data/processed/modelling/effects_ecoregion.rds")

#load ecoregion
wgb <- vect("data/processed/site_cond/wuchsgebiete_dissolved_counts.gpkg")

# merge with effects df
wgb.effect <- merge(wgb, ecoreg_df, by.x = "wg_bu", by.y = "ecoregion", all.x=T)

# store ecoregion with management effect
writeVector(wgb.effect, "data/processed/modelling/wuchsgebiete_dissolved_counts_effects_1o1.gpkg")




#################################################################################
# --- Simulation based model check --- 1:1 match models
#################################################################################

model <- fit_freq
#model <- fit_den
#model <- fit_size

original_data <- model_df 

# Simulate the data

set.seed(96)  # For reproducibility
simulated_data <- simulate(model, nsim = 100) 
simulated_data <- as.data.frame(simulated_data)

# add info to simulations
simulated_data$Research <- as.factor(original_data$Research)

# Reshape the data frame into long format
disturbance.sim_long <- pivot_longer(simulated_data, 
                                     cols = -c( Research ), 
                                     names_to = "simulation", 
                                     values_to = "value")

head(disturbance.sim_long)


# Calculate the metrics for each simulation

simulated_data <- disturbance.sim_long %>%
  group_by(simulation,Research) %>%
  summarize(mean_value = mean(value, na.rm = TRUE),
            min_value = min(value, na.rm = TRUE),
            max_value = max(value, na.rm = TRUE),
            sd = sd(value),
            n_zero = sum(value ==0)) %>%
  ungroup()


# calculate mean, min and max for the observed data!
# need to change to the respective metric in the summarize call!!!
observed_data <- original_data %>%
  group_by(Research) %>%
  summarize(mean_value = mean(mean.patch),# frequency density  max.patch mean.patch high.severity.rate 
            min_value = min(mean.patch), 
            max_value = max(mean.patch), 
            sd = sd(mean.patch),
            n_zero = sum(mean.patch ==0))


# add observation information to df for merging
observed_data$simulation <- "observed"
observed_data$mode <- "observed"

simulated_data$mode <- "simulated"

# Get the column order from observed_data
column_order <- names(observed_data)

# Rearrange the columns of simulated_data to match observed_data
simulated_data <- simulated_data %>%
  select(all_of(column_order))

# merge both dfs 
summary_obs_sim <- rbind(simulated_data, observed_data)

#create and reverse factor levels for plotting (so "observed is the last level)
summary_obs_sim$simulation <- as.factor(summary_obs_sim$simulation)
summary_obs_sim$simulation  <- factor(summary_obs_sim$simulation , levels = rev(levels(summary_obs_sim$simulation)))

# recode Reserach
summary_obs_sim$Research <- recode(as.factor(summary_obs_sim$Research),
                                   "1" = "Unmanaged",
                                   "0" = "Managed")
# Reverse factor levels
summary_obs_sim$Research <- factor(summary_obs_sim$Research, levels = rev(levels(summary_obs_sim$Research)))

summary_obs_sim_long <- pivot_longer(summary_obs_sim, cols = c(mean_value, min_value, sd, max_value,n_zero), names_to = "variable", values_to = "value")

My_Theme = theme(
  axis.title.x = element_text(size = 30),
  axis.text.x = element_text(size = 30),
  axis.text.y = element_text(size = 30),
  axis.title.y = element_text(size = 40),
  strip.text.x = element_text(size = 40),
  legend.key.height = unit(1, 'cm'),
  legend.key.width = unit(2, 'cm'),
  legend.title = element_text(size=40),
  legend.text = element_text(size=40))

plot <- ggplot(summary_obs_sim_long, aes(x = Research, y= value, fill=mode, color= mode)) +
  geom_boxplot() +
  labs(x = " ", y = "Value", fill="", color="") +
  scale_fill_manual(values=met.brewer("Kandinsky"))+
  scale_color_manual(values=met.brewer("Kandinsky"))+
  theme_minimal() +
  theme(legend.position = "top") +
  facet_wrap(~variable, scales = "free_y")+
  My_Theme
plot


# --- model check severity  ---

binom_model <- fit_severity
original_data <- dist.df_1o1f_long.sub.dist

nsim <- 100

# Simulate the data

set.seed(96)  # For reproducibility
simulated_data <- simulate(binom_model, nsim = nsim)  
simulated_data$Research <- as.factor(original_data$Research)


# ----- calculate disturbance rate

# Create an empty list to store results for each simulation
results_list <- vector("list", length = nsim)

# Each data frame contains columns 'dist.cells' and 'no.dist' for each simulation

num_simulations <- nsim
num_rows <-nrow(simulated_data)

# Initialize a matrix to store disturbance rate values for each row in each simulation
severity_rates <- matrix(0, nrow = num_rows, ncol = num_simulations)

for (sim_idx in 1:num_simulations) {
  sim <- simulated_data[[sim_idx]]
  
  # extract dist and no dist cells
  dist_severe <- sim[, 1]
  no_severe <- sim[, 2]
  
  for (row_idx in 1:num_rows) {
    #calculate annual disturbance rates
    severity_rate <- dist_severe[row_idx] / (dist_severe[row_idx] + no_severe[row_idx])
    severity_rates[row_idx, sim_idx] <- severity_rate
  }
}

disturbance.sim <- as.data.frame(severity_rates)
disturbance.sim$Research <- simulated_data$Research


# Reshape the data frame into long format
disturbance.sim_long <- pivot_longer(disturbance.sim, 
                                     cols = -c(Research ), 
                                     names_to = "simulation", 
                                     values_to = "value")

# Check the first few rows of the reshaped data frame
head(disturbance.sim_long)


# Calculate the metrics per year for each simulation

simulated_data <- disturbance.sim_long %>%
  group_by(simulation, Research) %>%
  summarize(mean_value = mean(value, na.rm = TRUE),
            min_value = min(value, na.rm = TRUE),
            max_value = max(value, na.rm = TRUE),
            sd = sd(value),
            n_zero = sum(value ==0)) %>%
  ungroup()


# ---- add observed values 

# calculate mean, min and max disturbance rate per year
observed_data <- original_data %>%
  group_by(Research) %>%
  summarize(mean_value = mean(high.severity.rate, na.rm = TRUE),
            min_value = min(high.severity.rate, na.rm = TRUE),
            max_value = max(high.severity.rate, na.rm = TRUE),
            sd = sd(high.severity.rate),
            n_zero = sum(high.severity.rate ==0))



# add observation information to df for merging
observed_data$simulation <- "observed"
observed_data$mode <- "observed"

simulated_data$mode <- "simulated"

# Get the column order from observed_data
column_order <- names(observed_data)

# Rearrange the columns of simulated_data to match observed_data
simulated_data <- simulated_data %>%
  select(all_of(column_order))

# merge both dfs 
summary_obs_sim <- rbind(simulated_data, observed_data)

#create and reverse factor levels for plotting (so "observed is the last level)
summary_obs_sim$simulation <- as.factor(summary_obs_sim$simulation)
summary_obs_sim$simulation  <- factor(summary_obs_sim$simulation , levels = rev(levels(summary_obs_sim$simulation)))

# recode Reserach
summary_obs_sim$Research <- recode(as.factor(summary_obs_sim$Research),
                                   "1" = "Unmanaged",
                                   "0" = "Managed")
# Reverse factor levels
summary_obs_sim$Research <- factor(summary_obs_sim$Research, levels = rev(levels(summary_obs_sim$Research)))

# Create a vector of colors for each level of V
v_colors <- c(rep("grey", length(unique(summary_obs_sim$simulation)) - 1), "red")


# density plot
ggplot(summary_obs_sim, aes(x = mean_value, group = simulation, color = simulation)) + 
  geom_density() +
  labs(x = "High severity rate", y = "Density") +
  theme_minimal() +
  scale_color_manual(values = v_colors, ) + 
  theme(legend.position = "none")


summary_obs_sim_long <- pivot_longer(summary_obs_sim, cols = c(mean_value, min_value, sd, max_value,n_zero), names_to = "variable", values_to = "value")


plot <- ggplot(summary_obs_sim_long, aes(x = Research, y= value, fill=mode, color= mode)) +
  geom_boxplot() +
  labs(x = " ", y = "Value", fill="", color="") +
  scale_fill_manual(values=met.brewer("Kandinsky"))+
  scale_color_manual(values=met.brewer("Kandinsky"))+
  theme_minimal() +
  theme(legend.position = "top") +
  facet_wrap(~variable, scales = "free_y") +
  My_Theme
plot

