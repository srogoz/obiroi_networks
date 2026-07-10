library(igraph)
library(intergraph)
library(networkLite)
library(network)
library(networkDynamic)
library(dplyr)
library(tidyr)
#plot prevalence
library(ggplot2)
library(tidyr)
library(ggrepel)
library(ggnewscale) #add colorscales

setwd("Desktop/EXP1_analysis")
source("global_parameters_exp1_inf.R")
antlevel_parameters<-c("strength_mean", "centrality", "burstiness", "local_efficiency", "waitingtime_mean", "clustering_local", "interactionlength_mean")
antlevel_parameters_nospace<-c("strengthmean", "centrality", "burstiness", "localefficiency", "waitingtimemean","waitingtimesd", "clusteringlocal", "interactionlength")
output_folder<-"network_parameters/exp1_inf/20min/10s_3min/1x/ANTLEVEL/"

#get files after infection
antlevel_files_after<-c("network_parameters/exp1_inf/20min/10s_3min/1x/ANTLEVEL_allintervalls_strength_mean24062026.csv",
                  "network_parameters/exp1_inf/20min/10s_3min/1x/ANTLEVEL_allintervalls_waitingtime_mean24062026.csv",
                  "network_parameters/exp1_inf/20min/10s_3min/1x/ANTLEVEL_allintervalls_waitingtime_sd24062026.csv",
                  "network_parameters/exp1_inf/20min/10s_3min/1x/ANTLEVEL_allintervalls_burstiness24062026.csv",
                  "network_parameters/exp1_inf/20min/10s_3min/1x/ANTLEVEL_allintervalls_interactionlength_mean24062026.csv",
                  "network_parameters/exp1_inf/20min/10s_3min/1x/ANTLEVEL_allintervalls_clustering_local24062026.csv",
                   "network_parameters/exp1_inf/20min/10s_3min/1x/ANTLEVEL_allintervalls_local_efficiency24062026.csv",
                  "network_parameters/exp1_inf/20min/10s_3min/1x/ANTLEVEL_allintervalls_centrality24062026.csv"
              )

antlevel_data_unmerged<-lapply(antlevel_files_after, read.csv)

#timepoint marker
names(antlevel_data_unmerged[[1]])<-paste0(names(antlevel_data_unmerged[[1]]), "_strengthmean")
names(antlevel_data_unmerged[[2]])<-paste0(names(antlevel_data_unmerged[[2]]), "_waitingtimemean")
names(antlevel_data_unmerged[[3]])<-paste0(names(antlevel_data_unmerged[[3]]), "_waitingtimesd")
names(antlevel_data_unmerged[[4]])<-paste0(names(antlevel_data_unmerged[[4]]), "_burstiness")
names(antlevel_data_unmerged[[5]])<-paste0(names(antlevel_data_unmerged[[5]]), "_interactionlength")
names(antlevel_data_unmerged[[6]])<-paste0(names(antlevel_data_unmerged[[6]]), "_clusteringlocal")
names(antlevel_data_unmerged[[7]])<-paste0(names(antlevel_data_unmerged[[7]]), "_localefficiency")
names(antlevel_data_unmerged[[8]])<-paste0(names(antlevel_data_unmerged[[8]]), "_centrality")




# #adding ant level
# antlevel_data_unmerged<-lapply(antlevel_data_unmerged, function(df){
#   df$ant<-expected_ants
#   return(df)
# })


antlevel_data_merged<-bind_cols(antlevel_data_unmerged)
antlevel_data_merged$ant<-expected_ants

antlevel_data <- antlevel_data_merged %>%
  pivot_longer(
    cols = -ant,
    names_to = "variable",
    values_to = "value"
  ) %>%
  separate(
    variable,
    into = c("colony", "time_interval", "parameter"),
    sep = "_"
  ) %>%
  mutate(
    time_interval = as.numeric(time_interval),
    treatment = sub("\\d+.*$", "", colony),
    
    genotype = case_when(
      treatment == "ba" ~ unlist(genotype_match[ant]),
      treatment == "a"  ~ "a",
      treatment == "b"  ~ "b"
    ),
    
    puremixed = case_when(
      treatment == "ba" ~ "m",
      TRUE              ~ "p"
    ),
    
    infection = if_else(
      ant == exposed_ants_dot[colony],
      "exposed",
      "nestmate"
    )
  )


#########################
#plot values for all ants
#########################
plots<-list()
for (param in antlevel_parameters_nospace){

plots[[param]]<-antlevel_data %>%
  filter(parameter == param) %>%
  ggplot(aes(x = colony, y = value)) +
  
  geom_boxplot(
    aes(color = treatment,
        fill = treatment),
    width = 0.5,
    outlier.shape = NA,
    alpha = 0.6
  ) +
  scale_color_manual(values = treatment_colors) +
  scale_fill_manual(values = treatment_colors) +
  ggnewscale::new_scale_color() +
  ggnewscale::new_scale_fill() +
  geom_point(
    aes(
      fill = genotype,
      color = infection
    ),
    shape = 21,
    size = 2,
    stroke = 1.5,
    alpha = 1,
    position = position_jitter(width = 0.12)
  ) +
  
  scale_fill_manual(values = anttypes_colors) +
  scale_color_manual(values = infection_status_ants) +
  
  labs(
    x = "Colony",
    y = paste0(param)
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
}


ALL_plots<-wrap_plots(plots, ncol = 2)+plot_layout(guides = 'collect', axes = "collect")
#highres 24, lowres 16, even_lower res 12
ggsave(
  filename = paste0(output_folder, "antlevelparameters_infection.png"),
  plot = ALL_boxplots,
  width = 14,
  height = 18,
  dpi = 300
)

##########
#before infection!!
#get files after infection
antlevel_files_after<-c("network_parameters/exp1_beforefungus/5_min/121_135/10s_3min_20minx1/ANTLEVEL_allintervalls_strength_mean23062026.csv", #STRENGTH MEAN
                        "network_parameters/exp1_beforefungus/5_min/121_135/10s_3min_20minx1/ANTLEVEL_allintervalls_waitingtime_mean23062026.csv", #wAITINGTIME MEAN
                        "network_parameters/exp1_beforefungus/5_min/121_135/10s_3min_20minx1/ANTLEVEL_allintervalls_waitingtime_sd23062026.csv",                       
                        "network_parameters/exp1_beforefungus/5_min/121_135/10s_3min_20minx1/ANTLEVEL_allintervalls_burstiness23062026.csv",
                        "network_parameters/exp1_beforefungus/5_min/121_135/10s_3min_20minx1/ANTLEVEL_allintervalls_interactionlength_mean23062026.csv",                       
                        "network_parameters/exp1_beforefungus/5_min/121_135/10s_3min_20minx1/ANTLEVEL_allintervalls_clustering_local23062026.csv",
                        "network_parameters/exp1_beforefungus/5_min/121_135/10s_3min_20minx1/ANTLEVEL_allintervalls_local_efficiency23062026.csv",
                        "network_parameters/exp1_beforefungus/5_min/121_135/10s_3min_20minx1/ANTLEVEL_allintervalls_centrality23062026.csv")

antlevel_data_unmerged<-lapply(antlevel_files_after, read.csv)

#timepoint marker
names(antlevel_data_unmerged[[1]])<-paste0(names(antlevel_data_unmerged[[1]]), "_strengthmean")
names(antlevel_data_unmerged[[2]])<-paste0(names(antlevel_data_unmerged[[2]]), "_waitingtimemean")
names(antlevel_data_unmerged[[3]])<-paste0(names(antlevel_data_unmerged[[3]]), "_waitingtimesd")
names(antlevel_data_unmerged[[4]])<-paste0(names(antlevel_data_unmerged[[4]]), "_burstiness")
names(antlevel_data_unmerged[[5]])<-paste0(names(antlevel_data_unmerged[[5]]), "_interactionlength")
names(antlevel_data_unmerged[[6]])<-paste0(names(antlevel_data_unmerged[[6]]), "_clusteringlocal")
names(antlevel_data_unmerged[[7]])<-paste0(names(antlevel_data_unmerged[[7]]), "_localefficiency")
names(antlevel_data_unmerged[[8]])<-paste0(names(antlevel_data_unmerged[[8]]), "_centrality")




# #adding ant level
# antlevel_data_unmerged<-lapply(antlevel_data_unmerged, function(df){
#   df$ant<-expected_ants
#   return(df)
# })


antlevel_data_merged<-bind_cols(antlevel_data_unmerged)
antlevel_data_merged$ant<-expected_ants

antlevel_data <- antlevel_data_merged %>%
  pivot_longer(
    cols = -ant,
    names_to = "variable",
    values_to = "value"
  ) %>%
  separate(
    variable,
    into = c("colony", "time_interval", "parameter"),
    sep = "_"
  ) %>%
  mutate(
    time_interval = as.numeric(time_interval),
    treatment = sub("\\d+.*$", "", colony),
    
    genotype = case_when(
      treatment == "ba" ~ unlist(genotype_match[ant]),
      treatment == "a"  ~ "a",
      treatment == "b"  ~ "b"
    ),
    
    puremixed = case_when(
      treatment == "ba" ~ "m",
      TRUE              ~ "p"
    ),
    
    infection = if_else(
      ant == exposed_ants_dot[colony],
      "exposed",
      "nestmate"
    )
  )


#########################
#plot values for all ants
#########################

plots<-list()
for (param in antlevel_parameters_nospace){
  
  plots[[param]]<-antlevel_data %>%
    filter(parameter == param) %>%
    ggplot(aes(x = colony, y = value)) +
    
    geom_boxplot(
      aes(color = treatment,
          fill = treatment),
      width = 0.5,
      outlier.shape = NA,
      alpha = 0.6
    ) +
    scale_color_manual(values = treatment_colors) +
    scale_fill_manual(values = treatment_colors) +
    ggnewscale::new_scale_color() +
    ggnewscale::new_scale_fill() +
    geom_point(
      aes(
        fill = genotype,
        color = infection
      ),
      shape = 21,
      size = 2,
      stroke = 1.5,
      alpha = 1,
      position = position_jitter(width = 0.12)
    ) +
    
    scale_fill_manual(values = anttypes_colors) +
    scale_color_manual(values = infection_status_ants) +
    
    labs(
      x = "Colony",
      y = paste0(param)
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}


ALL_plots<-wrap_plots(plots, ncol = 2)+plot_layout(guides = 'collect', axes = "collect")
#highres 24, lowres 16, even_lower res 12
ggsave(
  filename = paste0(output_folder, "antlevelparameters_preinfection.png"),
  plot = ALL_plots,
  width = 14,
  height = 18,
  dpi = 300
)


#######
#correlate speed of spread of colony
############################################
#GLMM on network parameter change 
############################################
mean_overtimes<-before_after_infection_merged%>%group_by(colony_name, globaltime,infection) %>%
  summarise(
    across(where(is.numeric), mean, na.rm = TRUE),
    .groups = "drop"
  )
#mean over before and after infection
mean_overtimes<-mean_overtimes%>%group_by(colony_name,infection) %>%
  summarise(
    across(where(is.numeric), mean, na.rm = TRUE),
    .groups = "drop"
  )

#long format for glmm delta_paramchange ~ treatment +(1|colony)

delta <- mean_overtimes %>%
  pivot_longer(
    cols = -c(colony_name, infection),
    names_to = "parameter",
    values_to = "mean_value"
  ) %>%
  pivot_wider(
    names_from = infection,
    values_from = mean_value
  ) %>%
  mutate(delta = before - after)

delta<-delta%>%mutate(
  treatment = sub("\\d+.*$", "", colony_name))
##NORMED delta vlaues
delta<-delta%>%mutate(delta_norm = delta/before)

delta_strength<-delta%>%filter(parameter=="strength_mean")
hist(delta_strength$delta)

##glmm for strength
###################################
model_delta_s<-glmmTMB(delta~treatment, 
                       data = delta_strength,
                       family = gaussian())

sim <- DHARMa::simulateResiduals(model_delta_s)
plot(sim)

pairs <- emmeans(model_delta_s, pairwise ~ treatment, type = "response")
stat_test <- pairs$contrasts
stat_test
#normed
hist(delta_strength$delta_norm)
model_delta_snorm<-glmmTMB(delta_norm~treatment, 
                           data = delta_strength,
                           family = gaussian())

sim <- DHARMa::simulateResiduals(model_delta_snorm)
plot(sim)


pairs <- emmeans(model_delta_snorm, pairwise ~ treatment, type = "response")
stat_test <- pairs$contrasts
stat_test

###########################
#conserving all different timepoints before
###########################
mean_overtimes<-before_after_infection_merged%>%group_by(colony_name, globaltime,infection) %>%
  summarise(
    across(where(is.numeric), mean, na.rm = TRUE),
    .groups = "drop"
  )
#mean over before and after infection
mean_overtimes<-mean_overtimes%>%group_by(colony_name,infection, globaltime) %>%
  summarise(
    across(where(is.numeric), mean, na.rm = TRUE),
    .groups = "drop"
  )

#long format for glmm delta_paramchange ~ treatment +(1|colony)

delta <- mean_overtimes %>%
  pivot_longer(
    cols = -c(colony_name, infection, globaltime),
    names_to = "parameter",
    values_to = "mean_value"
  )%>%
  pivot_wider(
    names_from = infection,
    values_from = mean_value
  ) %>%
  mutate(delta = before - after)

delta<-delta%>%mutate(
  treatment = sub("\\d+.*$", "", colony_name))
##NORMED delta vlaues
delta<-delta%>%mutate(delta_norm = delta/before)

delta_strength<-delta%>%filter(parameter=="strength_mean")
hist(delta_strength$delta)

##glmm for strength
###################################
model_delta_s<-glmmTMB(delta~treatment, 
                       data = delta_strength,
                       family = gaussian())

sim <- DHARMa::simulateResiduals(model_delta_s)
plot(sim)

pairs <- emmeans(model_delta_s, pairwise ~ treatment, type = "response")
stat_test <- pairs$contrasts
stat_test
#normed
hist(delta_strength$delta_norm)
model_delta_snorm<-glmmTMB(delta_norm~treatment, 
                           data = delta_strength,
                           family = gaussian())

sim <- DHARMa::simulateResiduals(model_delta_snorm)
plot(sim)


pairs <- emmeans(model_delta_snorm, pairwise ~ treatment, type = "response")
stat_test <- pairs$contrasts
stat_test


#get antlevel parameters of last interval before infection, with interaction length filter >10s
antlevel_files<-list()
antlevel_files[["strength_mean"]]<-c("network_parameters/exp1_beforefungus/5_min/121_135/10s_5min/ANTLEVEL_allintervalls_strength_mean16062026.csv",
                                     )

#mean over last 5 values
#get spreading speed from simulations in 10s interaction network
#overlay prevalences of all colonies in one treatment in one plot
prevalence_list<-get_colony_files("simulations/exp1_inf/10s_5min_5min/data/", selected_colonies)

prevalence_data_unmerged<-lapply(prevalence_list, readRDS)
prevalence_data_unmerged <- Map(function(df, name) {
  df %>%
    mutate(
      colony_name = name,
      treatment = str_extract(name, "^[^0-9]+")
    )
}, prevalence_data_unmerged, names(prevalence_data_unmerged))

prevalence_data_unmerged <- lapply(prevalence_data_unmerged, function(df) {
  df %>%
    mutate(
      treatment = str_extract(colony_name, "^[^0-9]+")
    )
})
#bind all into one df
prevalence_data_merged<-bind_rows(prevalence_data_unmerged)

saveRDS( prevalence_data_merged, paste0("simulations/exp1_inf/10s_5min_5min/data/prevalenceALL_10s.rds"))



#plot all in one mean
ggplot(prevalence_data_merged, 
       aes(x = time, y = prevalence, color = treatment)) +
  stat_summary(fun = mean, geom = "line", size = 1.2) +
  scale_color_manual(values = treatment_colors) +
  xlim(0,5000)
theme_minimal()

#############
#mean over colony and treatment
#############
#mean over colony, n = 16 seedants
prevalence_mean_colony<-prevalence_data_merged|>
  group_by(colony_name, time, treatment, genotype) |>
  summarise(
    mean_prev = mean(prevalence, na.rm = TRUE),
    sd_val   = sd(prevalence, na.rm = TRUE),
    n             = n(),
    se_val   = sd_val / sqrt(n),
    .groups = "drop"
  )



saveRDS(prevalence_mean_colony, paste0("simulations/colony_mean_ALL_2s.rds"))

#mean over treatment, n = 8 replicates per colony
prevalence_mean_treatment<-prevalence_mean_colony|>
  group_by(treatment, time) |>
  summarise(
    mean_preva = mean(mean_prev, na.rm = TRUE),
    n             = n(),
    sd_val   = sd(mean_prev, na.rm = TRUE),
    se_val   = sd_val / sqrt(n),
    .groups = "drop"
  )


saveRDS(prevalence_mean_treatment, paste0("simulations/treatment_mean_ALL_2s.rds"))

ggplot(prevalence_mean_treatment, 
       aes(x = time, y = mean_preva, color = treatment,
           fill  = treatment)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(
      ymin = mean_preva - se_val,
      ymax = mean_preva + se_val
    ),
    alpha = 0.25,
    color = NA
  ) +
  scale_color_manual(
    values = treatment_colors
  ) +
  scale_fill_manual(
    values = treatment_colors
  ) +
  xlim(0,12000)+
  labs(x = "time [frame]", y = " mean prevalence",
       color = "treatment", fill = "treatment",
       title = " deterministic SI-spread interaction(>10s) network  ") +
  theme_minimal()

