source("global_parameters_exp2_12hs.R")
antlevel_parameters<-c("strength_mean", "centrality", "burstiness", "local_efficiency", "waitingtime_mean", "clustering_local", "interactionlength_mean")


antlevel_files<-list()
antlevel_files[["strength_mean"]]<-c("network_parameters/exp2_12hs/5_mins/1_15/ANTLEVEL_allintervalls_strength_mean_10032026.csv",
                  "network_parameters/exp2_12hs/5_mins/16_30/ANTLEVEL_allintervalls_strength_mean10032026.csv",
                  "network_parameters/exp2_12hs/5_mins/106_120/ANTLEVEL_allintervalls_strength_mean10032026.csv",
                  "network_parameters/exp2_12hs/5_mins/121_135/ANTLEVEL_allintervalls_strength_mean10032026.csv",
                  "network_parameters/exp2_inf/5mins/ANTLEVEL_allintervalls_strength_mean10032026.csv" #infection
)

antlevel_files[["centrality"]]<-c("network_parameters/exp2_12hs/5_mins/1_15/ANTLEVEL_allintervalls_centrality_10032026.csv",
                                     "network_parameters/exp2_12hs/5_mins/16_30/ANTLEVEL_allintervalls_centrality10032026.csv",
                                     "network_parameters/exp2_12hs/5_mins/106_120/ANTLEVEL_allintervalls_centrality10032026.csv",
                                     "network_parameters/exp2_12hs/5_mins/121_135/ANTLEVEL_allintervalls_centrality10032026.csv",
                                     "network_parameters/exp2_inf/5mins/ANTLEVEL_allintervalls_centrality10032026.csv" #infection
)

antlevel_files[["clustering_local"]]<-c("network_parameters/exp2_12hs/5_mins/1_15/ANTLEVEL_allintervalls_clustering_local_10032026.csv",
                                  "network_parameters/exp2_12hs/5_mins/16_30/ANTLEVEL_allintervalls_clustering_local10032026.csv",
                                  "network_parameters/exp2_12hs/5_mins/106_120/ANTLEVEL_allintervalls_clustering_local10032026.csv",
                                  "network_parameters/exp2_12hs/5_mins/121_135/ANTLEVEL_allintervalls_clustering_local10032026.csv",
                                  "network_parameters/exp2_inf/5mins/ANTLEVEL_allintervalls_clustering_local10032026.csv" #infection
)

antlevel_files[["local_efficiency"]]<-c("network_parameters/exp2_12hs/5_mins/1_15/ANTLEVEL_allintervalls_local_efficiency_10032026.csv",
                                        "network_parameters/exp2_12hs/5_mins/16_30/ANTLEVEL_allintervalls_local_efficiency10032026.csv",
                                        "network_parameters/exp2_12hs/5_mins/106_120/ANTLEVEL_allintervalls_local_efficiency10032026.csv",
                                        "network_parameters/exp2_12hs/5_mins/121_135/ANTLEVEL_allintervalls_local_efficiency10032026.csv",
                                        "network_parameters/exp2_inf/5mins/ANTLEVEL_allintervalls_local_efficiency10032026.csv" #infection
)

antlevel_files[["waitingtime_mean"]]<-c("network_parameters/exp2_12hs/5_mins/1_15/ANTLEVEL_allintervalls_waitingtime_mean_10032026.csv",
                                        "network_parameters/exp2_12hs/5_mins/16_30/ANTLEVEL_allintervalls_waitingtime_mean10032026.csv",
                                        "network_parameters/exp2_12hs/5_mins/106_120/ANTLEVEL_allintervalls_waitingtime_mean10032026.csv",
                                        "network_parameters/exp2_12hs/5_mins/121_135/ANTLEVEL_allintervalls_waitingtime_mean10032026.csv",
                                        "network_parameters/exp2_inf/5mins/ANTLEVEL_allintervalls_waitingtime_mean10032026.csv" #infection
)

antlevel_files[["burstiness"]]<-c("network_parameters/exp2_12hs/5_mins/1_15/ANTLEVEL_allintervalls_burstiness_10032026.csv",
                                        "network_parameters/exp2_12hs/5_mins/16_30/ANTLEVEL_allintervalls_burstiness10032026.csv",
                                        "network_parameters/exp2_12hs/5_mins/106_120/ANTLEVEL_allintervalls_burstiness10032026.csv",
                                        "network_parameters/exp2_12hs/5_mins/121_135/ANTLEVEL_allintervalls_burstiness10032026.csv",
                                        "network_parameters/exp2_inf/5mins/ANTLEVEL_allintervalls_burstiness10032026.csv" #infection
)

antlevel_files[["interactionlength_mean"]]<-c("network_parameters/exp2_12hs/5_mins/1_15/ANTLEVEL_allintervalls_interactionlength_mean_10032026.csv",
                                        "network_parameters/exp2_12hs/5_mins/16_30/ANTLEVEL_allintervalls_interactionlength_mean10032026.csv",
                                        "network_parameters/exp2_12hs/5_mins/106_120/ANTLEVEL_allintervalls_interactionlength_mean10032026.csv",
                                        "network_parameters/exp2_12hs/5_mins/121_135/ANTLEVEL_allintervalls_interactionlength_mean10032026.csv",
                                        "network_parameters/exp2_inf/5mins/ANTLEVEL_allintervalls_interactionlength_mean10032026.csv" #infection
)

antlevel_data_list<-list()

for (param in antlevel_parameters){

antlevel_data_unmerged<-lapply(antlevel_files[[param]], read.csv)
#timepoint marker
names(antlevel_data_unmerged[[1]])<-paste0(names(antlevel_data_unmerged[[1]]), "_0-1.25h")
names(antlevel_data_unmerged[[2]])<-paste0(names(antlevel_data_unmerged[[2]]), "_1.25-2.5h")
names(antlevel_data_unmerged[[3]])<-paste0(names(antlevel_data_unmerged[[3]]), "_8.75-10h")
names(antlevel_data_unmerged[[4]])<-paste0(names(antlevel_data_unmerged[[4]]), "_10-11.25h")
names(antlevel_data_unmerged[[5]])<-paste0(names(antlevel_data_unmerged[[5]]), "_fungal exp")

#adding ant level
antlevel_data_unmerged<-lapply(antlevel_data_unmerged, function(df){
  df$ant<-expected_ants
  return(df)
})

antlevel_data_merged<-bind_rows(antlevel_data_unmerged)

#all data for all timepoints including infection
antlevel_data<- antlevel_data_merged%>%
  pivot_longer(
    cols = -ant,
    names_to = "colony_name",
    values_to = param
  )%>%
  separate(colony_name,
           into = c("colony_name", "time_interval", "globaltime"),
           sep = "_",
           remove = TRUE)%>%
  mutate(time_interval = as.numeric(time_interval))%>%
  mutate(
    treatment = sub("\\d+.*$", "", colony_name)   # extract 'a', 'b', or 'ba'
  )%>%
  mutate(
    anttype = case_when(
      treatment == "bA" ~ unlist(anttype_match[["bA"]][ant]),
      treatment == "bB" ~ unlist(anttype_match[["bB"]][ant]),# lookup in genotype_match
      treatment == "B" ~ "B",                            # all ants in pure B colony
      treatment == "b" ~ "b"                             # all ants in pure b colony
    )
  )%>%
  mutate(
    puremixed = case_when(
      treatment == "bA" ~ "m",
      treatment == "bB" ~ "m",
      treatment == "B" ~ "p",                           
      treatment == "b" ~ "p"                             
    )
  )
antlevel_data_list[[param]]<-antlevel_data
}

identical(antlevel_data_list[["centrality"]]$ant, antlevel_data_list[["burstiness"]]$ant)
identical(antlevel_data_list[["centrality"]]$ant, antlevel_data_list[["strength_mean"]]$ant)
identical(antlevel_data_list[["centrality"]]$ant, antlevel_data_list[["local_efficiency"]]$ant)
identical(antlevel_data_list[["centrality"]]$ant, antlevel_data_list[["waitingtime_mean"]]$ant)
identical(antlevel_data_list[["centrality"]]$ant, antlevel_data_list[["interactionlength_mean"]]$ant)
identical(antlevel_data_list[["centrality"]]$ant, antlevel_data_list[["clustering_local"]]$ant)

identical(antlevel_data_list[["centrality"]]$anttype, antlevel_data_list[["burstiness"]]$anttype)
identical(antlevel_data_list[["centrality"]]$anttype, antlevel_data_list[["strength_mean"]]$anttype)
identical(antlevel_data_list[["centrality"]]$anttype, antlevel_data_list[["local_efficiency"]]$anttype)
identical(antlevel_data_list[["centrality"]]$anttype, antlevel_data_list[["waitingtime_mean"]]$anttype)
identical(antlevel_data_list[["centrality"]]$anttype, antlevel_data_list[["interactionlength_mean"]]$anttype)
identical(antlevel_data_list[["centrality"]]$anttype, antlevel_data_list[["clustering_local"]]$anttype)

identical(antlevel_data_list[["centrality"]]$colony_name, antlevel_data_list[["burstiness"]]$colony_name)
identical(antlevel_data_list[["centrality"]]$colony_name, antlevel_data_list[["strength_mean"]]$colony_name)
identical(antlevel_data_list[["centrality"]]$colony_name, antlevel_data_list[["local_efficiency"]]$colony_name)
identical(antlevel_data_list[["centrality"]]$colony_name, antlevel_data_list[["waitingtime_mean"]]$colony_name)
identical(antlevel_data_list[["centrality"]]$colony_name, antlevel_data_list[["interactionlength_mean"]]$colony_name)
identical(antlevel_data_list[["centrality"]]$colony_name, antlevel_data_list[["clustering_local"]]$colony_name)


antlevel_data<-antlevel_data_list[["strength_mean"]]
antlevel_data$burstiness<-antlevel_data_list[["burstiness"]]$burstiness
antlevel_data$local_efficiency<-antlevel_data_list[["local_efficiency"]]$local_efficiency
antlevel_data$clustering_local<-antlevel_data_list[["clustering_local"]]$clustering_local
antlevel_data$waitingtime_mean<-antlevel_data_list[["waitingtime_mean"]]$waitingtime_mean
antlevel_data$interactionlength_mean<-antlevel_data_list[["interactionlength_mean"]]$interactionlength_mean
antlevel_data$centrality<-antlevel_data_list[["centrality"]]$centrality



#antleveldata mean over all timeintervals
#only first 30mins
antlevel_data<-antlevel_data%>%filter(time_interval<= "6")
antlevel_data_mean<-antlevel_data%>%group_by(colony_name, globaltime,ant, puremixed, anttype, treatment) %>%
  summarise(
    across(where(is.numeric), mean, na.rm = TRUE),
    .groups = "drop"
  )%>%select(-time_interval)

#plot colonies and single ants with anttypes
param<-"local_efficiency"

antlevel_data_mean_1<-antlevel_data_mean%>%filter(globaltime =="10-11.25h")
antlevel_data_mean_1<-antlevel_data_mean_1%>%filter(treatment =="bB")
ggplot(antlevel_data_mean_1, aes(x = colony_name, y = .data[[param]], fill = treatment)) +
  
  geom_jitter(aes(color = anttype), width = 0.1, alpha = 0.7, size = 2) +
  geom_boxplot(width = 0.3, alpha = 0.4) +
  theme_minimal() +
  scale_fill_manual(values = anttypes_colors) +
  scale_color_manual(values = treatment_colors) +
  labs(
    title = " ",
    x = "colony",
    y = param
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # rotate x labels for readability
# Save each plot as PNG
ant_plot<-ggplot(antlevel_data_mean,
                      aes(x = treatment,
                          y = .data[[param]],
                          color = anttype,
                          fill  = colony_name,
                          
                      )) +
  geom_boxplot(aes(fill = colony_name), alpha = 0.3)+
  
  geom_point(aes(color = anttype),size = 2)+

  # stat_summary(
  #   fun = mean,
  #   geom = "point",
  #   color = "seagreen"
  #   
  # )+
  scale_color_manual(
    values = anttypes_colors,
  ) +
  scale_fill_manual(
    values = treatment_colors
  ) +
  labs(
    x = "treatments",
    y = param,
    color = "anttype:",
    fill  = "treatment"
  ) +
  theme_minimal()

##
hist(infection_mean$interaction_length_mean)
model_lm<- glmmTMB(
  log(interaction_length_mean) ~ treatment,
  data = infection_mean, 
  #family = Gamma(link = log)# 'ensures is positive', right skewed (link = log)
  family = gaussian()
)

hist(antlevel_data_mean_1$strength_mean)
model_s<-glmmTMB(
log(strength_mean) ~ anttype + (1|colony_name) ,
#dispformula =~ 1 ,
 data = antlevel_data_mean_1,
family = gaussian()
#family = Gamma(link = "log")
  )

sim <- DHARMa::simulateResiduals(model_s)
plot(sim)

drop1(model_s, test = "Chisq")

emmeans(model_s, pairwise ~ anttype, type = "response")

#immediately move on to the deltas
