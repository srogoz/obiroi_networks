#load needed dependencies and functions
library(igraph)
library(intergraph)
library(networkLite)
library(network)
library(networkDynamic)
library(dplyr)
library(tidyr)
#plot prevalence
library(ggplot2)
library(stringr)

setwd("Desktop/EXP1_analysis")

#load initial parameters
source("global_parameters_exp2_inf.R")
source("function_collection.R")
source("simulation_functions.R")

#loop over all colonies
for (i in 1:length(selected_colonies)){
  
  
  colony_name<-selected_colonies[[i]]
  present_ants<-present_ants_list[[i]]
  cat("Running simulation for colony:", colony_name, "\n")
  framerate_col<-framerate[[colony_name]]
  lower_limit_interaction<-interaction_limit_sec*framerate_col
  higher_limit_interaction<-interaction_limit_min*60*framerate_col
  
  #intervall_list<-[]
  ###################
  
  #prepare network obj, interaction limits
  network_obj_raw <-readRDS(network_obj_list[[colony_name]])
  
  #pick time slice of network object
  missing_ants <- setdiff(expected_ants, present_ants)
  
  network_obj <- network_obj_raw[
    !network_obj_raw$head %in% missing_ants &
      !network_obj_raw$tail %in% missing_ants,
  ]
  
  number_intervalls<-round((12*60*60*framerate_col)/(interval_time*60*framerate_col),0)
  interval_size<-interval_time*60*framerate_col
  intervall_timepoints <- seq(
    from = 0,
    by   = interval_size,
    length.out = number_intervalls + 1
  )
  
  time_window_start<-intervall_timepoints[1]
  time_window_end<-intervall_timepoints[1+1]
  
  #limit to timeintervall
  network_obj_5<- network_obj[
    network_obj$onset < time_window_end &
      network_obj$terminus  < time_window_end,]
  
  started_but_not_ended<-sum(
    network_obj_5$onset<time_window_end &
      network_obj_5$terminus   > time_window_end)
  ##############
  #CONDITIONS
  ##############
  #truncate from both sides
  #network_obj_5 <-network_obj[(network_obj$duration>lower_limit_interaction & network_obj$duration<higher_limit_interaction),]
  #truncate and minimize from the high side. Still count interaction but does not keep counting after time
  # Cap durations at the upper limit
  network_obj_filter1 <-network_obj_5[(network_obj_5$duration>lower_limit_interaction),]
  network_obj_filter1$duration <- pmin(network_obj_filter1$duration, higher_limit_interaction)
  
  
  #transform into network dynamic object
  
  edge_df <- data.frame(onset=network_obj_filter1$onset, terminus=network_obj_filter1$terminus, tail=network_obj_filter1$tail, head=network_obj_filter1$head)
  #tail and head must be numeric
  
  #edge_df$tail<-match(edge_df$tail, present_ants)
  #edge_df$head<-match(edge_df$head, present_ants)
  #edge_df$onset<-as.integer(edge_df$onset)
  
  ##CHECK NUMBER OF ONE FRAME INTERACTIONS
  bad_edges <- filter(edge_df,(onset >= terminus))
  onset_diff <- diff(bad_edges$onset)
  ##transform into undirected network object
  edge_df_corrected <- edge_df %>% filter(onset < terminus)
  
  
  ##########write a function for simulating spread on a network
  #duration, whole length of simulation
  #spreading prob: 1
  #seed_ant: infected individual
  #edgelist structure: onset, terminus, head, tail
  
  edge_list<-edge_df_corrected
  
  
  
  #result <-spread_SI(edge_list=edge_df_corrected, spreading_prob=1, seed_ant = "PP", duration=50000)
  ##############
  ###for all ants in one colony
  ##############
  
  seedvariation_prevalence <- matrix(0, nrow = duration, ncol = length(present_ants))
  seedvariation_saturationtime<-matrix(0, nrow = 1, ncol = length(present_ants))
  colnames(seedvariation_saturationtime) <- present_ants
  colnames(seedvariation_prevalence) <- present_ants
  
  # Now loop over seeds
  for (ant in present_ants) {
    cat("Running simulation for seed:", ant, "\n")
    
    result <- spread_SI(edge_list = edge_df, 
                        spreading_prob = 1, 
                        seed_ant = ant, 
                        duration = duration,
                        present_ants = present_ants)
    
    # Fill this ant's column with its prevalence trajectory
    seedvariation_prevalence[, ant] <- result$prevalence
    seedvariation_saturationtime[ant] <- result$saturated
  }
  #convert to dataframe
  prevalence_df<-data.frame(time = 1:length(result$prevalence), seedvariation_prevalence)
  
  prevalence_long <- prevalence_df %>%
    pivot_longer(cols = -time, 
                 names_to = "seed_ant", 
                 values_to = "prevalence")
 
  prevalence_long<-prevalence_long%>%mutate(treatment = str_extract(colony_name, "^[^0-9]+"))
  
  prevalence_long<-prevalence_long%>%mutate(genotype =case_when(
    treatment == "bA" ~ unlist(anttype_match[["bA"]][seed_ant]),
    treatment == "bB" ~ unlist(anttype_match[["bB"]][seed_ant]),
    treatment == "b" ~ "b",                            # all ants in pure colony
    treatment == "B" ~ "B"     
  ))
  
  saveRDS( prevalence_long, paste0("simulations/exp2_inf/same_genotype_prevalence/",colony_name,"_prevalencesim_2s.rds"))
  
  max_duration_colony<-max(seedvariation_saturationtime)
  
  #all ants in one plot per colony
  prevalence_plot<-ggplot(prevalence_long, aes(x= time, y= prevalence, color= genotype, group = seed_ant))+
    geom_line(alpha = 0.8, size = 0.8) +
    scale_color_manual(values = anttypes_colors) +
    labs(x = "time", y = "prevalence",  
         title = paste0("SI-simulation: prevalence by seed ant\n colony: ", colony_name, "\n on min 3s interactions")) +
    xlim(0,max_duration_colony)+
    theme_minimal() +
    theme()
  
  ggsave(filename = paste0("simulations/exp2_inf/same_genotype_prevalence/" , 
                           colony_name , ".png"), plot = prevalence_plot)
  
  
}


#overlay prevalences of all colonies in one treatment in one plot
prevalence_list<-get_colony_files("simulations/exp2_inf/same_genotype_prevalence", selected_colonies)

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

saveRDS( prevalence_data_merged, paste0("simulations/exp2_inf/prevalenceALL_2s.rds"))

prevalence_data_merged<-readRDS("simulations/prevalenceALL_2s.rds")

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



saveRDS(prevalence_mean_colony, paste0("simulations/exp2_inf/colony_mean_ALL_2s.rds"))

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


saveRDS(prevalence_mean_treatment, paste0("simulations/exp2_inf/treatment_mean_ALL_2s.rds"))

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
  xlim(0,5000)+
  labs(x = "time [frame]", y = " mean prevalence",
       color = "treatment", fill = "treatment",
       title = " deterministic SI-model on 2s-long interaction network  ") +
  theme_minimal()


############
#check difference in prevalence between exposed and unexposed ant 
#mean over ants of different genotype
prevalence_mean_colony<-prevalence_data_merged|>
  group_by(colony_name, time, genotype, treatment) |> #add genotype of seedants in simulation to see if that has an effect in groups
  summarise(
    mean_prev = mean(prevalence, na.rm = TRUE),
    sd_val   = sd(prevalence, na.rm = TRUE),
    n             = n(),
    se_val   = sd_val / sqrt(n),
    .groups = "drop"
  )



saveRDS(prevalence_mean_colony, paste0("simulations/exp2_inf/split_genotype_prevalence/colony_mean_ALL_2s_byGENOTYPE.rds"))

#mean over treatment, n = 8 replicates per colony
prevalence_mean_gensplit<-prevalence_mean_colony|>
  group_by(treatment, genotype, time) |> #add genotype of seed - ant to see if that takes an effect in groups
  summarise(
    mean_preva = mean(mean_prev, na.rm = TRUE),
    n             = n(),
    sd_val   = sd(mean_prev, na.rm = TRUE),
    se_val   = sd_val / sqrt(n),
    .groups = "drop"
  )


saveRDS(prevalence_mean_gensplit, paste0("simulations/exp2_inf/split_genotype_prevalence/treatment_mean_ALL_2s_byGENOTYPE.rds"))

ggplot(prevalence_mean_gensplit, 
       aes(x = time, y = mean_preva, color = genotype,
           fill  = genotype)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(
      ymin = mean_preva - se_val,
      ymax = mean_preva + se_val
    ),
    alpha = 0.25,
    color = NA
  ) +
  facet_wrap(~ treatment) +  #subgraphs in different pannels for different genotype of seed_ant
  # scale_color_manual(
  #   values = treatment_colors
  # ) +
  # scale_fill_manual(
  #   values = treatment_colors
  # ) +
  xlim(0,5000)+
  labs(x = "time [frame]", y = " mean prevalence in genotype subgroups ",
       color = "subgroup", fill = "subgroup",
       title = " deterministic SI-sim on min 2s-long interaction network  ") +
  theme_minimal()


#########
#different seed ant on same genotype subgroup
ggplot(prevalence_mean_gensplit, 
       aes(x = time, y = mean_preva, 
           color = gen_split,
           fill = gen_split,
           group = genotype)) +
  
  geom_line(linewidth = 1) +
  
  geom_ribbon(
    aes(ymin = mean_preva - se_val,
        ymax = mean_preva + se_val),
    alpha = 0.25,
    color = NA
  ) +
  
  facet_wrap(~ gen_split) +  #sub pannels, with different levels of 
  
  scale_color_manual(values = treatment_colors) +
  scale_fill_manual(values = treatment_colors) +
  
  xlim(0, 5000) +
  
  labs(x = "time [frame]", 
       y = "mean prevalence",
       color = "subgroup",
       fill = "subgroup",
       title = "SI model: effect of seed genotype on subgroup prevalence") +
  
  theme_minimal()
#############
#BA
#############
# mixed colonies, spread differences in subgroups bB and bA
prevalence_list_A<-get_colony_files("simulations/exp2_inf/split_genotype_prevalence/A", bA_colonies)
prevalence_list_b<-get_colony_files("simulations/exp2_inf/split_genotype_prevalence/bfrombA", bA_colonies)

prevalence_data_unmerged_A<-lapply(prevalence_list_A, readRDS)
prevalence_data_unmerged_A <- Map(function(df, name) {
  df %>%
    mutate(
      colony_name = name,
      treatment = str_extract(name, "^[^0-9]+"),
      gen_split = "A",
     genotype = unlist(anttype_match[["bA"]][seed_ant])
    )
}, prevalence_data_unmerged_A, names(prevalence_data_unmerged_A))

prevalence_data_unmerged_A <- lapply(prevalence_data_unmerged_A, function(df) {
  df %>%
    mutate(
      treatment = str_extract(colony_name, "^[^0-9]+")
    )
})

prevalence_data_unmerged_b<-lapply(prevalence_list_b, readRDS)
prevalence_data_unmerged_b <- Map(function(df, name) {
  df %>%
    mutate(
      colony_name = name,
      treatment = str_extract(name, "^[^0-9]+"),
      gen_split = "b",
      genotype = unlist(anttype_match[["bA"]][seed_ant])
    )
}, prevalence_data_unmerged_b, names(prevalence_data_unmerged_b))

prevalence_data_unmerged_b <- lapply(prevalence_data_unmerged_b, function(df) {
  df %>%
    mutate(
      treatment = str_extract(colony_name, "^[^0-9]+")
    )
})


#bind all into one df
prevalence_data_merged_A<-bind_rows(prevalence_data_unmerged_A)
prevalence_data_merged_b<-bind_rows(prevalence_data_unmerged_b)

prevalence_data_merged_bA<-bind_rows(prevalence_data_merged_A, prevalence_data_merged_b)
saveRDS( prevalence_data_merged, paste0("simulations/exp2_inf/prevalenceALL_2s_anttype_split.rds"))

prevalence_data_merged_bA<-readRDS("simulations/prevalenceALL_2s_anttype_split.rds")

#plot all in one mean
ggplot(prevalence_data_merged_bA, 
       aes(x = time, y = prevalence, color = treatment)) +
  stat_summary(fun = mean, geom = "line", size = 1.2) +
  scale_color_manual(values = treatment_colors) +
  xlim(0,5000)
theme_minimal()

#############
#mean over colony and treatment
#############
#mean over colony, n = 16 seedants
prevalence_mean_colony<-prevalence_data_merged_bA|>
  group_by(colony_name, time, treatment, genotype, gen_split) |>
  summarise(
    mean_prev = mean(prevalence, na.rm = TRUE),
    sd_val   = sd(prevalence, na.rm = TRUE),
    n             = n(),
    se_val   = sd_val / sqrt(n),
    .groups = "drop"
  )



saveRDS(prevalence_mean_colony, paste0("simulations/exp2_inf/colony_mean_ALL_2s_anttypesplit.rds"))

#mean over treatment, n = 8 replicates per colony
prevalence_mean_treatment<-prevalence_mean_colony|>
  group_by(treatment, time, gen_split, genotype) |>
  summarise(
    mean_preva = mean(mean_prev, na.rm = TRUE),
    n             = n(),
    sd_val   = sd(mean_prev, na.rm = TRUE),
    se_val   = sd_val / sqrt(n),
    .groups = "drop"
  )


saveRDS(prevalence_mean_treatment, paste0("simulations/exp2_inf/treatment_mean_ALL_2s_anttypesplit.rds"))
prevalence_mean_treatment_bA<-readRDS("simulations/exp2_inf/treatment_mean_ALL_2s_anttypesplit.rds")

ggplot(prevalence_mean_treatment_bA, 
       aes(x = time, y = mean_preva,
           color = genotype, #genotype of the seed ant
           fill = genotype,#genotype of the seep ant
           )) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(
      ymin = mean_preva - se_val,
      ymax = mean_preva + se_val
    ),
    alpha = 0.25,
    color = NA
  ) +
  facet_wrap(~ gen_split)+
  scale_color_manual(
    values = anttypes_colors
  ) +
  scale_fill_manual(
    values = anttypes_colors
  ) +
  xlim(0,5000)+
  labs(x = "time [frame]", y = " mean prevalence",
       color = "genotype seedant", fill = "genotype seedant",
       title = " deterministic SI-model on 2s-long interaction network  ") +
  theme_minimal()


################
##look at simulation differnces between actually infected and uninfected individuals,
#for this need to add which ants were exposed in exp2 !!

prevalence_all_inf<-readRDS("simulations/exp2_inf/prevalenceALL_2s.rds")
#mark exposed individuals 
prevalence_all_inf<- prevalence_all_inf|>
  mutate(
    exposed_status = ifelse(
      seed_ant == exposed_ants[colony_name],
      "exposed",
      "nestmate"
    )
  )

#mean over genotype of seed ant in colony for single colonies
prevalence_inf_mean_colony_genotype<-prevalence_all_inf|>
  group_by(colony_name, time, treatment, genotype) |> #add genotype of seedants in simulation to see if that has an effect in groups
  summarise(
    mean_prev = mean(prevalence, na.rm = TRUE),
    sd_val   = sd(prevalence, na.rm = TRUE),
    n             = n(),
    se_val   = sd_val / sqrt(n),
    .groups = "drop"
  )

max_duration_colony<-5000

#plot comparison of different genotype of seed ants pervalence into  directory: /home/sarah/Desktop/EXP1_analysis/simulations/exp2_inf/same_genotype_prevalence/prevalence_plots
for (colony in selected_colonies){
  
  prevalence_plot<-prevalence_inf_mean_colony_genotype|>filter(colony_name == colony)|>ggplot( aes(x= time, y= mean_prev, color= genotype, group = genotype))+
    geom_line(alpha = 0.8, size = 0.8) +
    geom_ribbon(
      aes(ymin = mean_prev - se_val,
          ymax = mean_prev + se_val),
      alpha = 0.25,
      color = NA
    ) +
    scale_color_manual(values = anttypes_colors) +
    labs(x = "time", y = " mean prevalence (over exposed status)",  
         title = paste0("SI-sim: prevalence by seed ant\n colony: ", colony, "\n on min 2s interactions")) +
    xlim(0,max_duration_colony)+
    theme_minimal() +
    theme()
  ggsave(filename = paste0("simulations/exp2_inf/same_genotype_prevalence/prevalence_plots/genotype_mean_percolony/" , 
                           colony, ".png"), plot = prevalence_plot)
  
}


#mean over same treatment colony, 
prevalence_inf_mean_colony<-prevalence_all_inf|>
  group_by(colony_name, time, exposed_status, treatment) |> #add genotype of seedants in simulation to see if that has an effect in groups
  summarise(
    mean_prev = mean(prevalence, na.rm = TRUE),
    sd_val   = sd(prevalence, na.rm = TRUE),
    n             = n(),
    se_val   = sd_val / sqrt(n),
    .groups = "drop"
  )

#mean over same treatment colony, 
prevalence_inf_mean_<-prevalence_inf_mean_colony|>
  group_by(time, exposed_status, treatment) |> #add genotype of seedants in simulation to see if that has an effect in groups
  summarise(
    mean_preva = mean(mean_prev, na.rm = TRUE),
    sd_val   = sd(mean_prev, na.rm = TRUE),
    n             = n(),
    se_val   = sd_val / sqrt(n),
    .groups = "drop"
  )

prevalence_plot_exposed_treatment<-ggplot(prevalence_inf_mean_, aes(x= time, y= mean_preva, color= exposed_status, fill = exposed_status, group = exposed_status))+
  geom_line(alpha = 0.8, size = 0.8) +
  geom_ribbon(
    aes(ymin = mean_preva - se_val,
        ymax = mean_preva + se_val),
    alpha = 0.25,
    color = NA
  ) +
  scale_color_manual(values = exposed_colors) +
  scale_fill_manual(values = exposed_colors) +
  facet_wrap(~treatment)+
  labs(x = "time", y = "mean prevalence over 8 colonies (+- SE)",  
       title = paste0("")) +
  xlim(0,max_duration_colony)+
  theme_minimal() +
  theme()
ggsave(filename = paste0("simulations/exp2_inf/plots/exposed_unexposed_comparison/mean_allcolonies_alltreat.png"), plot = prevalence_plot_exposed_treatment)



max_duration_colony<-5000

#plot comparison of unexposed to exposed pervalence into  directory: exp2_inf/plots/exposed_unexposed_comparison
for (colony in selected_colonies){
  prevalence_plot<-prevalence_inf_mean_colony|>filter(colony_name == colony)|>ggplot( aes(x= time, y= mean_prev, color= exposed_status, group = exposed_status))+
    geom_line(alpha = 0.8, size = 0.8) +
    geom_ribbon(
      aes(ymin = mean_prev - se_val,
          ymax = mean_prev + se_val),
      alpha = 0.25,
      color = NA
    ) +
    scale_color_manual(values = exposed_colors) +
    labs(x = "time", y = " mean prevalence (over exposed status)",  
         title = paste0("SI-sim: prevalence by seed ant\n colony: ", colony, "\n on min 2s interactions")) +
    xlim(0,max_duration_colony)+
    theme_minimal() +
    theme()
  ggsave(filename = paste0("simulations/exp2_inf/plots/exposed_unexposed_comparison/" , 
                           colony, ".png"), plot = prevalence_plot)
  
}

##plot mean over all colonies
#mean over all colonies
prevalence_inf_a_mean_all<-prevalence_inf_a_mean|>
  group_by( time, exposed_status) |> #add genotype of seedants in simulation to see if that has an effect in groups
  summarise(
    mean_preva = mean(mean_prev, na.rm = TRUE),
    sd_val   = sd(mean_prev, na.rm = TRUE),
    n             = n(),
    se_val   = sd_val / sqrt(n),
    .groups = "drop"
  )

prevalence_plot_a<-ggplot(prevalence_inf_a_mean_all, aes(x= time, y= mean_preva, color= exposed_status, fill = exposed_status, group = exposed_status))+
  geom_line(alpha = 0.8, size = 0.8) +
  geom_ribbon(
    aes(ymin = mean_preva - se_val,
        ymax = mean_preva + se_val),
    alpha = 0.25,
    color = NA
  ) +
  scale_color_manual(values = exposed_colors) +
  scale_fill_manual(values = exposed_colors) +
  labs(x = "time", y = "prevalence",  
       title = paste0("a")) +
  xlim(0,max_duration_colony)+
  theme_minimal() +
  theme()
ggsave(filename = paste0("simulations/exp1_inf/plots/exposed_unexposed_comparison/mean_allcolonies_a.png"), plot = prevalence_plot)

#######
#bB
#######
# mixed colonies, spread differences in subgroups bB and bA

prevalence_list_B<-get_colony_files("simulations/exp2_inf/split_genotype_prevalence/B", bB_colonies)
prevalence_list_b<-get_colony_files("simulations/exp2_inf/split_genotype_prevalence/bfrombB", bB_colonies)

prevalence_data_unmerged_B<-lapply(prevalence_list_B, readRDS)
prevalence_data_unmerged_B <- Map(function(df, name) {
  df %>%
    mutate(
      colony_name = name,
      treatment = str_extract(name, "^[^0-9]+"),
      gen_split = "B",
      genotype = unlist(anttype_match[["bB"]][seed_ant])
    )
}, prevalence_data_unmerged_B, names(prevalence_data_unmerged_B))

prevalence_data_unmerged_B <- lapply(prevalence_data_unmerged_B, function(df) {
  df %>%
    mutate(
      treatment = str_extract(colony_name, "^[^0-9]+")
    )
})

prevalence_data_unmerged_b<-lapply(prevalence_list_b, readRDS)
prevalence_data_unmerged_b <- Map(function(df, name) {
  df %>%
    mutate(
      colony_name = name,
      treatment = str_extract(name, "^[^0-9]+"),
      gen_split = "b",
      genotype = unlist(anttype_match[["bB"]][seed_ant])
    )
}, prevalence_data_unmerged_b, names(prevalence_data_unmerged_b))

prevalence_data_unmerged_b <- lapply(prevalence_data_unmerged_b, function(df) {
  df %>%
    mutate(
      treatment = str_extract(colony_name, "^[^0-9]+")
    )
})


#bind all into one df
prevalence_data_merged_B<-bind_rows(prevalence_data_unmerged_B)
prevalence_data_merged_b<-bind_rows(prevalence_data_unmerged_b)

prevalence_data_merged_bB<-bind_rows(prevalence_data_merged_B, prevalence_data_merged_b)
saveRDS( prevalence_data_merged_bB, paste0("simulations/exp2_inf/prevalenceALL_2s_anttype_split_bB.rds"))

prevalence_data_merged_bB<-readRDS("simulations/prevalenceALL_2s_anttype_split.rds")

#plot all in one mean
ggplot(prevalence_data_merged_bB, 
       aes(x = time, y = prevalence, color = treatment)) +
  stat_summary(fun = mean, geom = "line", size = 1.2) +
  scale_color_manual(values = treatment_colors) +
  xlim(0,5000)
theme_minimal()

#############
#mean over colony and treatment
#############
#mean over colony, n = 16 seedants
prevalence_mean_colony<-prevalence_data_merged_bB|>
  group_by(colony_name, time, treatment, genotype, gen_split) |>
  summarise(
    mean_prev = mean(prevalence, na.rm = TRUE),
    sd_val   = sd(prevalence, na.rm = TRUE),
    n             = n(),
    se_val   = sd_val / sqrt(n),
    .groups = "drop"
  )



saveRDS(prevalence_mean_colony, paste0("simulations/exp2_inf/colony_mean_ALL_2s_anttypesplit_bB.rds"))

#mean over treatment, n = 8 replicates per colony
prevalence_mean_treatment<-prevalence_mean_colony|>
  group_by(treatment, time, gen_split, genotype) |>
  summarise(
    mean_preva = mean(mean_prev, na.rm = TRUE),
    n             = n(),
    sd_val   = sd(mean_prev, na.rm = TRUE),
    se_val   = sd_val / sqrt(n),
    .groups = "drop"
  )


saveRDS(prevalence_mean_treatment, paste0("simulations/exp2_inf/treatment_mean_ALL_2s_anttypesplit_bB.rds"))
prevalence_mean_treatment_bA<-readRDS("simulations/exp2_inf/treatment_mean_ALL_2s_anttypesplit_bB.rds")

ggplot(prevalence_mean_treatment, 
       aes(x = time, y = mean_preva,
           color = genotype, #genotype of the seed ant
           fill = genotype,#genotype of the seep ant
       )) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(
      ymin = mean_preva - se_val,
      ymax = mean_preva + se_val
    ),
    alpha = 0.25,
    color = NA
  ) +
  facet_wrap(~ gen_split)+
  scale_color_manual(
    values = anttypes_colors
  ) +
  scale_fill_manual(
    values = anttypes_colors
  ) +
  xlim(0,5000)+
  labs(x = "time [frame]", y = " mean prevalence",
       color = "genotype seedant", fill = "genotype seedant",
       title = " deterministic SI-model on 2s-long interaction network  ") +
  theme_minimal()


################
##look at simulation differnces between actually infected and uninfected individuals,
#for this need to add which ants were exposed in exp2 !!

prevalence_all_inf<-readRDS("simulations/exp2_inf/prevalenceALL_2s.rds")
#mark exposed individuals 
prevalence_all_inf<- prevalence_all_inf|>
  mutate(
    exposed_status = ifelse(
      seed_ant == exposed_ants[colony_name],
      "exposed",
      "nestmate"
    )
  )

#mean over genotype of seed ant in colony for single colonies
prevalence_inf_mean_colony_genotype<-prevalence_all_inf|>
  group_by(colony_name, time, treatment, genotype) |> #add genotype of seedants in simulation to see if that has an effect in groups
  summarise(
    mean_prev = mean(prevalence, na.rm = TRUE),
    sd_val   = sd(prevalence, na.rm = TRUE),
    n             = n(),
    se_val   = sd_val / sqrt(n),
    .groups = "drop"
  )

max_duration_colony<-5000

#plot comparison of different genotype of seed ants pervalence into  directory: /home/sarah/Desktop/EXP1_analysis/simulations/exp2_inf/same_genotype_prevalence/prevalence_plots
for (colony in selected_colonies){
  
  prevalence_plot<-prevalence_inf_mean_colony_genotype|>filter(colony_name == colony)|>ggplot( aes(x= time, y= mean_prev, color= genotype, group = genotype))+
    geom_line(alpha = 0.8, size = 0.8) +
    geom_ribbon(
      aes(ymin = mean_prev - se_val,
          ymax = mean_prev + se_val),
      alpha = 0.25,
      color = NA
    ) +
    scale_color_manual(values = anttypes_colors) +
    labs(x = "time", y = " mean prevalence (over exposed status)",  
         title = paste0("SI-sim: prevalence by seed ant\n colony: ", colony, "\n on min 2s interactions")) +
    xlim(0,max_duration_colony)+
    theme_minimal() +
    theme()
  ggsave(filename = paste0("simulations/exp2_inf/same_genotype_prevalence/prevalence_plots/genotype_mean_percolony/" , 
                           colony, ".png"), plot = prevalence_plot)
  
}


#mean over same treatment colony, 
prevalence_inf_mean_colony<-prevalence_all_inf|>
  group_by(colony_name, time, exposed_status, treatment) |> #add genotype of seedants in simulation to see if that has an effect in groups
  summarise(
    mean_prev = mean(prevalence, na.rm = TRUE),
    sd_val   = sd(prevalence, na.rm = TRUE),
    n             = n(),
    se_val   = sd_val / sqrt(n),
    .groups = "drop"
  )


max_duration_colony<-5000
#plot comparison of unexposed to exposed pervalence into  directory: exp2_inf/plots/exposed_unexposed_comparison
for (colony in selected_colonies){
  prevalence_plot<-prevalence_inf_mean_colony|>filter(colony_name == colony)|>ggplot( aes(x= time, y= mean_prev, color= exposed_status, group = exposed_status))+
    geom_line(alpha = 0.8, size = 0.8) +
    geom_ribbon(
      aes(ymin = mean_prev - se_val,
          ymax = mean_prev + se_val),
      alpha = 0.25,
      color = NA
    ) +
    scale_color_manual(values = exposed_colors) +
    labs(x = "time", y = " mean prevalence (over exposed status)",  
         title = paste0("SI-sim: prevalence by seed ant\n colony: ", colony, "\n on min 2s interactions")) +
    xlim(0,max_duration_colony)+
    theme_minimal() +
    theme()
  ggsave(filename = paste0("simulations/exp2_inf/plots/exposed_unexposed_comparison/" , 
                           colony, ".png"), plot = prevalence_plot)
  
}

##plot mean over all colonies
#mean over all colonies
prevalence_inf_a_mean_all<-prevalence_inf_a_mean|>
  group_by( time, exposed_status) |> #add genotype of seedants in simulation to see if that has an effect in groups
  summarise(
    mean_preva = mean(mean_prev, na.rm = TRUE),
    sd_val   = sd(mean_prev, na.rm = TRUE),
    n             = n(),
    se_val   = sd_val / sqrt(n),
    .groups = "drop"
  )

prevalence_plot_a<-ggplot(prevalence_inf_a_mean_all, aes(x= time, y= mean_preva, color= exposed_status, fill = exposed_status, group = exposed_status))+
  geom_line(alpha = 0.8, size = 0.8) +
  geom_ribbon(
    aes(ymin = mean_preva - se_val,
        ymax = mean_preva + se_val),
    alpha = 0.25,
    color = NA
  ) +
  scale_color_manual(values = exposed_colors) +
  scale_fill_manual(values = exposed_colors) +
  labs(x = "time", y = "prevalence",  
       title = paste0("a")) +
  xlim(0,max_duration_colony)+
  theme_minimal() +
  theme()
ggsave(filename = paste0("simulations/exp1_inf/plots/exposed_unexposed_comparison/mean_allcolonies_a.png"), plot = prevalence_plot)



