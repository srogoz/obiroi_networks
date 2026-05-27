##get prevalence for seperate genotypes
##################
#PARAMETERS
##################
#functional
interaction_limit_min<-2
interaction_limit_sec<-3
interval_time<-90 # in minutes
duration<- 50000 #simulation-duration
source_folder_edgelist<-"edge_lists/exp1_12hs/"
network_obj_list<-get_colony_files(source_folder_edgelist, selected_colonies)

#########
#organisational

folder_path_networkpara <-"simulations/"
output_folder <- "simulations"
timestamp <- format(Sys.time(), "%d%m%Y")
############SI-simulations-function

spread_SI_anttype_subsets <-function(edge_list, spreading_prob, seed_ant, duration, present_ants){
  #SI-spreading process on timeprdered edgelist with a seed ant
  
  
  
  
  infection_hist<-matrix(0,ncol=length(present_ants), nrow=duration)
  a_hist<-matrix(0,ncol=length(present_ants), nrow=duration)
  colnames(a_hist)<-present_ants
  b_hist<-matrix(0,ncol=length(present_ants), nrow=duration)
  colnames(b_hist)<-present_ants
  colnames(infection_hist)<-present_ants
  #SI, no recovery infect seed ant
  
  infection_hist[1:duration, seed_ant]<-1
  
  #assign genotye of seedant
  seed_idx <- match(seed_ant, present_ants)
  gen <- genotype_match[[seed_ant]]
  
  if (gen == "a") {
    a_hist[1:duration, seed_idx] <- 1
  } else if (gen == "b") {
    b_hist[1:duration, seed_idx] <- 1
  }
  
  #order interactions by time
  timeordered <- edge_list[order(edge_list$onset),]
  number_interactions<-nrow(timeordered)
  
  
  for (i in seq_len(number_interactions)){
    time_step = timeordered$onset[i]
    
    #respect simulation duration
    if (time_step > duration) next
    
    #interaction-partners
    partner1<-timeordered$tail[i]
    partner2<-timeordered$head[i]
    
    #collect infected ants
    infected_ants<-present_ants[infection_hist[time_step, ] == 1]
    
    
    #if both infected, skip
    
    if (partner2 %in% infected_ants && partner1 %in% infected_ants){
      next
    }
    
    # transmission from 1 to 2 with prob
    if (partner1 %in% infected_ants && !(partner2 %in% infected_ants)){
      if (runif(1) < spreading_prob) {
        # once infected, stays infected for all future times
        seed_idx <- match(partner2, present_ants)
        
        infection_hist[time_step:duration, seed_idx] <-1
        if(genotype_match[[partner2]] == "a"){
          a_hist[time_step:duration, seed_idx] <-1
            }else{
              b_hist[time_step:duration, seed_idx] <-1
            }
      }
    }
    #transmission from 2 to 1
    if (partner2 %in% infected_ants && !(partner1 %in% infected_ants)){
      if(runif(1) < spreading_prob){
        seed_idx <- match(partner1, present_ants)
        infection_hist[time_step:duration, seed_idx ] <-1
        if(genotype_match[[partner1]] == "a"){
          a_hist[time_step:duration, seed_idx] <-1
        }else{
          b_hist[time_step:duration, seed_idx] <-1
        }
      }
    }
    #end simulation if everyone is infected  
    if (length(infected_ants) == length(present_ants)) break
    saturated<-timeordered$terminus[i]
    
  }
  
  prevalence<-rowSums(infection_hist)
  prevalence_a<-rowSums(a_hist)
  prevalence_b<-rowSums(b_hist)
  return(list(infection_hist = infection_hist, 
              prevalence_a= prevalence_a,
              prevalence_b= prevalence_b,
              prevalence = prevalence,
              saturated = saturated)
  )
  
}

###########
#loop over all colonies
for (i in 1:length(mixed_colonies)){
  
  
  colony_name<-mixed_colonies[[i]]
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
  #collect all prevalences over all seedants
  seedvariation_prevalence <- matrix(0, nrow = duration, ncol = length(present_ants))
  seedvariation_saturationtime<-matrix(0, nrow = 1, ncol = length(present_ants))
  colnames(seedvariation_saturationtime) <- present_ants
  colnames(seedvariation_prevalence) <- present_ants
  #split prevalences by genotype, to check id genotype clusters saturate faster within their own genotype
  seedvariation_prevalence_a<-matrix(0, nrow = duration, ncol = length(present_ants))
  colnames(seedvariation_prevalence_a) <- present_ants
  seedvariation_prevalence_b<-matrix(0, nrow = duration, ncol = length(present_ants))
  colnames(seedvariation_prevalence_b) <- present_ants
  
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
    seedvariation_prevalence_a[, ant] <- result$prevalence_a
    seedvariation_prevalence_b[, ant] <- result$prevalence_b
  }
  #convert to dataframe
  prevalence_df<-data.frame(time = 1:length(result$prevalence), seedvariation_prevalence)
  prevalence_a_df<-data.frame(time = 1:length(result$prevalence), seedvariation_prevalence_a)
  prevalence_b_df<-data.frame(time = 1:length(result$prevalence), seedvariation_prevalence_b)
  
  prevalence_a_long <- prevalence_a_df %>%
    pivot_longer(cols = -time, 
                 names_to = "seed_ant", 
                 values_to = "prevalence")
  
  prevalence_b_long<- prevalence_b_df %>%
    pivot_longer(cols = -time, 
                 names_to = "seed_ant", 
                 values_to = "prevalence")
  
  prevalence_long <- prevalence_df %>%
    pivot_longer(cols = -time, 
                 names_to = "seed_ant", 
                 values_to = "prevalence")
  
  prevalence_long<-prevalence_long%>%mutate(genotype =case_when(
    treatment_match[colony_name] == "3" ~ unlist(genotype_match[seed_ant])
  ))
  
  prevalence_a_long<-prevalence_a_long%>%mutate(genotype =case_when(
    treatment_match[colony_name] == "3" ~ unlist(genotype_match[seed_ant])
  ))
  prevalence_b_long<-prevalence_b_long%>%mutate(genotype =case_when(
    treatment_match[colony_name] == "3" ~ unlist(genotype_match[seed_ant])
  ))
  
  saveRDS( prevalence_long, paste0("simulations/exp1_12hs/split_genotype_prevalence/",colony_name,"_prevalencesim_2s.rds"))
  saveRDS( prevalence_a_long, paste0("simulations/exp1_12hs/split_genotype_prevalence/",colony_name,"_prevalencesim_2s_a.rds"))
  saveRDS( prevalence_b_long, paste0("simulations/exp1_12hs/split_genotype_prevalence/",colony_name,"_prevalencesim_2s_b.rds"))
  
  max_duration_colony<-max(seedvariation_saturationtime)
  
  #all ants in one plot per colony
  prevalence_a_plot<-ggplot(prevalence_a_long, aes(x= time, y= prevalence, color= genotype, group = seed_ant))+
    geom_line(alpha = 0.8, size = 0.8) +
    scale_color_manual(values = anttypes_colors) +
    labs(x = "time", y = "prevalence genotype a ",  
         title = paste0("SI-simulation: prevalence by seed ant\n colony on genotype a : ", colony_name, "\n on min 2s interactions")) +
    xlim(0,max_duration_colony)+
    theme_minimal() +
    theme()
  
  ggsave(filename = paste0("simulations/exp1_12hs/split_genotype_prevalence/" , 
                           colony_name , "_a.png"), plot = prevalence_a_plot)
  
  prevalence_b_plot<-ggplot(prevalence_b_long, aes(x= time, y= prevalence, color= genotype, group = seed_ant))+
    geom_line(alpha = 0.8, size = 0.8) +
    scale_color_manual(values = anttypes_colors) +
    labs(x = "time", y = "prevalence genotype a ",  
         title = paste0("SI-simulation: prevalence by seed ant\n colony on genotype b : ", colony_name, "\n on min 2s interactions")) +
    xlim(0,max_duration_colony)+
    theme_minimal() +
    theme()
  
  ggsave(filename = paste0("simulations/exp1_12hs/split_genotype_prevalence/" , 
                           colony_name , "_b.png"), plot = prevalence_b_plot)
  
  
  
  prevalence_plot<-ggplot(prevalence_long, aes(x= time, y= prevalence, color= genotype, group = seed_ant))+
    geom_line(alpha = 0.8, size = 0.8) +
    scale_color_manual(values = anttypes_colors) +
    labs(x = "time", y = "prevalence",  
         title = paste0("SI-simulation: prevalence by seed ant\n colony: ", colony_name, "\n on min 3s interactions")) +
    xlim(0,max_duration_colony)+
    theme_minimal() +
    theme()
  
  ggsave(filename = paste0("simulations/exp1_12hs/split_genotype_prevalence/" , 
                           colony_name , ".png"), plot = prevalence_plot)
  
  
}

#check speed of transmission on subset of genotype/ anttype groups

#overlay prevalences of all colonies in one treatment in one plot
prevalence_a_list<-get_colony_files("simulations/exp1_12hs/split_genotype_prevalence/prevalence_a/", mixed_colonies)
prevalence_b_list<-get_colony_files("simulations/exp1_12hs/split_genotype_prevalence/prevalence_b/", mixed_colonies)

prevalence_data_a_unmerged<-lapply(prevalence_a_list, readRDS)
prevalence_data_a_unmerged <- Map(function(df, name) {
  df %>%
    mutate(
      colony_name = name,
      treatment = str_extract(name, "^[^0-9]+")
    )
}, prevalence_data_a_unmerged, names(prevalence_data_a_unmerged))

prevalence_data_a_unmerged <- lapply(prevalence_data_a_unmerged, function(df) {
  df %>%
    mutate(
      treatment = str_extract(colony_name, "^[^0-9]+")
    )
})


#bind all into one df
prevalence_data_a_merged<-bind_rows(prevalence_data_a_unmerged)
prevalence_data_a_merged<-prevalence_data_a_merged%>%mutate(gen_split = "a")

#check-in for b
prevalence_data_b_unmerged<-lapply(prevalence_b_list, readRDS)
prevalence_data_b_unmerged <- Map(function(df, name) {
  df %>%
    mutate(
      colony_name = name,
      treatment = str_extract(name, "^[^0-9]+")
    )
}, prevalence_data_b_unmerged, names(prevalence_data_b_unmerged))

prevalence_data_b_unmerged <- lapply(prevalence_data_b_unmerged, function(df) {
  df %>%
    mutate(
      treatment = str_extract(colony_name, "^[^0-9]+")
    )
})


#bind all into one df
prevalence_data_b_merged<-bind_rows(prevalence_data_b_unmerged)
prevalence_data_b_merged<-prevalence_data_b_merged%>%mutate(gen_split = "b")

#merge genotype ab and a subgroup spread
prevalence_data_merged<-bind_rows(prevalence_data_a_merged, prevalence_data_b_merged)
saveRDS( prevalence_data_merged, paste0("simulations/exp1_12hs/split_genotype_prevalence/prevalence_gensplitALL_2s.rds"))

prevalence_data_merged<-readRDS("simulations/prevalenceALL_2s.rds")

#plot all in one mean
ggplot(prevalence_data_merged, 
       aes(x = time, y = prevalence, color = treatment)) +
  stat_summary(fun = mean, geom = "line", size = 1.2) +
  scale_color_manual(values = treatment_colors) +
  xlim(0,5000)
theme_minimal()

#############
#mean over colony and treatment for gen_split
#############
#mean over colony, n = 16 seedants
prevalence_mean_colony<-prevalence_data_merged|>
  group_by(colony_name, time, gen_split, genotype) |> #add genotype of seedants in simulation to see if that has an effect in groups
  summarise(
    mean_prev = mean(prevalence, na.rm = TRUE),
    sd_val   = sd(prevalence, na.rm = TRUE),
    n             = n(),
    se_val   = sd_val / sqrt(n),
    .groups = "drop"
  )



saveRDS(prevalence_mean_colony, paste0("simulations/colony_mean_ALL_2s.rds"))

#mean over treatment, n = 8 replicates per colony
prevalence_mean_gensplit<-prevalence_mean_colony|>
  group_by(gen_split, genotype, time) |> #add genotype of seed - ant to see if that takes an effect in groups
  summarise(
    mean_preva = mean(mean_prev, na.rm = TRUE),
    n             = n(),
    sd_val   = sd(mean_prev, na.rm = TRUE),
    se_val   = sd_val / sqrt(n),
    .groups = "drop"
  )


saveRDS(prevalence_mean_treatment, paste0("simulations/treatment_mean_ALL_2s.rds"))

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
  facet_wrap(~ gen_split) +  #subgraphs in different pannels for different genotype of seed_ant
  scale_color_manual(
    values = treatment_colors
  ) +
  scale_fill_manual(
    values = treatment_colors
  ) +
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

##for mixed colonies check seed ant effect
prevalence_a_seperated<-prevalence_data_merged%>%filter(treatment == "ba")

prevalence_a_seperated<-prevalence_a_seperated|>group_by(colony_name, genotype, time) |>
  summarise(
    mean_prev = mean(prevalence, na.rm = TRUE),
    sd_val   = sd(prevalence, na.rm = TRUE),
    n             = n(),
    se_val   = sd_val / sqrt(n),
    .groups = "drop"
  )

prevalence_a_seperated_mean<-prevalence_a_seperated|>group_by(genotype, time) |>
  summarise(
    mean_preva = mean(mean_prev, na.rm = TRUE),
    sd_val   = sd(mean_prev, na.rm = TRUE),
    n             = n(),
    se_val   = sd_val / sqrt(n),
    .groups = "drop"
  )

##plot mixed colonies mean over all colonies all a and seperated all b seed ants
ggplot(prevalence_a_seperated_mean, 
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
  scale_color_manual(
    values = treatment_colors
  ) +
  scale_fill_manual(
    values = treatment_colors
  ) +
  xlim(0,5000)+
  labs(x = "time [frame]", y = " mean prevalence",
       color = "genotype", fill = "genotype",
       title = " deterministic SI-simulation on 2s-long interaction network \n mean over mixed colonies sepereated by genotype of seed ant  ") +
  theme_minimal()

#add complete treatment ba mixed! 

#colony_level
for (colony in mixed_colonies){
  mixed_prevalence_plot<-prevalence_mean_colony|>filter(colony_name == colony)|>
    ggplot( aes(x = time, y = mean_prev, color = gen_split,
                fill  = gen_split)) +
    geom_line(linewidth = 1) +
    geom_ribbon(
      aes(
        ymin = mean_prev - se_val,
        ymax = mean_prev + se_val
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
    labs(x = "time [frame]", y = " mean prevalence for genotype subgroups",
         color = "treatment", fill = "treatment",
         title = " deterministic SI-model on 2s-long interaction network  ") +
    theme_minimal()
  
  ggsave(filename = paste0("simulations/exp1_12hs/split_genotype_prevalence/plots/" , 
                           colony , ".png"), plot = mixed_prevalence_plot)
  
}

prevalence_plot<-ggplot(prevalence_data_merged, aes(x= time, y= prevalence, color= treatment))+
  geom_line(alpha = 0.6, size = 0.8) +
  scale_color_manual(values = treatment_colors) +
  labs(x = "time", y = "prevalence",  
       title = paste0("SI- model: prevalence by seed Ant\n colony: ", colony_name)) +
  theme_minimal() +
  theme()

ggsave(filename = paste0("simulations/prevalence_" , 
                         colony_name , ".png"), plot = mixed_prevalence_plot)


###plot infected and uninfected dynamics together per treatment
#infected:
prevalence_mean_inf<-readRDS("simulations/treatment_mean_ALLINF.rds")
prevalence_mean_inf<-prevalence_mean_treatment
prevalence_mean_inf$infection<-"after"
prevalence_mean_uninf<-readRDS("simulations/exp1_12hs/data/treatment_mean_ALL.rds")
prevalence_mean_uninf$infection<-"before"

#
prevalence_mean_total<-bind_rows(prevalence_mean_inf, prevalence_mean_uninf)

a<-prevalence_mean_total%>%filter(treatment == "a")
b<-prevalence_mean_total%>%filter(treatment == "b")
ba<-prevalence_mean_total%>%filter(treatment == "ba")

plot_a<-ggplot(a, 
               aes(x = time, y = mean_prev, color = infection,
                   fill  = treatment)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(
      ymin = mean_prev - sd_val,
      ymax = mean_prev + sd_val
    ),
    alpha = 0.25,
    color = NA
  ) +
  scale_color_manual(
    values = infection_status
  ) +
  scale_fill_manual(
    values = infection_status
  ) +
  xlim(0,5000)+
  labs(x = "time [s]", y = " mean prevalence",
       color = "treatment", fill = "treatment",
       title = " deterministic SI-model for treatment a\n before and after fungal exposure ") +
  
  theme_minimal()

plot_b<-ggplot(b, 
               aes(x = time, y = mean_prev, color = infection,
                   fill  = treatment)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(
      ymin = mean_prev - sd_val,
      ymax = mean_prev + sd_val
    ),
    alpha = 0.25,
    color = NA
  ) +
  scale_color_manual(
    values = infection_status
  ) +
  scale_fill_manual(
    values = infection_status
  ) +
  xlim(0,5000)+
  labs(x = "time [s]", y = " mean prevalence",
       color = "treatment", fill = "treatment",
       title = " deterministic SI-model for treatment b\n before and after fungal exposure ") +
  
  theme_minimal()

plot_ba<-ggplot(ba, 
                aes(x = time, y = mean_prev, color = infection,
                    fill  = treatment)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(
      ymin = mean_prev - sd_val,
      ymax = mean_prev + sd_val
    ),
    alpha = 0.25,
    color = NA
  ) +
  scale_color_manual(
    values = infection_status
  ) +
  scale_fill_manual(
    values = infection_status
  ) +
  xlim(0,5000)+
  labs(x = "time [s]", y = " mean prevalence",
       color = "treatment", fill = "treatment",
       title = " deterministic SI-model for treatment ba\n before and after fungal exposure ") +
  
  theme_minimal()

library(patchwork)
grid_plot<-(plot_a| plot_b| plot_ba) +
  plot_layout( axes = "collect" , guides = "collect") 

ggsave(paste0("simulations/exp1_12hs/plots/mean_before_after.png"), 
       plot =grid_plot,
       width = 16,
       height = 5,
       dpi = 300)



