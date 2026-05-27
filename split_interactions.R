#number of interactions in mixed colonies 
#EXP2 mixed interactions from bB and bA
#edgelist
################
setwd("Desktop/EXP1_analysis")
source("function_collection.R")
source("global_parameters_exp2_inf.R")

source_folder_edgelist<-"edge_lists/exp2_12hs/"
folder_path_networkpara <-"network_parameters/exp2_inf/5mins/"
output_folder <- "network_parameter_plots/exp2_inf/5mins/time_mean/interactions_seperated"
timestamp <- format(Sys.time(), "%d%m%Y")

#######
# for function all interavalls
interval_time<-5
number_intervalls<-15
time_between<-121
#upper and lower interaction limit
interaction_limit_min<-2
interaction_limit_sec<-1
#get network parameters 
network_obj_list<-get_colony_files(source_folder_edgelist, bB_colonies)
for (i in 1:length(bB_colonies)){
  
  
  colony_name<-selected_colonies[[i]]
  colony_name<-bB_colonies[[1]]
  present_ants<-present_ants_list[[1]]
  cat("number interactions", colony_name, "\n")
  framerate_col<-framerate[[colony_name]]
  lower_limit_interaction<-interaction_limit_sec*framerate_col
  higher_limit_interaction<-interaction_limit_min*60*framerate_col
  
  #intervall_list<-[]
  ###################
  
  #prepare network obj, interaction limits
  network_obj_raw <-readRDS(network_obj_list[[colony_name]])
  
 
  missing_ants <- setdiff(expected_ants, present_ants)
  
  network_obj <- network_obj_raw[
    !network_obj_raw$head %in% missing_ants &
      !network_obj_raw$tail %in% missing_ants,
  ]
  
  #pick time slice of network object from 12hs
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
  
  edge_list$genotype_head<-
  edge_list$genotype_head<-tail
  edge_list$interactiontype<-paste0()
  
  # ADD genotype of a genotype_head, genotype_tail, interactiontype = bb, BB, bB
  
  
  
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
