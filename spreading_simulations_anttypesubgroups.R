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

spread_SI <-function(edge_list, spreading_prob, seed_ant, duration, present_ants){
  #SI-spreading process on timeprdered edgelist with a seed ant
  
  
  
  
  infection_hist<-matrix(0,ncol=length(present_ants), nrow=duration)
  a_hist<-matrix(0,ncol=length(present_ants), nrow=duration)
  colnames(a_hist)<-present_ants
  b_hist<-matrix(0,ncol=length(present_ants), nrow=duration)
  colnames(b_hist)<-present_ants
  colnames(infection_hist)<-present_ants
  #SI, no recovery infect seed ant
  
  infection_hist[1:duration, seed_ant]<-1
  
  
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
        if(genotype_match[[partner2]] == "a"){
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
    a_prevalence[, ant] <- result$prevalence
  }
  #convert to dataframe
  prevalence_df<-data.frame(time = 1:length(result$prevalence), seedvariation_prevalence)
  
  prevalence_long <- prevalence_df %>%
    pivot_longer(cols = -time, 
                 names_to = "seed_ant", 
                 values_to = "prevalence")
  
  prevalence_long<-prevalence_long%>%mutate(genotype =case_when(
    treatment_match[colony_name] == "3" ~ unlist(genotype_match[seed_ant]),
    treatment_match[colony_name] == "1" ~ "a",                            # all ants in pure colony
    treatment_match[colony_name] == "2" ~ "b"     
  ))
  
  saveRDS( prevalence_long, paste0("simulations/",colony_name,"_prevalencesim_3s.rds"))
  
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
  
  ggsave(filename = paste0("simulations/prevalence_3s_cor_" , 
                           colony_name , ".png"), plot = prevalence_plot)
  
  
}
