library(sna)
library(tsna)
#library(ndtv)
library(reticulate)
library(tidyr)
library(igraph)
library(intergraph)
library(networkLite)
library(network)
library(networkDynamic)
library(dplyr)
library(ggplot2)
library(viridisLite)
library(viridis)
library(glmmTMB)
library(rstatix)
library(emmeans)
library(ggpubr)

# Load .npy file function
np <- import("numpy")

################
setwd("Desktop/EXP1_analysis")
source("function_collection.R")
source("global_parameters_exp1_12hs.R")
source_folder_edgelist<-"edge_lists/exp1_12hs/"
folder_path_networkpara <-"network_parameters/exp1_12hs/5_mins/120_135/"
output_folder <- "network_parameter_plots/exp1_12hs/5_mins/"
timestamp <- format(Sys.time(), "%d%m%Y")

#set truncation time for infected 75mins
interval_time<-5
number_intervalls<-15
#unpper and lower interaction limit
interaction_limit_min<-2
interaction_limit_sec<-1
#get network parameters
network_obj_list<-get_colony_files(source_folder_edgelist, selected_colonies)

####get network parameters for timepoints of interest ###################################
#antlevel parameter collection
strength_antlevel_mean<-data.frame(row.names = expected_ants)
strength_antlevel_sd<-data.frame(row.names = expected_ants)
local_efficiency<-data.frame(row.names = expected_ants)
waiting_time_ant_mean_df<-data.frame(row.names = expected_ants)
waiting_time_ant_sd_df<-data.frame(row.names = expected_ants)
interactionlength_ant_mean<-data.frame(row.names = expected_ants)
burstiness_ant_df<-data.frame(row.names = expected_ants)
clustering_local_df<-data.frame(row.names = expected_ants)
centrality_antlevel_df<-data.frame(row.names = expected_ants)

results_networkparameters_allintervalls<-data.frame(colony_name = NA,
                                                    time_interval= NA ,
                                                    number_ants=NA,
                                                    total_number_interactions=NA,
                                                    total_sum_of_interactions= NA,
                                                    strength_mean = NA,
                                                    strength_sd= NA,
                                                    interaction_length_mean= NA,
                                                    interaction_length_sd= NA,
                                                    waiting_time_mean=NA,
                                                    waiting_time_sd=NA,
                                                    burstiness=NA,
                                                    mean_distance = NA,
                                                    density_aggr=NA,
                                                    global_eff=NA,
                                                    assortativity=NA,
                                                    clustering_global= NA,
                                                    started_but_not_ended=NA)

for (i in 1:length(selected_colonies)){
  
  colony_name<-selected_colonies[[i]]
  framerate_col<-framerate[[colony_name]]
  present_ants<-present_ants_list[[colony_name]]
  duration<-interval_time*60*framerate_col
  lower_limit_interaction<-interaction_limit_sec*framerate_col
  higher_limit_interaction<-interaction_limit_min*60*framerate_col
  
  ###################
  
  #prepare network obj, interaction limits
  network_obj <-readRDS(network_obj_list[[colony_name]])
  
 ##adjust intervallsize for other timelines
  
  interval_size<-interval_time*60*framerate_col
  intervall_timepoints <- seq(
    from = 0,
    by   = interval_size,
    #number intervalls +1 (+time between)
    length.out = number_intervalls + 121
  )
  
  for (j in 1:number_intervalls){
    #start and stop
    #j, add same length if later timeintervalls
    time_window_start<-intervall_timepoints[j+120]
    #j+1 add same length if later timeintervalls
    time_window_end<-intervall_timepoints[j+121]
    
    #limit to timeintervall
    network_obj_5<- network_obj[
      network_obj$onset >= time_window_start &
        network_obj$terminus  <= time_window_end,]
    
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
    network_obj_5 <-network_obj_5[(network_obj_5$duration>lower_limit_interaction),]
    network_obj_5$duration <- pmin(network_obj_5$duration, higher_limit_interaction)
    
    
    #remove (empty) data of missing ants
    missing_ants <- setdiff(expected_ants, present_ants)
    
    network_obj_5 <- network_obj_5[
      !network_obj_5$head %in% missing_ants &
        !network_obj_5$tail %in% missing_ants,
    ]
    
    #force empty entries of ants who are in theory present, but dont have any interactions in this time period
    
    #count how many interaction and percentages of frames ommited
    interaction_control <- setNames(
      data.frame(matrix(NA, nrow = length(selected_colonies), ncol = 6)),
      c("sum_interactions", "sum_frames", "frames_upper", 
        "frames_lower", "interactions_upper", "interactions_lower")
    )
    
    rownames(interaction_control) <- selected_colonies
    
    sum_interactionframes<-sum(network_obj$duration)
    number_interactions<-length(network_obj$duration)
    lower_cut_framesum<-sum(network_obj[(network_obj$duration<=lower_limit_interaction),]$duration)
    lower_cut_interactionumber<-length(network_obj[(network_obj$duration<=lower_limit_interaction),]$duration)
    lower_cut_framepercentage<-round(lower_cut_framesum/ (sum_interactionframes/100),2)
    lower_cut_interactionumberpercentage<-round(lower_cut_interactionumber/ (number_interactions/100),2)
    
    
    upper_cut_framesum<-sum(network_obj[(network_obj$duration>=higher_limit_interaction),]$duration)
    upper_cut_interactionumber<-length(network_obj[(network_obj$duration>=higher_limit_interaction),]$duration)
    upper_cut_framepercentage<-round(upper_cut_framesum/ (sum_interactionframes/100),2)
    upper_cut_interactionumberpercentage<-round(upper_cut_interactionumber/ (number_interactions/100),2)
    
    sum_discarded_frames<-upper_cut_framepercentage+lower_cut_framepercentage
    affected_interactions<-lower_cut_interactionumberpercentage+upper_cut_interactionumberpercentage
    #########
    interaction_control[colony_name,"sum_interactions"]<-affected_interactions
    interaction_control[colony_name, "sum_frames"]<-sum_discarded_frames
    interaction_control[colony_name, "frames_upper"]<-upper_cut_framepercentage
    interaction_control[colony_name, "frames_lower"]<-lower_cut_framepercentage
    interaction_control[colony_name, "interactions_upper"]<-upper_cut_interactionumberpercentage
    interaction_control[colony_name, "interactions_lower"]<-lower_cut_interactionumberpercentage
    
    
    ##
    #plot histogram of interaction durations
    
    # p_hist <- ggplot(network_obj, aes(duration)) +
    #   geom_histogram(
    #     aes(fill = after_stat(
    #       ifelse(
    #         x < lower_limit_interaction, "below",
    #         ifelse(x > higher_limit_interaction, "above", "inside")
    #       )
    #     )),
    #     bins = 30,
    #     color = "grey"
    #   )  +
    #   scale_fill_manual(
    #     values = c(
    #       "below" = "grey",
    #       "inside" = "magenta",
    #       "above" = "seagreen2"
    #     ),
    #     labels = c(
    #       below  = paste0("< 1s (removed)\n(", lower_cut_framepercentage ," /",lower_cut_interactionumberpercentage, " )"),
    #       inside = paste0("unaffected"),
    #       above  = paste0("< 2 min (truncated)\n(", upper_cut_framepercentage, " /" ,upper_cut_interactionumberpercentage, " )")
    #     ),
    #     name = "Thresholds \n(% of frames, \ninteractions affected)"
    #   ) +
    #   scale_y_log10() +   #log-scale on y-axis
    #   geom_text(
    #     stat = "bin",
    #     bins = 40,
    #     aes(label = after_stat(count)),
    #     vjust = -0.5,
    #     size = 3
    #   ) +
    #   theme_minimal() +
    #   labs(fill = paste0("Bin midpoint > ", lower_limit_interaction),
    #        y = "Count (log scale)")+
    #   ggtitle(paste0("Interactionlength histogram for ", colony_name))+
    #   annotate(
    #     "text",
    #     x = Inf, y = Inf,
    #     hjust = 1.1, vjust = 2,
    #     label = paste0(" total frames :" ,sum_discarded_frames ," %", "\ntotal interactions: " ,affected_interactions ," %" ),
    #     size = 4
    #   ) +
    #   theme(
    #     # smaller labels
    #     legend.key.height = unit(2.2, "lines"),        # more vertical spacing
    #     # more horizontal spacing
    #   )
    # 
    # ggsave(paste0("network_parameter_plots/antlevel_1_truncated_2min/interaction_histograms/hist_interactions_",colony_name,".png"), 
    #        plot = p_hist,
    #        width = 6,
    #        height = 5,
    #        dpi = 300)
    # 
    # 
    #########
    #present_ants<-unique(c(network_obj_5$tail,network_obj_5$head))
    #problem with present ants being not properly assigned
    #present_ants<-unique(c(network_obj_5$tail,network_obj_5$head))
    ######################
    #on ant level interactiontime, burstiness, waitingtime_mean, waitingtime_sd
    
    interactionlength_means_vec <- rep(NA_real_, length(expected_ants))
    names(interactionlength_means_vec) <- expected_ants
    waiting_time_ant_mean_vec <- rep(NA_real_, length(expected_ants))
    names(waiting_time_ant_mean_vec) <- expected_ants
    waiting_time_ant_sd_vec <- rep(NA_real_, length(expected_ants))
    names(waiting_time_ant_sd_vec) <- expected_ants
    burstiness_ant_vec <- rep(NA_real_, length(expected_ants))
    names(burstiness_ant_vec) <- expected_ants
    
    for (ant_name in present_ants) {
      df_ant <- network_obj_5[
        network_obj_5$head == ant_name | network_obj_5$tail == ant_name, 
      ]
      
      interactionlength_means_vec[ant_name] <- mean(df_ant$duration)
      
      wt <- diff(sort(df_ant$onset))
      wt <- wt[wt > 0] / framerate_col
      
      waiting_time_ant_mean_vec[ant_name] <- mean(wt)
      waiting_time_ant_sd_vec[ant_name] <- sd(wt)
      
      burstiness_ant_vec[ant_name] <- (waiting_time_ant_sd_vec[ant_name] - waiting_time_ant_mean_vec[ant_name]) /
        (waiting_time_ant_sd_vec[ant_name] + waiting_time_ant_mean_vec[ant_name])
    }
    
    #assign to global file
    waiting_time_ant_mean_df[[paste0(colony_name,"_",j)]]<-interactionlength_means_vec
    waiting_time_ant_sd_df[[paste0(colony_name,"_",j)]]<-waiting_time_ant_sd_vec
    interactionlength_ant_mean[[paste0(colony_name,"_",j)]]<-interactionlength_means_vec
    burstiness_ant_df[[paste0(colony_name,"_",j)]]<-burstiness_ant_vec
    
    #initialize columnnames with numbers for networks
    network_obj_5$tail<-match(network_obj_5$tail, expected_ants)
    network_obj_5$head<-match(network_obj_5$head, expected_ants)
    
    #cretae aggregated network with matching interaction-time filters
    prelim_agg_network<-aggregate_from_edgelist(network_obj_5)
    
    #force empty entries of ants who are in theory present, but dont have any interactions in this time period to have zero values! 
    present_ants_numbers<-match(present_ants,expected_ants)
    aggregated_network <- matrix(
      0,
      nrow = length(present_ants_numbers),
      ncol = length(present_ants_numbers),
      dimnames = list(present_ants_numbers, present_ants_numbers)
    )
    aggregated_network[rownames(prelim_agg_network), colnames(prelim_agg_network)] <- prelim_agg_network
    
    
    
    #matrix_dim<- max (dim(aggregated_network))
    #aggregated_network[,matrix_dim]<-t(aggregated_network[matrix_dim,])
    #aggregated_network[matrix_dim,]<-0
    
    #flexible dependency of matrix dimension
    #aggregated_network[,16]<-t(aggregated_network[16,])
    #aggregated_network[16,]<-0
    #aggregated network of interactions with hs as unit corrected for 
    network_aggregated_symm<-aggregated_network
    network_aggregated_symm[lower.tri(network_aggregated_symm)]<-t(aggregated_network)[lower.tri(aggregated_network)]
    #network_aggregated symm transformed into hours as unit, correcting for the framerate difference
    network_aggregated_symm<-round(network_aggregated_symm/((framerate_col)),2) 
    #into igraph object
    graph_aggregated<-graph_from_adjacency_matrix(network_aggregated_symm, weighted = TRUE, mode = "undirected", add.colnames = NULL)
    
    #net_density <- get_density(network_obj, present_ants)
    #density_temporal_mean<-mean(net_density)
    #density_temporal_sd<-sd(net_density)
    
    ##################################################
    #AGGREGATED MEASURES
    number_ants<-length(present_ants)
    
    total_number_interactions<-length(network_obj_5$duration)
    #results_networkparameters$total_number_interactions<-total_number_interactions
    total_sum_of_interactions<-(sum(network_obj_5$duration))/framerate_col
    #frame change, 10 frames and 5 frames
    density_aggr<-total_sum_of_interactions/(((number_ants*(number_ants-1))/2)*duration)
    
    #colony level, from aggregated network
    strength_mean<-mean(strength(graph_aggregated))
    strength_sd<-sd(strength(graph_aggregated))
    
    #antlevel measure
    strength_collection<-c(strength(graph_aggregated))
    strength_collection <- setNames(
      c(strength(graph_aggregated)),
      present_ants
    )
    
    strength_full <- setNames(
      rep(NA_real_, length(expected_ants)),
      expected_ants
    )
    # fill in values for ants that exist
    strength_full[names(strength_collection)] <- strength_collection
    
    #assign to collection
    strength_antlevel_mean[[paste0(colony_name,"_",j)]]<-strength_full
    
    V(graph_aggregated)$strength <- strength(graph_aggregated, weights = E(graph_aggregated)$weight)
    ass_ants<-assortativity(graph_aggregated, types1   = V(graph_aggregated)$strength,
                            directed = FALSE
    )
    
    
    #clustering coefficient global and local
    
    clustering_global<-transitivity(
      graph_aggregated,
      type ="global",
      vids = NULL,
      weights = E(graph_aggregated)$weight,
      isolates = c("NaN")
    )
    
    clustering_local<-transitivity(
      graph_aggregated,
      type ="barrat",
      vids = NULL,
      weights = E(graph_aggregated)$weight,
      isolates = c("NaN")
    )
    
    #antlevel measure
    clustering_collection<-c(clustering_local)
    clustering_collection <- setNames(
      c(clustering_local),
      present_ants
    )
    
    clustering_full <- setNames(
      rep(NA_real_, length(expected_ants)),
      expected_ants
    )
    # fill in values for ants that exist
    clustering_full[names(clustering_collection)] <- clustering_collection
    
    #assign to collection
    clustering_local_df[[paste0(colony_name,"_",j)]]<-clustering_full
    
    ####Local and global Efficiency
    #colony level:
    global_eff<-global_efficiency(graph_aggregated, weights = E(graph_aggregated)$weight)
    
    #ant level
    local_eff<-local_efficiency(
      graph_aggregated,
      vids = V(graph_aggregated),
      weights = E(graph_aggregated)$weight,
      mode = c("all")
    )
    
    #correct for small values
    global_eff<-round(global_eff*10000,2)
    local_eff<-round(local_eff*10000,2)
    
    #name local efficiency
    local_eff_collection <- setNames(
      c(local_eff),
      present_ants
    )
    
    local_eff_full <- setNames(
      rep(NA_real_, length(expected_ants)),
      expected_ants
    )
    
    
    # fill in values for ants that exist
    local_eff_full[names(local_eff_collection)] <- local_eff_collection
    
    #assign to collection
    local_efficiency[[paste0(colony_name,"_",j)]]<-local_eff_full
    
    #colony level: mean distance of total graph
    mean_distance<-mean_distance(
      graph_aggregated,
      weights = E(graph_aggregated)$weight,
      unconnected = TRUE,
      details = FALSE
      
    )
    
    # antlevel eigenvector-centrallity
    centrality<-eigen_centrality(
      graph_aggregated,
      directed = FALSE,
      scale = TRUE,
      weights = E(graph_aggregated)$weight,
      options = arpack_defaults()
    )
    
    #restore ant order and save externally
    centrality_collectiion<-centrality$vector
    
    centrality_collection<-setNames(
      centrality$vector,
      present_ants
    )
    
    #pad to all expected ants
    centrality_full <- setNames(
      rep(NA_real_, length(expected_ants)),
      expected_ants
    )
    
    # fill in values for ants that exist
    centrality_full[names(centrality_collection)] <- centrality_collection
    
    #assign to extrernal collection
    centrality_antlevel_df[[paste0(colony_name,"_",j)]]<-centrality_full
    
    ########
    #clustering_fast and greedy
    
    #cluster_fast_greedy(graph_aggregated, weights = E(graph_aggregated)$weight)
    #cluster_louvain(graph_aggregated, weights = E(graph_aggregated)$weight, resolution = 0.5)
    #cluster_optimal(graph_aggregated, weights = E(graph_aggregated)$weight)
    
    #TEMPORAL MEASURES
    #colony level measure interaction durations
    interaction_length_mean<-mean(network_obj_5$duration)/framerate_col
    interaction_length_sd<-sd(network_obj_5$duration)/framerate_col
    
    #colonylevel measure waiting time mean and sd
    waiting_times<-data.frame(wt=diff(sort(network_obj_5$onset)) )
    wait <- waiting_times %>%
      filter(wt > 0)
    waiting_times <- waiting_times[waiting_times$wt > 0, , drop = FALSE]
    waiting_times_5<-data.frame(wt=diff(sort(network_obj_5$onset)))
    waiting_times_5<-waiting_times_5[waiting_times_5$wt>0,,drop = FALSE]/framerate_col
    max_wt<- max(wait$wt)
    max_wt_5<- max(waiting_times_5$wt)
    
    waiting_time_mean<-mean(unlist(waiting_times_5))
    waiting_time_sd<-sd(unlist(waiting_times_5))
    #burstyness
    burstiness<-(waiting_time_sd-waiting_time_mean)/(waiting_time_sd+waiting_time_mean)
    
    #################
    #collect and write result_networkparameters
    results_networkparameters_allintervalls<-results_networkparameters_allintervalls %>%
      add_row(colony_name = colony_name,
              time_interval= j,
              number_ants = number_ants, 
              total_number_interactions = total_number_interactions,
              total_sum_of_interactions = total_sum_of_interactions,
              strength_mean= strength_mean ,
              strength_sd = strength_sd,
              density_aggr = density_aggr,
              interaction_length_mean = interaction_length_mean,
              interaction_length_sd = interaction_length_sd,
              waiting_time_mean = waiting_time_mean,
              waiting_time_sd = waiting_time_sd,
              assortativity=ass_ants,
              global_eff=global_eff,
              mean_distance = mean_distance,
              burstiness = burstiness,
              clustering_global= clustering_global
              
      )
    #add timestamp and folder and save
    write.csv(results_networkparameters_allintervalls,
              file = paste0(folder_path_networkpara, colony_name, "_network_parameters_allintervalls_", timestamp, ".csv"),
              row.names = FALSE
    )
    #progress check
    print(paste0(colony_name, " done", j ))
    
  }
  
}

##########################
#save antlevel data
write.csv(clustering_local_df,
          file = paste0(folder_path_networkpara, "ANTLEVEL_allintervalls_clustering_local", timestamp, ".csv"),
          row.names = FALSE)
write.csv(strength_antlevel_mean,
          file = paste0(folder_path_networkpara, "ANTLEVEL_allintervalls_strength_mean", timestamp, ".csv"),
          row.names = FALSE)
write.csv(local_efficiency,
          file = paste0(folder_path_networkpara, "ANTLEVEL_allintervalls_local_efficiency", timestamp, ".csv"),
          row.names = FALSE)
write.csv(waiting_time_ant_mean_df,
          file = paste0(folder_path_networkpara, "ANTLEVEL_allintervalls_waitingtime_mean", timestamp, ".csv"),
          row.names = FALSE)
write.csv(waiting_time_ant_sd_df,
          file = paste0(folder_path_networkpara, "ANTLEVEL_allintervalls_waitingtime_sd", timestamp, ".csv"),
          row.names = FALSE)
write.csv(interactionlength_ant_mean,
          file = paste0(folder_path_networkpara, "ANTLEVEL_allintervalls_interactionlength_mean", timestamp, ".csv"),
          row.names = FALSE)
write.csv(burstiness_ant_df,
          file = paste0(folder_path_networkpara, "ANTLEVEL_allintervalls_burstiness", timestamp, ".csv"),
          row.names = FALSE)
write.csv(centrality_antlevel_df,
          file = paste0(folder_path_networkpara, "ANTLEVEL_allintervalls_centrality", timestamp, ".csv"),
          row.names = FALSE)
##########################

##############################################
#plot time component


#read data
#first 75min
results_networkparameters_allintervalls<-read.csv(file = "network_parameters/exp1_12hs/5_mins/1_15/ba16-8_network_parameters_allintervalls_26022026.csv")
#2nd 75 min
results_networkparameters_allintervalls<-read.csv(file = "network_parameters/exp1_12hs/5_mins/16_30/ba16-8_network_parameters_allintervalls_26022026.csv")
#after 10 hs
results_networkparameters_allintervalls<-read.csv(file = "network_parameters/exp1_12hs/5_mins/106_120/ba16-8_network_parameters_allintervalls_27022026.csv")
#last 75 mins after 11hs
results_networkparameters_allintervalls<-read.csv(file = "network_parameters/exp1_12hs/5_mins/120_135/ba16-8_network_parameters_allintervalls_27022026.csv")

#remove NA line, first line
results_networkparameters_allintervalls <-
  results_networkparameters_allintervalls[-1, ]
#sort data for glmm, add treatment line 
df_networkparameters_allintervalls <- results_networkparameters_allintervalls %>%
  mutate(
    treatment = sub("\\d+.*$", "", colony_name)   # extract 'a', 'b', or 'ba'
  )







ggplot(df_networkparameters_allintervalls,
       aes(x = time_interval,
           y = strength_mean,
           color = colony_name,
           group = treatment)) +
  geom_point(size = 2) +
  scale_color_manual(
    values = colony_colors
  ) +
  geom_line(
    alpha = 0.5
  )+
  labs(
    x = "Time interval",
    y = "Mean strength",
    color = "Treatment"
  ) +
  theme_minimal()

####################
####################
#time development per treatment and per colony_colors

treatments<-c("a", "b", "ba")
folder_path_timeplots<-"network_parameter_plots/exp1_12hs/5_mins/time_continuum/overlays/"
for (param in parameters){
  for (treat in treatments){
    
    param_plot<-df_networkparameters_allintervalls |>
      filter(treatment == treat) |>
      ggplot(
        aes(x = time_interval,
            y = .data[[param]],
            color = colony_name,
            group = colony_name)) +
      geom_point(size = 2) +
      scale_color_manual(
        values = colony_colors
      ) +
      geom_line(
        alpha = 0.5
      )+
      labs(
        x = "time [5min intervalls]",
        y = param,
        color = "Treatment"
      ) +
      theme_minimal()
    
    ggsave(
      filename = paste0(folder_path_timeplots, "/", param,"_",treat, ".png"),
      plot = param_plot,
      width = 10,
      height = 5,
      dpi = 300
    )
  }
}

####################
#overlay 
#mean over all colonies for every timepoint, use single upper and lower bounds and geomribbon to indicate variance
####################

for (param in parameters){
  summary_df <- df_networkparameters_allintervalls |>
    group_by(treatment, time_interval) |>
    summarise(
      mean_val = mean(.data[[param]], na.rm = TRUE),
      sd_val   = sd(.data[[param]], na.rm = TRUE),
      n             = n(),
      se_val   = sd_val / sqrt(n),
      .groups = "drop"
    )
  
  overlay_plot<-ggplot(summary_df,
                       aes(x = time_interval,
                           y = mean_val,
                           color = treatment,
                           fill  = treatment,
                           group = treatment)) +
    geom_line(linewidth = 1) +
    geom_ribbon(
      aes(
        ymin = mean_val - se_val,
        ymax = mean_val + se_val
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
    labs(
      x = "time [5min intervalls]",
      y = param,
      color = "Treatment",
      fill  = "Treatment"
    ) +
    theme_minimal()
  
  ggsave(
    filename = paste0(folder_path_timeplots, "/", param,"_","overlay.png"),
    plot = overlay_plot,
    width = 10,
    height = 5,
    dpi = 300
  )
  
}
###################
#compare network before and after infection
###################
library(patchwork)
#read data
 files<-c("network_parameters/exp1_12hs/5_mins/1_15/ba16-8_network_parameters_allintervalls_26022026.csv", #first 75mins
          "network_parameters/exp1_12hs/5_mins/16_30/ba16-8_network_parameters_allintervalls_26022026.csv", #2nd 75mins
          "network_parameters/exp1_12hs/5_mins/106_120/ba16-8_network_parameters_allintervalls_27022026.csv", #after 10hs,
          "network_parameters/exp1_12hs/5_mins/120_135/ba16-8_network_parameters_allintervalls_27022026.csv",#last 75mins after 11 hs 
          "network_parameters/exp1_inf/5mins/ba16-8_network_parameters_allintervalls_27022026.csv"#infection data
          )

before_afterinfection<-lapply(files, read.csv)

before_afterinfection<-lapply(before_afterinfection, function(df){
  df<-df[-1,]#remove NA in first line
  
  df<-df%>% mutate(
    treatment = sub("\\d+.*$", "", colony_name))   # extract 'a', 'b', or 'ba'
    
  df<-df%>%filter(time_interval != 15)
  
})

## add time point marker to certain files to keep them apart
before_afterinfection[[1]]<-before_afterinfection[[1]]%>%mutate(globaltime = "0-1.25h")
before_afterinfection[[2]]<-before_afterinfection[[2]]%>%mutate(globaltime = "1.25-2.5h")
before_afterinfection[[3]]<-before_afterinfection[[3]]%>%mutate(globaltime = "8.75h-10h")
before_afterinfection[[4]]<-before_afterinfection[[4]]%>%mutate(globaltime = "10h-11.25h")
before_afterinfection[[5]]<-before_afterinfection[[5]]%>%mutate(globaltime = "fungus exp")

#introduce marker for before and after infection
before_afterinfection[[1]]<-before_afterinfection[[1]]%>%mutate(infection = "before")
before_afterinfection[[2]]<-before_afterinfection[[2]]%>%mutate(infection = "before")
before_afterinfection[[3]]<-before_afterinfection[[3]]%>%mutate(infection = "before")
before_afterinfection[[4]]<-before_afterinfection[[4]]%>%mutate(infection = "before")
before_afterinfection[[5]]<-before_afterinfection[[5]]%>%mutate(infection = "after")

#merge all
before_after_infection_merged<-bind_rows(before_afterinfection)

##filter datasets by treatment
before_after_infection_a<-before_after_infection_merged %>% filter(treatment == "a")
before_after_infection_b<-before_after_infection_merged %>% filter(treatment == "b")
before_after_infection_ba<-before_after_infection_merged %>% filter(treatment == "ba")

# look at each treatment seperately before and after infection for each parameter
global_times<-c("0-1.25h", "1.25-2.5h","8.75h-10h", "10h-11.25h", "fungus exp")
global_times_colors_a<-c("0-1.25h" = "cyan", 
                       "1.25-2.5h" = "slateblue2",
                       "8.75h-10h" = "blue",
                       "10h-11.25h" = "dodgerblue1",
                       "fungus exp"= "springgreen1")


global_times_colors_b<-c("0-1.25h" = "firebrick2", 
                       "1.25-2.5h" = "red4",
                       "8.75h-10h" = "tomato2",
                       "10h-11.25h" = "darkorange",
                       "fungus exp"= "springgreen1")

global_times_colors_ba<-c("0-1.25h" = "darkviolet", 
                       "1.25-2.5h" = "magenta",
                       "8.75h-10h" = "plum2",
                       "10h-11.25h" = "orchid",
                       "fungus exp"= "springgreen1")

for (param in parameters){
 
  summary_df <- before_after_infection_ba|>
    group_by(globaltime, time_interval) |>
    summarise(
      mean_val = mean(.data[[param]], na.rm = TRUE),
      sd_val   = sd(.data[[param]], na.rm = TRUE),
      n             = n(),
      se_val   = sd_val / sqrt(n),
      .groups = "drop"
    )
  
  overlay_plot<-ggplot(summary_df,
                       aes(x = time_interval,
                           y = mean_val,
                           color = globaltime,
                           fill  = globaltime,
                           group = globaltime)) +
    geom_line(linewidth = 1) +
    geom_ribbon(
      aes(
        ymin = mean_val - se_val,
        ymax = mean_val + se_val
      ),
      alpha = 0.25,
      color = NA
    ) +
    scale_color_manual(
      values = global_times_colors
    ) +
    scale_fill_manual(
      values = global_times_colors
    ) +
    coord_cartesian(ylim = c(100, 1300))+
    labs(
      x = "time [5min intervalls]",
      y = param,
      color = "Time",
      fill  = "Time"
    ) +
    theme_minimal()
  
  ggsave(
    filename = paste0(folder_path_timeplots, "/comparison/", param,"_","overlaybeforeafterinf_ba.png"),
    plot = overlay_plot,
    width = 10,
    height = 5,
    dpi = 300
  )
  
}

###########################
#same plot of comparisons all in one
###########################
folder_path_timeplots<-"network_parameter_plots/exp1_12hs/5_mins/time_continuum/overlays/"

for (param in parameters){

  summary_a <- before_after_infection_a|>
    group_by(globaltime, time_interval) |>
    summarise(
      mean_val = mean(.data[[param]], na.rm = TRUE),
      sd_val   = sd(.data[[param]], na.rm = TRUE),
      n             = n(),
      se_val   = sd_val / sqrt(n),
      .groups = "drop"
    )
  summary_b <- before_after_infection_b|>
    group_by(globaltime, time_interval) |>
    summarise(
      mean_val = mean(.data[[param]], na.rm = TRUE),
      sd_val   = sd(.data[[param]], na.rm = TRUE),
      n             = n(),
      se_val   = sd_val / sqrt(n),
      .groups = "drop"
    )
  summary_ba <- before_after_infection_ba|>
    group_by(globaltime, time_interval) |>
    summarise(
      mean_val = mean(.data[[param]], na.rm = TRUE),
      sd_val   = sd(.data[[param]], na.rm = TRUE),
      n             = n(),
      se_val   = sd_val / sqrt(n),
      .groups = "drop"
    )
  
  min_range<- min(c(summary_a$mean_val-summary_a$se_val, summary_b$mean_val-summary_b$se_val, summary_ba$mean_val-summary_ba$se_val))
  
  
  max_range<- max(c(summary_a$mean_val+summary_a$se_val, summary_b$mean_val+summary_b$se_val, summary_ba$mean_val+summary_ba$se_val))
 
  
  comparisonplot_a<-ggplot(summary_a,
                       aes(x = time_interval,
                           y = mean_val,
                           color = globaltime,
                           fill  = globaltime,
                           group = globaltime)) +
    geom_line(linewidth = 1) +
    geom_ribbon(
      aes(
        ymin = mean_val - se_val,
        ymax = mean_val + se_val
      ),
      alpha = 0.25,
      color = NA
    ) +
    scale_color_manual(
      values = global_times_colors_a
    ) +
    scale_fill_manual(
      values = global_times_colors_a
    ) +
    coord_cartesian(ylim = c(min_range, max_range))+
    labs(
      x = "time increments [5 min]",
      y = param,
      title = "a",
      color = "extracted from:",
      fill  = "extracted from:"
    ) +
    theme_minimal()+
    theme(
      legend.position = c(0.98, 0.98),     # inside top-right
      legend.justification = c(1, 1),      # anchor legend at its top-right
       legend.text  = element_text(size = 10),
      legend.key.size = unit(0.6, "cm")    # size of legend keys
    )
  
  comparisonplot_b<-ggplot(summary_b,
                           aes(x = time_interval,
                               y = mean_val,
                               color = globaltime,
                               fill  = globaltime,
                               group = globaltime)) +
    geom_line(linewidth = 1) +
    geom_ribbon(
      aes(
        ymin = mean_val - se_val,
        ymax = mean_val + se_val
      ),
      alpha = 0.25,
      color = NA
    ) +
    scale_color_manual(
      values = global_times_colors_b
    ) +
    scale_fill_manual(
      values = global_times_colors_b
    ) +
    coord_cartesian(ylim = c(min_range, max_range))+
    labs(
      x = "time increments [5 min]",
      title = "b",
      y = param,
      color = "extracted from:",
      fill  = "extracted from:"
    ) +
    theme_minimal()+
    theme(
      legend.position = c(0.98, 0.98),     # inside top-right
      legend.justification = c(1, 1),      # anchor legend at its top-right
      legend.text  = element_text(size = 10),
      legend.key.size = unit(0.6, "cm")    # size of legend keys
    )
  
  comparisonplot_ba<-ggplot(summary_ba,
                           aes(x = time_interval,
                               y = mean_val,
                               color = globaltime,
                               fill  = globaltime,
                               group = globaltime)) +
    geom_line(linewidth = 1) +
    geom_ribbon(
      aes(
        ymin = mean_val - se_val,
        ymax = mean_val + se_val
      ),
      alpha = 0.25,
      color = NA
    ) +
    scale_color_manual(
      values = global_times_colors_ba
    ) +
    scale_fill_manual(
      values = global_times_colors_ba
    ) +
    coord_cartesian(ylim = c(min_range, max_range))+
    labs(
      x = "time increments [5 min]",
      y = param,
      title = "ba",
      color = "extracted from:",
      fill  = "extracted from:"
    ) +
    theme_minimal()+
    theme(
      legend.position = c(0.98, 0.98),     # inside top-right
      legend.justification = c(1, 1),      # anchor legend at its top-right
      legend.text  = element_text(size = 10),
      legend.key.size = unit(0.6, "cm")    # size of legend keys
    )
  #combine into one plot with same legend
  
  grid_plot<-(comparisonplot_a| comparisonplot_b| comparisonplot_ba) +
    plot_layout( axes = "collect" )
  # & theme(legend.position = "bottom")
  
  ggsave(
    filename = paste0(folder_path_timeplots, "/comparison/", param,"_","overlaybeforeafterinf_ALL.png"),
    plot = grid_plot,
    width = 12,
    height = 5,
    dpi = 300
  )
  
}


#############################
#############################
library(lme4)

hist (before_after_infection_a$strength_mean)

model<-lmer(data = before_after_infection_a, formula = strength_mean ~ globaltime + (1|colony_name))
model<-glm

emmeans(model, pairwise ~ globaltime)

####################################
#glmm compare after infection period across treatments
####################################
 infection<-before_after_infection_merged%>%filter(globaltime =="fungus exp")
for (param in parameters){
  hist (before_after_infection_merged$param )
}
#strength_mean
hist (before_after_infection_merged$strength_mean )
model_g<- glmmTMB(
strength_mean ~ treatment + (1 | colony_name),
 dispformula =~ 1,
  data = before_after_infection_merged, 
  #family = gaussian()
  family = Gamma(link = "log")
)

model_g<- glmmTMB(
  log(strength_mean) ~ treatment + (1 | colony_name),
  dispformula =~ treatment,
  data = before_after_infection_merged, 
  family = gaussian()
  
)

drop1(model_g, test = "Chisq")

sim <- DHARMa::simulateResiduals(model_g)
plot(sim)


AIC(model_g)

emm_trt <- emmeans(model_gt, ~ genotype*puremixed)
pairs(emm_trt, adjust = "tukey")

emmeans(model_g, pairwise ~ treatment)
emmeans(model_g, pairwise ~ time_interval | treatment)

#########
#glmm compare only infection period 
#########
infection<-before_after_infection_merged%>%filter(globaltime =="fungus exp")

for (param in parameters){
  png(file=paste0("network_parameter_plots/exp1_inf/5_mins/time_continuum/overlays/hist/",param, ".png"),
      width=600, height=350)
  hist(infection[[param]], xlab = param, col = "purple", main= "")
  dev.off()
}

#strength mean
model_g<- glmmTMB(
  strength_mean ~ treatment + (1 | colony_name),
  dispformula =~ 1,
  data = infection, 
  #family = gaussian()
  family = Gamma(link = "log")
)

sim <- DHARMa::simulateResiduals(model_g)
plot(sim)

emmeans(model_g, pairwise ~ treatment)

#strength sd
hist(infection$strength_sd)
model_g<- glmmTMB(
  strength_sd ~ treatment + (1 | colony_name),
  dispformula =~ 1,
  data = infection, 
  #family = gaussian()
  family = Gamma(link = "log")
)

sim <- DHARMa::simulateResiduals(model_g)
plot(sim)

emmeans(model_g, pairwise ~ treatment)

#assortativity
hist(infection$assortativity)
model_g<- glmmTMB(
  assortativity ~ treatment + (1 | colony_name),
  #dispformula =~ treatment,
  data = infection, 
  #family = gaussian()
  family = t_family()
)

sim <- DHARMa::simulateResiduals(model_g)
plot(sim)

emmeans(model_g, pairwise ~ treatment)

#burstiness
hist(infection$burstiness)
model_g<- glmmTMB(
  burstiness ~ treatment + (1 | colony_name),
  dispformula =~ 1,
  data = infection, 

  family = t_family()
)

sim <- DHARMa::simulateResiduals(model_g)
plot(sim)

emmeans(model_g, pairwise ~ treatment)
# global effectiveness
infection$global_eff_corr<-infection$global_eff/10000
hist(infection$global_eff_corr)
model_g<- glmmTMB(
  global_eff_corr ~ treatment + (1 | colony_name),
  #dispformula =~ 1,
  data = infection, 
  
  family = gaussian()
)

model_b<- glmmTMB(
  global_eff_corr ~ treatment + (1 | colony_name),
  #dispformula =~ 1,
  data = infection, 
  
  family = beta_family()
)

sim <- DHARMa::simulateResiduals(model_b)
plot(sim)

AIC(model_g, model_b)

emmeans(model_b, pairwise ~ treatment)

#density aggr
hist(infection$density_aggr)
model_g<- glmmTMB(
density_aggr ~ treatment + (1 | colony_name),
  #dispformula =~ 1,
  data = infection, 
  
  family = Gamma(link = "log")
)

sim <- DHARMa::simulateResiduals(model_g)
plot(sim)



emmeans(model_g, pairwise ~ treatment)
#mean distance
hist(infection$mean_distance)
model_g<- glmmTMB(
  mean_distance ~ treatment + (1 | colony_name),
 
  data = infection, 
  
  family = Gamma(link = "log")
)

sim <- DHARMa::simulateResiduals(model_g)
plot(sim)



emmeans(model_g, pairwise ~ treatment)

#interaction length mean
hist(infection$interaction_length_mean)
model_g<- glmmTMB(
  interaction_length_mean ~ treatment + (1 | colony_name),
  dispformula =~ 1,
  data = infection, 
  
  family = t_family()
)

sim <- DHARMa::simulateResiduals(model_g)
plot(sim)



emmeans(model_g, pairwise ~ treatment)
#interaction length_ sd
hist(infection$interaction_length_sd)
model_g<- glmmTMB(
  interaction_length_sd ~ treatment + (1 | colony_name),
  dispformula =~ 1,
  data = infection, 
  #family = gaussian()
  family = t_family()
)

sim <- DHARMa::simulateResiduals(model_g)
plot(sim)



emmeans(model_g, pairwise ~ treatment)

#waiting time mean
hist(infection$waiting_time_mean)
model_g<- glmmTMB(
  waiting_time_mean ~ treatment + (1 | colony_name),
  # dispformula =~ 1,
  data = infection, 
  #family = gaussian()
  family = Gamma(link = "log")
)

sim <- DHARMa::simulateResiduals(model_g)
plot(sim)



emmeans(model_g, pairwise ~ treatment)
# waiting time sd
hist(infection$waiting_time_sd)
model_g<- glmmTMB(
  waiting_time_sd ~ treatment + (1 | colony_name),
  # dispformula =~ 1,
  data = infection, 
  #family = gaussian()
  family = Gamma(link = "log")
)

sim <- DHARMa::simulateResiduals(model_g)
plot(sim)



emmeans(model_g, pairwise ~ treatment)
#global_clustering

hist(infection$clustering_global)
infection$clustering_global_corr<-infection$clustering_global-0.0005
model_g<- glmmTMB(
  clustering_global_corr ~ treatment + (1 | colony_name),
  # dispformula =~ 1,
  data = infection, 
  family = t_family()
  #family = Gamma(link = "log")
)

sim <- DHARMa::simulateResiduals(model_g)
plot(sim)



emmeans(model_g, pairwise ~ treatment)

################################
#difference between treatments no fungal infection
################################
for (param in parameters){
png(file= paste0("network_parameter_plots/exp1_12hs/5_mins/time_continuum/overlays/hist/", param, ".png"),
    width=600, height=350)
hist(no_infection[[param]], xlab = param , col = "purple", main="")
dev.off()
}

no_infection<-before_after_infection_merged%>%filter(globaltime != "fungal exp")
#mean over 5 min time interval
no_infection_mean<-no_infection%>%group_by(colony_name, globaltime, treatment)%>%summarise(
  across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
  .groups = "drop"
)

for (param in parameters){
  png(file= paste0("network_parameter_plots/exp1_12hs/5_mins/time_mean/hist/", param, ".png"),
      width=600, height=350)
  hist(no_infection[[param]], xlab = param , col = "purple", main="")
  dev.off()
}
#strength  mean over mean of 
hist(no_infection_mean$strength_mean)
#g1 bigger AIC
model_g1<- glmmTMB(
  log(strength_mean) ~ treatment + (1 | colony_name),
  #dispformula =~globaltime*,
  data = no_infection_mean, 
  family = gaussian()
  #family = Gamma(link = "log")
)

drop1(model_g1, test = "Chisq")

sim <- DHARMa::simulateResiduals(model_g1)
plot(sim)

emmeans(model_g1, pairwise ~ treatment)

model_g2<- glmmTMB(
  strength_mean ~ treatment + (1 | colony_name),
  #ziformula = ~1,
  #dispformula =~treatment,
  data = no_infection_mean, 
  family = Gamma(link = log)
  #family = Gamma(link = "log")
)

drop1(model_g2, test = "Chisq")
AIC(model_g1, model_g2)

sim <- DHARMa::simulateResiduals(model_g2)
plot(sim)

sim <- DHARMa::simulateResiduals(model_g)
plot(sim)

AIC(model_g, model_g1, model_g2)

#strength_sd
hist(no_infection_mean$strength_sd)
model_g<- glmmTMB(
  log(strength_sd) ~ treatment + (1| colony_name),
 #ziformula = ~1,
   dispformula =~ treatment,
  data = no_infection_mean, 
  # family = Gamma(link = "log")
  family = gaussian()
)
sim <- DHARMa::simulateResiduals(model_g)
plot(sim)

model_g1<- glmmTMB(
  strength_sd ~ treatment + (1| colony_name),
  #ziformula = ~1,
  dispformula =~ treatment,
  data = no_infection_mean, 
  # family = Gamma(link = "log")
  family = gaussian()
)

sim <- DHARMa::simulateResiduals(model_g1)
plot(sim)
AIC(model_g, model_g1)

drop1(model_g, test = "Chisq")
emmeans(model_g, pairwise ~ treatment)


#assortativity

hist(no_infection_mean$assortativity)
model_g<- glmmTMB(
  assortativity ~ treatment + (1 | colony_name),
  #dispformula =~ treatment,
  data = no_infection_mean, 
  #family = gaussian()
  family = t_family()
)

sim <- DHARMa::simulateResiduals(model_g)
plot(sim)

model_g1<- glmmTMB(
  assortativity ~ treatment + (1 | colony_name),
  #dispformula =~ treatment,
  data = no_infection_mean, 
  family = gaussian()
  
)

sim <- DHARMa::simulateResiduals(model_g1)
plot(sim)
drop1(model_g, test= "Chisq")
drop1(model_g1, test= "Chisq")
AIC(model_g, model_g1)
emmeans(model_g, pairwise ~ treatment)

#burstiness
hist(no_infection_mean$burstiness)
model_g<- glmmTMB(
  burstiness ~ treatment + (1 | colony_name),
  # ziformula = ~1,
  # dispformula =~ 1,
  data = no_infection_mean, 
  
  family = gaussian()
)

sim <- DHARMa::simulateResiduals(model_g)
plot(sim)

drop1(model_g, test = "Chisq")
emmeans(model_g, pairwise ~ treatment)
# global effectiveness
no_infection_mean$global_eff_corr<-no_infection_mean$global_eff/10000
hist(no_infection_mean$global_eff_corr)
model_g<- glmmTMB(
  global_eff_corr ~ treatment + (1 | colony_name),
  #dispformula =~ 1,
  data = no_infection_mean, 
  
  family = gaussian()
)
sim <- DHARMa::simulateResiduals(model_g)
plot(sim)

model_g1<- glmmTMB(
  global_eff_corr ~ treatment + (1 | colony_name),
  data = no_infection_mean, 
  family = t_family()
)

sim <- DHARMa::simulateResiduals(model_g1)
plot(sim)
AIC(model_g1, model_g)

drop1(model_g, test = "Chisq")
emmeans(model_g, pairwise ~treatment)


AIC(model_g, model_b)

emmeans(model_b, pairwise ~ treatment)

#density aggr
hist(no_infection_mean$density_aggr)
model_g<- glmmTMB(
  density_aggr ~ treatment + (1|colony_name),
  data = no_infection_mean, 
 
  family = Gamma(link = "log")
)

sim <- DHARMa::simulateResiduals(model_g)
plot(sim)

drop1(model_g, test = "Chisq")



emmeans(model_g, pairwise ~ treatment)

#mean distance
hist(no_infection_mean$mean_distance)
model_g<- glmmTMB(
  mean_distance ~ treatment + (1 | colony_name),
  data = no_infection_mean, 
  family = Gamma(link = "log")
)

drop1(model_g, test = "Chisq")



emmeans(model_g, pairwise ~ treatment)

#interaction length mean
hist(no_infection_mean$interaction_length_mean)
model_g<- glmmTMB(
 log( interaction_length_mean) ~ treatment + (1 | colony_name),

  data = no_infection_mean, 
  
  family = gaussian()
)

sim <- DHARMa::simulateResiduals(model_g)
plot(sim)
drop1(model_g, test= "Chisq")



#interaction length_ sd
hist(no_infection_mean$interaction_length_sd)
#works wayy better, lower AIC
model_g<- glmmTMB(
  log(interaction_length_sd) ~ treatment + (1 | colony_name),
  data = no_infection_mean, 
  family = gaussian()
)

sim <- DHARMa::simulateResiduals(model_g)
plot(sim)

drop1(model_g, test= "Chisq")

#waiting time mean
hist(no_infection_mean$waiting_time_mean)
model_g<- glmmTMB(
  waiting_time_mean ~ treatment + (1 | colony_name),
   #dispformula =~ 1,
  # ziformula = ~0,
  data = no_infection_mean,
  family = Gamma(link = "log")
 
)

sim <- DHARMa::simulateResiduals(model_g)
plot(sim)

drop1(model_g, test = "Chisq")
emmeans(model_g, pairwise ~ treatment)
# waiting time sd
hist(no_infection_mean$waiting_time_sd)
model_g<- glmmTMB(
  waiting_time_sd ~ treatment + (1 | colony_name),
  data = no_infection_mean, 
  family = Gamma(link="log")
)

sim <- DHARMa::simulateResiduals(model_g)
plot(sim)

drop1(model_g, test= "Chisq")

emmeans(model_g, pairwise ~ treatment)
#global_clustering
no_infection_mean$clustering_global_corr<-no_infection_mean$clustering_global-0.0005
hist(no_infection_mean$clustering_global)

model_g<- glmmTMB(
  clustering_global ~ treatment + (1 | colony_name),
  data = no_infection_mean, 
  family = gaussian()
)

sim <- DHARMa::simulateResiduals(model_g)
plot(sim)


drop1(model_g, test = "Chisq")
emmeans(model_g, pairwise ~ treatment)
############
#comparing per treatment before and after infection
############
#add infection marker

before_after_infection_mean<-before_after_infection_merged%>%group_by(globaltime, colony_name, treatment, infection)%>%summarise(across(where(is.numeric),~mean(.x, na.rm =TRUE)), .groups = "drop")
#pre_post_comparison$infection<-ifelse(pre_post_comparison_a$globaltime == "fungus exp", "post", "pre")

for (param in parameters){
  png(file= paste0("network_parameter_plots/exp1_12hs/5_mins/time_mean/hist_infection/", param, ".png"),
      width=600, height=350)
  hist(before_after_infection_mean[[param]], xlab = param , col = "purple", main="")
  dev.off()
}
##strength mean
hist(before_after_infection_mean$strength_mean)
model_g<- glmmTMB(
  log(strength_mean) ~ infection*treatment + (1 | colony_name),
  data = before_after_infection_mean,
  family = gaussian()
)

sim <- DHARMa::simulateResiduals(model_g)
plot(sim)


drop1(model_g, test = "Chisq")
emmeans(model_g, pairwise ~ infection|treatment, type = "response")

#strength_sd

hist(before_after_infection_mean$strength_sd)
model_g<- glmmTMB(
  log(strength_sd) ~ infection*treatment + (1 | colony_name),
  data = before_after_infection_mean,
  family = gaussian()
)

model_g1<- glmmTMB(
  strength_sd ~ infection*treatment + (1 | colony_name),
  data = before_after_infection_mean,
  family = gaussian()
)

sim <- DHARMa::simulateResiduals(model_g1)
plot(sim)

AIC(model_g, model_g1)

drop1(model_g, test = "Chisq")
emmeans(model_g, pairwise ~ infection|treatment, type = "response")
## assortativity
hist(before_after_infection_mean$assortativity)
model_g<- glmmTMB(
  assortativity ~ infection*treatment + (1 | colony_name),
  data = before_after_infection_mean,
  family = gaussian()
)
sim <- DHARMa::simulateResiduals(model_g)
plot(sim)

model_g1<- glmmTMB(
  assortativity ~ infection+treatment + (1 | colony_name),
  data = before_after_infection_mean,
  family = gaussian()
)

sim <- DHARMa::simulateResiduals(model_g1)
plot(sim)



drop1(model_g, test = "Chisq")
drop1(model_g1, test = "Chisq")
emmeans(model_g1, pairwise ~ infection|treatment, type = "response")
###burstiness
hist(before_after_infection_mean$burstiness)
model_g<- glmmTMB(
 burstiness ~ infection*treatment + (1 | colony_name),
  data = before_after_infection_mean,
  family = gaussian()
)
sim <- DHARMa::simulateResiduals(model_g)
plot(sim)

model_g1<- glmmTMB(
  burstiness ~ infection+treatment + (1 | colony_name),
  data = before_after_infection_mean,
  family = gaussian()
)

sim <- DHARMa::simulateResiduals(model_g1)
plot(sim)

drop1(model_g, test = "Chisq")
drop1(model_g1, test = "Chisq")
emmeans(model_g1, pairwise ~ infection|treatment, type = "response")

#global_effectiveness
hist(before_after_infection_mean$global_eff_corr)
before_after_infection_mean$global_eff_corr<-before_after_infection_mean$global_eff/10000
model_g<- glmmTMB(
 global_eff_corr ~ infection*treatment + (1 | colony_name),
  data = before_after_infection_mean,
  family = gaussian()
)
sim <- DHARMa::simulateResiduals(model_g)
plot(sim)

drop1(model_g, test = "Chisq")

emmeans(model_g, pairwise ~ infection|treatment, type = "response")
## density aggr
hist(before_after_infection_mean$density_aggr)

model_g<- glmmTMB(
 density_aggr ~ infection*treatment + (1 | colony_name),
  data = before_after_infection_mean,
  family = Gamma(link="log")
)
sim <- DHARMa::simulateResiduals(model_g)
plot(sim)

drop1(model_g, test = "Chisq")

emmeans(model_g, pairwise ~ infection|treatment, type = "response")

#waiting time mean

hist(before_after_infection_mean$waiting_time_mean)

model_g<- glmmTMB(
  waiting_time_mean ~ infection*treatment + (1 | colony_name),
  data = before_after_infection_mean,
  dispformula = ~treatment,
  family = Gamma(link="log")
)
sim <- DHARMa::simulateResiduals(model_g)
plot(sim)

drop1(model_g, test = "Chisq")

model_g1<- glmmTMB(
  waiting_time_mean ~ infection+treatment + (1 | colony_name),
  data = before_after_infection_mean,
  dispformula = ~treatment,
  family = Gamma(link="log")
)
sim <- DHARMa::simulateResiduals(model_g1)
plot(sim)
drop1(model_g1, test = "Chisq")
emmeans(model_g1, pairwise ~ infection|treatment, type = "response")

#waiting time sd
hist(before_after_infection_mean$waiting_time_sd)

model_g<- glmmTMB(
  waiting_time_sd ~ infection*treatment + (1 | colony_name),
  data = before_after_infection_mean,
  dispformula = ~infection,
  #family = gaussian()
  family = Gamma(link= "log")
)
sim <- DHARMa::simulateResiduals(model_g)
plot(sim)

drop1(model_g, test = "Chisq")

model_g1<- glmmTMB(
  waiting_time_mean ~ infection+treatment + (1 | colony_name),
  data = before_after_infection_mean,
  dispformula = ~treatment,
  family = Gamma(link="log")
)
sim <- DHARMa::simulateResiduals(model_g1)
plot(sim)
drop1(model_g1, test = "Chisq")
emmeans(model_g1, pairwise ~ infection|treatment, type = "response")

#mean distance
hist(before_after_infection_mean$mean_distance)

model_g2<- glmmTMB(
  mean_distance ~ infection*treatment + (1 | colony_name),
  data = before_after_infection_mean,
  dispformula = ~ infection,
  
   family = Gamma(link="log")
)
sim <- DHARMa::simulateResiduals(model_g2)
plot(sim)

model_g<- glmmTMB(
  log(mean_distance) ~ infection*treatment + (1 | colony_name),
  data = before_after_infection_mean,
  dispformula = ~ infection,
  family = gaussian()
  # family = Gamma(link="log")
)
sim <- DHARMa::simulateResiduals(model_g)
plot(sim)

AIC(model_g, model_g2)

drop1(model_g, test = "Chisq")

model_g1<- glmmTMB(
  log(mean_distance) ~ infection+treatment + (1 | colony_name),
  data = before_after_infection_mean,
  dispformula = ~ infection,
  family = gaussian()
  # family = Gamma(link="log")
)
sim <- DHARMa::simulateResiduals(model_g1)
plot(sim)
drop1(model_g1, test = "Chisq")

emmeans(model_g1, pairwise ~ infection|treatment, type = "response")

#global clustering
hist(before_after_infection_mean$clustering_global)

model_g<- glmmTMB(
 clustering_global ~ infection*treatment + (1 | colony_name),
  data = before_after_infection_mean,
  dispformula = ~ infection,
  
  family = Gamma(link="log")
)
sim <- DHARMa::simulateResiduals(model_g)
plot(sim)
model_g1<- glmmTMB(
  clustering_global ~ infection*treatment + (1 | colony_name),
  data = before_after_infection_mean,
  
  
  family = beta_family()
)
sim <- DHARMa::simulateResiduals(model_g1)
plot(sim)

AIC(model_g, model_g1)
drop1(model_g1, test = "Chisq")

model_g<- glmmTMB(
  clustering_global ~ infection+treatment + (1 | colony_name),
  data = before_after_infection_mean,
  
  
  family = beta_family()
)
sim <- DHARMa::simulateResiduals(model_g)
plot(sim)
drop1(model_g, test = "Chisq")

#interaction mean
hist(before_after_infection_mean$interaction_length_mean)

model_g<- glmmTMB(
  log(interaction_length_mean) ~ infection*treatment + (1 | colony_name),
  data = before_after_infection_mean,
  family = gaussian
)
sim <- DHARMa::simulateResiduals(model_g)
plot(sim)

drop1(model_g, test = "Chisq")

emmeans(model_g, pairwise ~ infection|treatment, type = "response")

#interaction sd

hist(before_after_infection_mean$interaction_length_sd)

model_g<- glmmTMB(
  log(interaction_length_sd) ~ infection*treatment + (1 | colony_name),
  data = before_after_infection_mean,
  family = gaussian
)
sim <- DHARMa::simulateResiduals(model_g)
plot(sim)

drop1(model_g, test = "Chisq")

emmeans(model_g, pairwise ~ infection|treatment, type = "response")


##plot mean parameters that work on:
infection_status<-c("before" = "gray75",
                    "after" = "springgreen")
for (param in parameters){

  boxplot_means<-ggplot(before_after_infection_mean,
                       aes(x = treatment,
                           y = .data[[param]],
                           color = colony_name,
                           fill  = treatment,
                           group = treatment
                           )) +
    geom_boxplot(aes(fill = treatment), alpha = 0.3)+
  geom_point(aes(color = colony_name), size = 3)+
   
    scale_color_manual(
      values = colony_colors
    ) +
    scale_fill_manual(
      values = treatment_colors
    ) +
   labs(
      x = "treatments",
      y = param,
      color = "colony:",
      fill  = "treatment"
    ) +
    theme_minimal()
  
  ggsave(
    filename = paste0("network_parameter_plots/exp1_12hs/5_mins/time_mean/raw_means/colorby_colony/" , param, ".png"),
    plot = boxplot_means,
    width = 10,
    height = 5,
    dpi = 300
  )
  
}

#####color by infection period

for (param in parameters){
  
  boxplot_means<-ggplot(before_after_infection_mean,
                        aes(x = treatment,
                            y = .data[[param]],
                            color = infection,
                            fill  = treatment,
                            group = treatment
                        )) +
    geom_boxplot(aes(fill = treatment), alpha = 0.3)+
    geom_point(aes(color = infection), size = 2)+
    
    scale_color_manual(
      values = infection_status
    ) +
    scale_fill_manual(
      values = treatment_colors
    ) +
    labs(
      x = "treatments",
      y = param,
      color = "infection:",
      fill  = "treatment"
    ) +
    theme_minimal()
  
  ggsave(
    filename = paste0("network_parameter_plots/exp1_12hs/5_mins/time_mean/raw_means/colorby_infection/" , param, ".png"),
    plot = boxplot_means,
    width = 10,
    height = 5,
    dpi = 300
  )
  
}
#

for (param in parameters){
  
  boxplot_means<-ggplot(before_after_infection_mean,
                        aes(x = treatment,
                            y = .data[[param]],
                            color = infection,
                            fill  = treatment
                        )) +
    geom_boxplot(aes(fill = treatment), alpha = 0.3)+
    geom_point(aes(color = infection), size = 2)+
    
    scale_color_manual(
      values = infection_status
    ) +
    scale_fill_manual(
      values = treatment_colors
    ) +
    labs(
      x = "treatments",
      y = param,
      color = "infection:",
      fill  = "treatment"
    ) +
    theme_minimal()
  
  ggsave(
    filename = paste0("network_parameter_plots/exp1_12hs/5_mins/time_mean/raw_means/separateby_infection/" , param, ".png"),
    plot = boxplot_means,
    width = 10,
    height = 5,
    dpi = 300
  )
  
}



################################
#antlevel data plotting
################################

antlevel_files<-c("network_parameters/exp1_12hs/5_mins/1_15/ANTLEVEL_allintervalls_strength_mean11022026.csv",
                  "network_parameters/exp1_12hs/5_mins/16-30/ANTLEVEL_allintervalls_strength_mean11022026.csv",
                  "network_parameters/exp1_12hs/5_mins/106_120/ANTLEVEL_allintervalls_strength_mean11022026.csv",
                  "network_parameters/exp1_12hs/5_mins/120_135//ANTLEVEL_allintervalls_strength_mean11022026.csv",
                  "network_parameters/exp1_inf/5mins/ANTLEVEL_allintervalls_strength_mean10022026.csv" #infection
                  )

antlevel_data_unmerged<-lapply(antlevel_files, read.csv)

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
library(tidyr)

antlevel_data<- antlevel_data_merged%>%
  pivot_longer(
    cols = -ant,
    names_to = "colony_name",
    values_to = "strength_mean"
  )%>%
 separate(colony_name,
           into = c("colony", "time_interval", "globaltime"),
           sep = "_",
           remove = TRUE)%>%
  mutate(time_interval = as.numeric(time_interval))%>%
  mutate(
    treatment = sub("\\d+.*$", "", colony)   # extract 'a', 'b', or 'ba'
  )%>%
  mutate(
    genotype = case_when(
      treatment == "ba" ~ unlist(genotype_match[ant]),  # lookup in genotype_match
      treatment == "a" ~ "a",                            # all ants in pure a colony
      treatment == "b" ~ "b"                             # all ants in pure b colony
    )
  )%>%
  mutate(
    puremixed = case_when(
      treatment == "ba" ~ "m",  # lookup in genotype_match
      treatment == "a" ~ "p",                            # all ants in pure a colony
      treatment == "b" ~ "p"                             # all ants in pure b colony
    )
  )


#########################
#plot values for all ants
#########################
library(ggrepel)
global_times_colors_a<-c("0-1.25h" = "cyan", 
                         "1.25-2.5h" = "slateblue2",
                         "8.75-10h" = "blue",
                         "10-11.25h" = "dodgerblue1",
                         "fungal exp"= "springgreen1")


output_folder <- "network_parameter_plots/antlevel_1_truncated_2min/"

ants_a <- antlevel_data%>%filter(treatment == "b")
halflist_a<-c("b16.1", "b16.2", "b16.3", "b16.4")
ants_a <- antlevel_data%>%filter(colony %in% halflist_a)
ants_a<-ants_a%>%filter(time_interval !=15)

ants_a <- ants_a %>%
  arrange(colony, ant) %>%
  mutate(colony_ant = interaction(ant,colony, sep = "_"))
ants_a<-ants_a%>%mutate(infected = colony_ant %in% exposed_ants_string)

ggplot(ants_a,
       aes(x = colony_ant,
           y = strength_mean,
           color = globaltime)) +
  geom_point(
    aes(fill = globaltime),
    shape = 21,             # circle with outline
    color = ifelse(ants_a$infected, "red","white"),  # outline
    size = 3,
    stroke = 1.5,           # thickness of outline
    alpha = 0.7
  ) +
  theme_minimal() +
  scale_color_manual(values = global_times_colors_a) +
  labs(
    x = "Ant (within colony)",
    y = "strength_mean [s]"
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggplot(ants_a, aes(x = colony, y = strength_mean, fill = treatment)) +
  
  geom_jitter(aes(color = genotype), width = 0.1, alpha = 0.7, size = 2) +
  geom_boxplot(width = 0.3, outlier.shape = NA, color = "black", alpha = 0.4) +
  theme_minimal() +
  scale_fill_manual(values = global_times_colors_a) +
  scale_color_manual(values = global_times_colors_a) +
  labs(
    title = " ",
    x = "colony",
    y = "strength_mean [s]"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # rotate x labels for readability

ggplot(ants_a, aes(x = colony, y = strength_mean, color= globaltime, group = ant)) +
  geom_point(
    aes(color = globaltime),
    position = position_jitterdodge(
      jitter.width = 0.1,
      dodge.width = 1
    ),
    size = 1,
    alpha = 0.7
  )+
  theme_minimal() +
  scale_fill_manual(values = global_times_colors_a) +
  scale_color_manual(values = global_times_colors_a) +
  labs(
    x = "Colony",
    y = "strength_mean [s]"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# Save each plot as PNG
ggsave(
  filename = paste0(output_folder, "/", "ANTLEVEL_meaninteractionlength.png"),
  plot = param_plot,
  width = 6,
  height = 5,
  dpi = 300
)

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

