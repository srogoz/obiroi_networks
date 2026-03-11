###time sequence for infected in 5 min intervalls
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
source("global_parameters_exp2_12hs.R")

source_folder_edgelist<-"edge_lists/exp2_12hs/"
folder_path_networkpara <-"network_parameters/exp2_12hs/5_mins/121_135/"
output_folder <- "network_parameter_plots/exp2_12hs/5_mins/"
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
network_obj_list<-get_colony_files(source_folder_edgelist, selected_colonies)
#######################################
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
                                                    clustering_global= NA,
                                                    assortativity=NA,
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
    length.out = number_intervalls + time_between + 1
  )
  
  for (j in 1:number_intervalls){
    #start and stop
    #j, add same length if later timeintervalls
    time_window_start<-intervall_timepoints[j+ time_between]
    #j+1 add same length if later timeintervalls
    time_window_end<-intervall_timepoints[j+1+time_between]
    
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
    global_eff<-global_eff
    local_eff<-local_eff
    
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


##############################################
#plot time component


#read data
#first 75min
results_networkparameters_allintervalls<-read.csv(file = "network_parameters/exp1_12hs/5_mins/1_15/ba16-8_network_parameters_allintervalls_11022026.csv")
#2nd 75 min
results_networkparameters_allintervalls<-read.csv(file = "network_parameters/exp1_12hs/5_mins/16_30/ba16-8_network_parameters_allintervalls_11022026.csv")
#after 10 hs
results_networkparameters_allintervalls<-read.csv(file = "network_parameters/exp1_12hs/5_mins/106_120/ba16-8_network_parameters_allintervalls_11022026.csv")
#last 75 mins after 11hs
results_networkparameters_allintervalls<-read.csv(file = "network_parameters/exp1_12hs/5_mins/120_135/ba16-8_network_parameters_allintervalls_11022026.csv")

#remove NA line, first line
results_networkparameters_allintervalls <-
  results_networkparameters_allintervalls[-1, ]
#sort data for glmm, add treatment line 
df_networkparameters_allintervalls <- results_networkparameters_allintervalls %>%
  mutate(
    treatment = sub("\\d+.*$", "", colony_name)   # extract 'a', 'b', or 'ba'
  )


#plot time development of parameters

treatment_colors<-c(
  "a" = "#05e0fc",
  "b" = "#fc0536",
  "ba" = "#db05fc"
)

colony_colors <-c(
  "a16-1" = "blue",
  "a16-2" = "royalblue",
  "a16-3" = "slateblue2",
  "a16-4" = "navyblue",
  "a16-5" = "lightblue2",
  "a16-6" = "deepskyblue1",
  "a16-7"="skyblue2",
  "a16-8"="cyan",
  
  "b16-1" = "firebrick2",
  "b16-2" = "red4",
  "b16-3"= "tomato2",
  "b16-4" = "indianred",
  "b16-5" = "darkorange2",
  "b16-6" = "sienna1",
  "b16-7"="coral",
  "b16-8"="lightsalmon",
  
  "ba16-1"= "darkviolet",
  "ba16-2"="orchid2",
  "ba16-3"="darkmagenta",
  "ba16-4"="magenta1", 
  "ba16-5"="deeppink1",
  "ba16-6"="plum2",
  "ba16-7"="violet",
  "ba16-8"="hotpink2"
)


exposed_ants <-c(
  "a16-1" = "OP",
  "a16-2" = "BP",
  "a16-3" = "OO",
  "a16-4" = "BO",
  "a16-5" = "OG",
  "a16-6" = "BG",
  "a16-7"= "OB",
  "a16-8"= "BB",
  "b16-1" = "PP",
  "b16-2" = "GP",
  "b16-3"= "PO",
  "b16-4" = "GO",
  "b16-5" = "PG",
  "b16-6" = "GG",
  "b16-7"= "PB",
  "b16-8"= "GB",
  "ba16-1"= "OP",
  "ba16-2"="PP",
  "ba16-3"="OG",
  "ba16-4"="GO", 
  "ba16-5"="GG",
  "ba16-6"="PB",
  "ba16-7"="GP",
  "ba16-8"="BB"
)


exposed_ants_string <-c(
  "OP_a16.1",
  "BP_a16.2",
  "OO_a16.3",
  "BO_a16.4",
  "OG_a16.5",
  "BG_a16.6",
  "OB_a16.7",
  "BB_a16.8",
  "PP_b16.1",
  "GP_b16.2",
  "PO_b16.3",
  "GO_b16.4",
  "PG_b16.5",
  "GG_b16.6",
  "PB_b16.7",
  "GB_b16.8",
  "OP_ba16.1",
  "PP_ba16.2",
  "OG_ba16.3",
  "GO_ba16.4", 
  "GG_ba16.5",
  "PB_ba16.6",
  "GP_ba16.7",
  "BB_ba16.8"
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

treatments<-c("B", "b", "bB", "bA")
folder_path_timeplots<-"network_parameter_plots/exp2_12hs/5mins/time_continuum/overlays/"
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
files<-c("network_parameters/exp2_12hs/5_mins/1_15/bA16-8_network_parameters_allintervalls_10032026.csv", #first 75mins
         "network_parameters/exp2_12hs/5_mins/16_30/bA16-8_network_parameters_allintervalls_10032026.csv", #2nd 75mins
         "network_parameters/exp2_12hs/5_mins/106_120/bA16-8_network_parameters_allintervalls_10032026.csv", #after 10hs,
         "network_parameters/exp2_12hs/5_mins/121_135/bA16-8_network_parameters_allintervalls_10032026.csv",#last 75mins after 11 hs 
         "network_parameters/exp2_inf/5mins/bA16-8_network_parameters_allintervalls_10032026.csv"#infection data
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
before_afterinfection[[3]]<-before_afterinfection[[3]]%>%mutate(globaltime = "8.75-10h")
before_afterinfection[[4]]<-before_afterinfection[[4]]%>%mutate(globaltime = "10-11.25h")
before_afterinfection[[5]]<-before_afterinfection[[5]]%>%mutate(globaltime = "fungal exp")

#introduce marker for before and after infection
before_afterinfection[[1]]<-before_afterinfection[[1]]%>%mutate(infection = "before")
before_afterinfection[[2]]<-before_afterinfection[[2]]%>%mutate(infection = "before")
before_afterinfection[[3]]<-before_afterinfection[[3]]%>%mutate(infection = "before")
before_afterinfection[[4]]<-before_afterinfection[[4]]%>%mutate(infection = "before")
before_afterinfection[[5]]<-before_afterinfection[[5]]%>%mutate(infection = "after")

#merge all
before_after_infection_merged<-bind_rows(before_afterinfection)

##filter datasets by treatment
before_after_infection_b<-before_after_infection_merged %>% filter(treatment == "b")
before_after_infection_B<-before_after_infection_merged %>% filter(treatment == "B")
before_after_infection_bB<-before_after_infection_merged %>% filter(treatment == "bB")
before_after_infection_bA<-before_after_infection_merged %>% filter(treatment == "bA")

# look at each treatment separately before and after infection for each parameter

for (param in parameters){
  
  summary_df <- before_after_infection_b|>
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
      values = global_times_colors_b
    ) +
    scale_fill_manual(
      values = global_times_colors_b
    ) +
   # coord_cartesian(ylim = c(100, 1300))+
    labs(
      x = "time [5min intervalls]",
      y = param,
      color = "Time",
      fill  = "Time"
    ) +
    theme_minimal()
  
  ggsave(
    filename = paste0(folder_path_timeplots, "comparison/", param,"_","overlaybeforeafterinf_b.png"),
    plot = overlay_plot,
    width = 10,
    height = 5,
    dpi = 300
  )
  
}

###########################
#same plot of comparisons all in one
###########################

folder_path_timeplots<-"network_parameter_plots/exp2_12hs/5mins/time_continuum/overlays/"
for (param in parameters){
  summary_b <- before_after_infection_b|>
    group_by(globaltime, time_interval) |>
    summarise(
      mean_val = mean(.data[[param]], na.rm = TRUE),
      sd_val   = sd(.data[[param]], na.rm = TRUE),
      n             = n(),
      se_val   = sd_val / sqrt(n),
      .groups = "drop"
    )
  summary_B <- before_after_infection_B|>
    group_by(globaltime, time_interval) |>
    summarise(
      mean_val = mean(.data[[param]], na.rm = TRUE),
      sd_val   = sd(.data[[param]], na.rm = TRUE),
      n             = n(),
      se_val   = sd_val / sqrt(n),
      .groups = "drop"
    )
  summary_bB <- before_after_infection_bB|>
    group_by(globaltime, time_interval) |>
    summarise(
      mean_val = mean(.data[[param]], na.rm = TRUE),
      sd_val   = sd(.data[[param]], na.rm = TRUE),
      n             = n(),
      se_val   = sd_val / sqrt(n),
      .groups = "drop"
    )
  summary_bA <- before_after_infection_bA|>
    group_by(globaltime, time_interval) |>
    summarise(
      mean_val = mean(.data[[param]], na.rm = TRUE),
      sd_val   = sd(.data[[param]], na.rm = TRUE),
      n             = n(),
      se_val   = sd_val / sqrt(n),
      .groups = "drop"
    )
  
  min_range<- min(c(summary_B$mean_val-summary_B$se_val, summary_bA$mean_val-summary_bA$se_val, summary_bB$mean_val-summary_bB$se_val, summary_b$mean_val-summary_b$se_val))
  
  
  max_range<- max(c(summary_bA$mean_val+summary_bA$se_val, summary_b$mean_val+summary_b$se_val, summary_bB$mean_val+summary_bB$se_val, summary_B$mean_val+summary_B$se_val))
  
  
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
  
  
  comparisonplot_B<-ggplot(summary_B,
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
      values = global_times_colors_B
    ) +
    scale_fill_manual(
      values = global_times_colors_B
    ) +
    coord_cartesian(ylim = c(min_range, max_range))+
    labs(
      x = "time increments [5 min]",
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
  
  
  comparisonplot_bB<-ggplot(summary_bB,
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
      values = global_times_colors_bB
    ) +
    scale_fill_manual(
      values = global_times_colors_bB
    ) +
    coord_cartesian(ylim = c(min_range, max_range))+
    labs(
      x = "time increments [5 min]",
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
  
  
  comparisonplot_bA<-ggplot(summary_bA,
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
      values = global_times_colors_bA
    ) +
    scale_fill_manual(
      values = global_times_colors_bA
    ) +
    coord_cartesian(ylim = c(min_range, max_range))+
    labs(
      x = "time increments [5 min]",
      y = param,
      color = "extracted from:",
      fill  = "extracted from:"
    ) +
    theme_minimal()+
    theme(
      legend.position = c(0.98, 0.98),     # inside top-right
      legend.justification = c(1, 1),      # anchor legend at its top-right
      legend.text  = element_text(size = 8),
      legend.key.size = unit(0.6, "cm")    # size of legend keys
    )
  
  #combine into one plot with same legend
  
  # grid_plot<-(comparisonplot_b| comparisonplot_B| comparisonplot_bB|comparisonplot_bA) +
  #   plot_layout(ncol=2,nrow=2, guides = "collect", axes = "collect_y" ) &
  #   theme(legend.position = "bottom")
  # 
  grid_plot<-(comparisonplot_b| comparisonplot_B| comparisonplot_bB|comparisonplot_bA) +
    plot_layout( axes = "collect" ) 
   
  
  ggsave(
    filename = paste0(folder_path_timeplots, "/comparison/", param,"_","overlaybeforeafterinf_ALL.png"),
    plot = grid_plot,
    width = 30,
    height = 5,
    dpi = 300
  )
  
  
  
}

#########
#plot only infection period overlay for all 4 treatments
#########
folder_path_timeplots<-"network_parameter_plots/exp2_inf/5mins/time_mean/"
infection_data<-before_after_infection_merged%>%filter(infection == "after")

for (param in parameters){

  summary_df <- infection_data |>
    group_by(treatment, time_interval) |>
    summarise(
      mean_val = mean(.data[[param]], na.rm = TRUE),
      sd_val   = sd(.data[[param]], na.rm = TRUE),
      n             = n(),
      se_val   = sd_val / sqrt(n),
      .groups = "drop"
    )
  
  overlay_inf_plot<-ggplot(summary_df,
                       aes(x = time_interval,
                           y = mean_val,
                           color = treatment,
                           fill  = treatment
                          )) +
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
      x = "time increments [5min]",
      y = param,
      color = "mean over 8 \n colony replicates",
     fill= "sd" ) +
    theme_minimal()
  
  ggsave(
    filename = paste0(folder_path_timeplots, param,"_","infection_overlay.png"),
    plot = overlay_inf_plot,
    width = 10,
    height = 5,
    dpi = 300
  )
  
}

#boxplots of values

for (param in parameters){
  

 
  boxplot_inf_plot<-ggplot(infection_data,
                           aes(x = treatment,
                               y = .data[[param]],
                               color = colony_name,
                               fill  = treatment,
                               group = treatment
                           )) +
    geom_boxplot(alpha= 0.6)+
   geom_jitter()+
    scale_color_manual(
      values = colony_colors
    ) +
    scale_fill_manual(
      values = treatment_colors
    ) +
    labs(
      x = "treatments",
      y = param,
      color = "15 x 5 min per \n colony",
      fill= "treatment"
     ) +
    theme_minimal()
  
  ggsave(
    filename = paste0("network_parameter_plots/exp2_inf/5mins/time_mean/boxplot/", param,"_","infection_boxplot.png"),
    plot = boxplot_inf_plot,
    width = 10,
    height = 5,
    dpi = 300
  )
  
}

#boxplots of means
infection_mean<-infection_data%>%group_by(colony_name, treatment)%>%summarise(
  across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
  .groups = "drop"
)

for (param in parameters){
  png(file= paste0("network_parameter_plots/exp2_inf/5mins/time_mean/hist_boxplot_means/", param, ".png"),
      width=600, height=350)
  hist(infection_mean[[param]], xlab = param , col = "purple", main="")
  dev.off()
}

# plot boxplot of means
for (param in parameters){

boxplot_inf_mean_plot<-ggplot(infection_mean,
                         aes(x = treatment,
                             y = .data[[param]],
                             color = colony_name,
                             fill  = treatment,
                             group = treatment
                         )) +
  geom_boxplot(alpha= 0.6)+
  geom_jitter()+
  scale_color_manual(
    values = colony_colors
  ) +
  scale_fill_manual(
    values = treatment_colors
  ) +
  labs(
    x = "treatments",
    y = param,
    color = "1 x 5 min mean per \n colony",
    fill= "treatment"
  ) +
  theme_minimal()

ggsave(
  filename = paste0("network_parameter_plots/exp2_inf/5mins/time_mean/boxplots_means/", param,"_","infection_boxplot.png"),
  plot = boxplot_inf_mean_plot,
  width = 10,
  height = 5,
  dpi = 300
)

}

#############################
#glmm
#############################
hist(infection_mean$strength_mean)
model_g<- glmmTMB(
  log(strength_mean) ~ treatment,
  data = infection_mean, 
  #family = Gamma(link = log)# 'ensures is positive', right skewed (link = log)
  family = gaussian()
)

sim <- DHARMa::simulateResiduals(model_g)
plot(sim)

drop1(model_g, test = "Chisq")

emmeans(model_g, pairwise ~ treatment )
##density_aggr
hist(infection_mean$density_aggr)
model_d<- glmmTMB(
  log(density_aggr) ~ treatment,
  data = infection_mean, 
  #family = Gamma(link = log)# 'ensures is positive', right skewed (link = log)
  family = gaussian()
)

sim <- DHARMa::simulateResiduals(model_d)
plot(sim)

drop1(model_d, test = "Chisq")

emmeans(model_d, pairwise ~ treatment )

#interactionlength_mean
hist(infection_mean$interaction_length_mean)
model_lm<- glmmTMB(
  log(interaction_length_mean) ~ treatment,
  data = infection_mean, 
  #family = Gamma(link = log)# 'ensures is positive', right skewed (link = log)
  family = gaussian()
)

sim <- DHARMa::simulateResiduals(model_lm)
plot(sim)

drop1(model_lm, test = "Chisq")

emmeans(model_lm, pairwise ~ treatment )

#interaction length sd
hist(infection_mean$interaction_length_sd)
model_lsd<- glmmTMB(
  interaction_length_sd ~ treatment,
  data = infection_mean, 
  #family = Gamma(link = log)# 'ensures is positive', right skewed (link = log)
  family = gaussian()
)

sim <- DHARMa::simulateResiduals(model_lsd)
plot(sim)

drop1(model_lsd, test = "Chisq")

emmeans(model_lsd, pairwise ~ treatment )

#global_eff
hist(infection_mean$global_eff)
model_ge<- glmmTMB(
  global_eff ~ treatment,
  data = infection_mean, 
  #family = Gamma(link = log)# 'ensures is positive', right skewed (link = log)
  family = gaussian()
)

sim <- DHARMa::simulateResiduals(model_ge)
plot(sim)

drop1(model_ge, test = "Chisq")

#clustering coeff global
hist(infection_mean$clustering_global)
model_cg<- glmmTMB(
  clustering_global ~ treatment,
  data = infection_mean, 
  #family = Gamma(link = log)# 'ensures is positive', right skewed (link = log)
  family = gaussian()
)
sim <- DHARMa::simulateResiduals(model_cg)
plot(sim)

drop1(model_cg, test = "Chisq")

#mean distance

hist(infection_mean$mean_distance)
model_md<- glmmTMB(
  mean_distance ~ treatment,
  data = infection_mean, 
  dispformula = ~1,
  #family = Gamma(link = log)# 'ensures is positive', right skewed (link = log)
  family = Gamma(link = "log")
)
sim <- DHARMa::simulateResiduals(model_md)
plot(sim)

drop1(model_md, test = "Chisq")

#waiting_time_mean
hist(infection_mean$waiting_time_mean)
model_wtm<- glmmTMB(
  log(waiting_time_mean) ~ treatment,
  data = infection_mean, 
  #family = Gamma(link = log)# 'ensures is positive', right skewed (link = log)
  family = gaussian()
)

model_wtm2<- glmmTMB(
  waiting_time_mean ~ treatment,
  data = infection_mean, 
  #family = Gamma(link = log)# 'ensures is positive', right skewed (link = log)
  family = Gamma(link = "log")
)
sim <- DHARMa::simulateResiduals(model_wtm2)
plot(sim)
 AIC(model_wtm, model_wtm2)
drop1(model_wtm, test = "Chisq")

# interaction sum
hist(infection_mean$total_sum_of_interactions)
model_si<- glmmTMB(
  total_sum_of_interactions ~ treatment,
  data = infection_mean, 
  #family = Gamma(link = log)# 'ensures is positive', right skewed (link = log)
  family = gaussian()
)
sim <- DHARMa::simulateResiduals(model_si)
plot(sim)
drop1(model_si, test = "Chisq")

#total number interactions
hist(infection_mean$total_number_interactions)
model_ni<- glmmTMB(
  total_number_interactions ~ treatment,
  data = infection_mean, 
  #family = Gamma(link = log)# 'ensures is positive', right skewed (link = log)
  family = gaussian()
)
sim <- DHARMa::simulateResiduals(model_ni)
plot(sim)
drop1(model_ni, test = "Chisq")

##########
#comparison between treatments
##########
no_infection_mean<-before_after_infection_merged%>%filter(infection != "after")%>%group_by(colony_name, treatment, globaltime)%>%summarise(
  across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
  .groups = "drop"
)

for (param in parameters){
  png(file= paste0("network_parameter_plots/exp2_12hs/5mins/time_continuum/overlays/comparison/hist/", param, ".png"),
      width=600, height=350)
  hist(no_infection_mean[[param]], xlab = param , col = "purple", main="")
  dev.off()
}

#strength_mean
hist(no_infection_mean$strength_mean)
model_g<- glmmTMB(
  strength_mean ~ treatment + (1|colony_name),
  data = no_infection_mean, 
  dispformula = ~treatment,
  #family = Gamma(link = log)# 'ensures is positive', right skewed (link = log)
  family = gaussian()
)

sim <- DHARMa::simulateResiduals(model_g)
plot(sim)

drop1(model_g, test = "Chisq")

emmeans(model_g, pairwise ~ treatment )

#density aggr
hist(no_infection_mean$density_aggr)
model_d<- glmmTMB(
  density_aggr ~ treatment + (1|colony_name),
  data = no_infection_mean, 
  dispformula = ~treatment,
  #family = Gamma(link = log)# 'ensures is positive', right skewed (link = log)
  family = gaussian()
)

sim <- DHARMa::simulateResiduals(model_d)
plot(sim)

drop1(model_d, test = "Chisq")

emmeans(model_g, pairwise ~ treatment )

#global_effectiveness

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
#ANALYSING pre and post infection differences on colony level

############################################
#mean over all time intervalls for all colonies
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

prep_delta <- mean_overtimes %>%
  pivot_longer(
    cols = -c(colony_name, infection, globaltime),
    names_to = "parameter",
    values_to = "mean_value")
  
prep_delta_1<-prep_delta%>%filter(globaltime=="until 1.25h")
prep_delta_2<-prep_delta%>%filter(globaltime=="1.25-2.5h")
prep_delta_3<-prep_delta%>%filter(globaltime=="10h-11.25h")
prep_delta_4<-prep_delta%>%filter(globaltime=="8.75h-10h")

prep_delta_inf<-prep_delta%>%filter(globaltime=="1.25h after pathogen exposure")

prep_delta_1<-prep_delta_1%>%mutate(delta = mean_value- prep_delta_inf$mean_value)
prep_delta_1<-prep_delta_1%>%mutate(delta_norm = (mean_value- prep_delta_inf$mean_value)/mean_value)

prep_delta_2<-prep_delta_2%>%mutate(delta = mean_value- prep_delta_inf$mean_value)
prep_delta_2<-prep_delta_2%>%mutate(delta_norm = (mean_value- prep_delta_inf$mean_value)/mean_value)

prep_delta_3<-prep_delta_3%>%mutate(delta = mean_value- prep_delta_inf$mean_value)
prep_delta_3<-prep_delta_3%>%mutate(delta_norm = (mean_value- prep_delta_inf$mean_value)/mean_value)

prep_delta_4<-prep_delta_4%>%mutate(delta = mean_value- prep_delta_inf$mean_value)
prep_delta_4<-prep_delta_4%>%mutate(delta_norm = (mean_value- prep_delta_inf$mean_value)/mean_value)

#
prep_delta_all<-bind_rows(prep_delta_1, prep_delta_2, prep_delta_3, prep_delta_4)

prep_delta_all<-prep_delta_all%>%mutate(treatment = sub("\\d+.*$", "", colony_name))

delta_strength<-prep_delta_all%>%filter(parameter=="strength_mean"
  )


hist(delta_strength$delta_norm)

##glmm for strength
###################################
model_delta_s<-glmmTMB(delta_norm~treatment + (1|colony_name), 
                       data = delta_strength,
                       dispformula = ~1,
                       family = gaussian())

sim <- DHARMa::simulateResiduals(model_delta_s)
plot(sim)

pairs <- emmeans(model_delta_s, pairwise ~ treatment, type = "response")
stat_test <- pairs$contrasts
stat_test

########
model_delta_s<-glmmTMB(delta~treatment + (1|colony_name), 
                       data = delta_strength,
                       dispformula = ~1,
                       family = gaussian())

sim <- DHARMa::simulateResiduals(model_delta_s)
plot(sim)

pairs <- emmeans(model_delta_s, pairwise ~ treatment, type = "response")
stat_test <- pairs$contrasts
stat_test

#for density
hist(delta_density$delta)
delta_density<-delta%>%filter(parameter=="density_aggr")
model_delta_d<-glmmTMB(delta_norm~treatment, 
                       data = delta_density,
                       family = gaussian())

sim <- DHARMa::simulateResiduals(model_delta_d)
plot(sim)

pairs <- emmeans(model_delta_d, pairwise ~ treatment, type = "response")
stat_test <- pairs$contrasts
stat_test