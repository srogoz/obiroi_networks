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
setwd("Desktop/EXP1_analysis_BROKEN")
source("function_collection.R")



#######

expected_ants= c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP")

selected_colonies<-list(
  "a16-1", "a16-2", "a16-3", "a16-4","a16-5","a16-6","a16-7", "a16-8",
  "b16-1", "b16-2", "b16-3", "b16-4", "b16-5", "b16-6", "b16-7","b16-8",
  "ba16-1","ba16-2", "ba16-3", "ba16-4", "ba16-5",  "ba16-6","ba16-7","ba16-8"
)
adjmatrix_list<-get_colony_files(source_folder_adjmatrix, selected_colonies)

present_ants_list<-list( "a16-1" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "a16-2" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "a16-3" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "a16-4" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "a16-5" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "a16-6" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "a16-7" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "a16-8" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "b16-1" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "b16-2" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "b16-3" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "b16-4" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "b16-5" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO"),
                         "b16-6" = c("BB", "BG", "BO", "BP", "GB", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "b16-7" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "b16-8" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "ba16-1" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OO", "PB", "PG", "PO", "PP"),
                         "ba16-2" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "ba16-3" = c("BB", "BG", "BO", "BP", "GB", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "ba16-4" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "ba16-5" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "ba16-6" = c("BB", "BG", "BO", "BP", "GB", "GO", "OB", "OG", "OP", "PB", "PG", "PO", "PP"),
                         "ba16-7" = c("BB", "BG", "GB", "GG", "GO", "GP", "OB", "OG", "OP", "PB", "PG", "PO", "PP"),
                         "ba16-8" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP")
)

framerate<-list("a16-1" = 10,
                "a16-2" = 10,
                "a16-3" = 10,
                "a16-4" = 10,
                "a16-5" = 10,
                "a16-6" = 10,
                "b16-1" = 10,
                "b16-2" = 10,
                "b16-3"= 10,
                "b16-4" = 10,
                "b16-5" = 10,
                "b16-6" = 10,
                "ba16-1"= 10,
                "ba16-2"=10,
                "ba16-3"=10,
                "ba16-4"=10, 
                "ba16-5"=10,
                "ba16-6"=10,
                "a16-7"=5,
                "a16-8"=5,
                "b16-7"=5,
                "b16-8"=5,
                "ba16-7"=5,
                "ba16-8"=5)


treatment_match<-list("a16-1" = "1",
                      "a16-2" = "1",
                      "a16-3" = "1",
                      "a16-4" = "1",
                      "a16-5" = "1",
                      "a16-6" = "1",
                      "a16-7" = "1",
                      "a16-8" = "1",
                      "b16-1" = "2", 
                      "b16-2" = "2",
                      "b16-3" = "2", 
                      "b16-4" = "2", 
                      "b16-5" = "2", 
                      "b16-6" = "2",
                      "b16-7" = "2",
                      "b16-8" = "2",
                      "ba16-1" = "3", 
                      "ba16-2" = "3",
                      "ba16-3" = "3",
                      "ba16-4" = "3",
                      "ba16-5" = "3",
                      "ba16-6" = "3",
                      "ba16-7" = "3",
                      "ba16-8" = "3"
)
genotype_match<-list("BB" = "b",
                     "BG" = "b",
                     "BO" = "b",
                     "BP" = "b",
                     "GB" = "b",
                     "GG" = "b", 
                     "GO" = "b",
                     "GP" = "b", 
                     "OB" = "a",
                     "OG" = "a",
                     "OO" = "a",
                     "OP" = "a",
                     "PB" = "a",
                     "PG" = "a",
                     "PO" = "a",
                     "PP" = "a")

parameters <-c("total_number_interactions", "total_sum_of_interactions", "strength_mean", "strength_sd", "interaction_length_mean", "interaction_length_sd", "waiting_time_mean", "waiting_time_sd", "burstiness", "density_aggr", "global_eff", "assortativity", "mean_distance")

#############################
#CREATE EDGELISTS
#############################
#check number of frames

for (i in 1:length(selected_colonies)){
  network<-np$load(adjmatrix_list[[selected_colonies[[i]]]])
  print(dim(network)[3])
}
########################
#set truncation time for infected 75mins



interval_time<-5
number_intervalls<-15


####################

#unpper and lower interaction limit

interaction_limit_min<-2
interaction_limit_sec<-1

source_folder_edgelist<-"edge_lists/exp1_inf/75mins/"
folder_path_networkpara <-"network_parameters/exp1_inf/5mins/"
output_folder <- "network_parameter_plots/exp1_inf/5mins/"
timestamp <- format(Sys.time(), "%d%m%Y")

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
  
  
  interval_size<-interval_time*60*framerate_col
  intervall_timepoints <- seq(
    from = 0,
    by   = interval_size,
    length.out = number_intervalls + 1
  )
  
  for (j in 1:number_intervalls){
    #start and stop
    time_window_start<-intervall_timepoints[j]
    time_window_end<-intervall_timepoints[j+1]
    
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
                                        burstiness = burstiness
                                        
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
results_networkparameters_allintervalls<-read.csv(file = "network_parameters/exp1_inf/5mins/ba16-8_network_parameters_allintervalls_10022026.csv")
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
time development per treatment and per colony_colors
treatments<-c("a", "b", "ba")
folder_path_timeplots<-"network_parameter_plots/exp1_inf/5_mins/time_continuum/"
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
#antlevel data
###################
#plot infected individual in another color



