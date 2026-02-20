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

#EXP1_12HS
#list of adj-matrix files
source_folder_adjmatrix<-"adj_matrix/exp1_inf/"
folder_edgelist<-"edge_lists/exp1_inf/"

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

expected_ants= c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP")
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

truncationtime_5f<-22500
truncationtime_10f<-45000


#30 mins
truncationtime_10f_30<-18000
truncationtime_5f_30<-9000


for (i in 1:length(selected_colonies)){
  
  colony_name<-selected_colonies[[i]]
  network<-np$load(adjmatrix_list[[selected_colonies[[i]]]])
  
  colnames(network)<-expected_ants
  row.names(network)<-expected_ants
  ###########
  #create edgelist from interaction array
  # transform into network object, undirected
  #reduce  all networks to 12hs: 432000 frames, starting from the end, account for difference in framerate
  
  if (framerate[[selected_colonies[[i]]]]==10){
    network_cut<-network[,,1:truncationtime_10f]}
  else
  { network_cut<-network[,,1:truncationtime_5f]}
  
  
  network_obj <-make_edgelist(network_cut, expected_ants)
  
  #save edgelist
  timetag<-format(Sys.time(), "%d%m%Y")
  saveRDS(network_obj,
          paste0(folder_edgelist ,colony_name, "_edgelist_", timetag, ".rds"))
  
  print(paste0(colony_name, " done"))
  
}


##################
#PARAMETERS
##################
#functional
interaction_limit_min<-2
interaction_limit_sec<-1
#########
#organisational
source_folder_edgelist<-"edge_lists/exp1_inf/"
folder_path_networkpara <-"network_parameters/exp1_inf/"
output_folder <- "network_parameter_plots/exp1_inf/"
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

for (i in 1:length(selected_colonies)){
  
  colony_name<-selected_colonies[[i]]
  framerate_col<-framerate[[colony_name]]
  present_ants<-present_ants_list[[colony_name]]
  lower_limit_interaction<-interaction_limit_sec*framerate_col
  higher_limit_interaction<-interaction_limit_min*60*framerate_col
  
  ###################
  
  #prepare network obj, interaction limits
  network_obj <-readRDS(network_obj_list[[colony_name]])
  
  ##############
  #CONDITIONS
  ##############
  #truncate from both sides
  #network_obj_5 <-network_obj[(network_obj$duration>lower_limit_interaction & network_obj$duration<higher_limit_interaction),]
  #truncate and minimize from the high side. Still count interaction but does not keep counting after time
  # Cap durations at the upper limit
  network_obj_5 <-network_obj[(network_obj$duration>lower_limit_interaction),]
  network_obj_5$duration <- pmin(network_obj_5$duration, higher_limit_interaction)
  
  #remove (empty) data of missing ants
  missing_ants <- setdiff(expected_ants, present_ants)
  
  network_obj_5 <- network_obj_5[
    !network_obj_5$head %in% missing_ants &
      !network_obj_5$tail %in% missing_ants,
  ]
  
  
  
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
  waiting_time_ant_mean_df[[colony_name]]<-interactionlength_means_vec
  waiting_time_ant_sd_df[[colony_name]]<-waiting_time_ant_sd_vec
  interactionlength_ant_mean[[colony_name]]<-interactionlength_means_vec
  burstiness_ant_df[[colony_name]]<-burstiness_ant_vec
  
  #initialize columnnames with numbers for networks
  network_obj_5$tail<-match(network_obj_5$tail, expected_ants)
  network_obj_5$head<-match(network_obj_5$head, expected_ants)
  
  #cretae aggregated network with matching interaction-time filters
  aggregated_network<-aggregate_from_edgelist(network_obj_5)
  
  
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
  duration<-12*60*60*framerate_col
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
  strength_antlevel_mean[[colony_name]]<-strength_full
  
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
  local_efficiency[[colony_name]]<-local_eff_full
  
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
  centrality_antlevel_df[[colony_name]]<-centrality_full
  
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
  results_networkparameters<-data.frame(colony_name = colony_name,
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
  write.csv(results_networkparameters,
            file = paste0(folder_path_networkpara, colony_name, "_network_parameters_", timestamp, ".csv"),
            row.names = FALSE
  )
  #progress check
  print(paste0(colony_name, " done"))
  
}

##########################
#save antlevel data
write.csv(strength_antlevel_mean,
          file = paste0(folder_path_networkpara, "ANTLEVEL_strength_mean", timestamp, ".csv"),
          row.names = FALSE)
write.csv(local_efficiency,
          file = paste0(folder_path_networkpara, "ANTLEVEL_local_efficiency", timestamp, ".csv"),
          row.names = FALSE)
write.csv(waiting_time_ant_mean_df,
          file = paste0(folder_path_networkpara, "ANTLEVEL_waitingtime_mean", timestamp, ".csv"),
          row.names = FALSE)
write.csv(waiting_time_ant_sd_df,
          file = paste0(folder_path_networkpara, "ANTLEVEL_waitingtime_sd", timestamp, ".csv"),
          row.names = FALSE)
write.csv(interactionlength_ant_mean,
          file = paste0(folder_path_networkpara, "ANTLEVEL_interactionlength_mean", timestamp, ".csv"),
          row.names = FALSE)
write.csv(burstiness_ant_df,
          file = paste0(folder_path_networkpara, "ANTLEVEL_burstiness", timestamp, ".csv"),
          row.names = FALSE)
write.csv(centrality_antlevel_df,
          file = paste0(folder_path_networkpara, "ANTLEVEL_centrality", timestamp, ".csv"),
          row.names = FALSE)
##########################

##########################
##%%read network parameters
source_folder_networkpara <-"network_parameters/exp1_inf/"
folder_path_networkplots <-"network_parameter_plots/exp1_inf/"
folder_path_qqplots<-"network_parameter_plots/exp1_inf/STATS/qqplots/"

treatment_colors<-c(
  "a" = "#05e0fc",
  "b" = "#fc0536",
  "ba" = "#db05fc"
)

# overview parameters
parameters <-c("total_number_interactions", "total_sum_of_interactions", "strength_mean", "strength_sd", "interaction_length_mean", "interaction_length_sd", "waiting_time_mean", "waiting_time_sd", "burstiness", "density_aggr", "global_eff", "assortativity", "mean_distance")
parameters_units <-list("total_number_interactions" = "",
                        "total_sum_of_interactions" = "[s]",
                        "strength_mean" = "[s]",
                        "strength_sd" = "[s]",
                        "interaction_length_mean" = "[s]",
                        "interaction_length_sd"= "[s]",
                        "waiting_time_mean" = "[s]",
                        "waiting_time_sd" = "[s]",
                        "burstiness" = "",
                        "density_aggr" = "",
                        "global_eff" = "",
                        "assortativity" = "",
                        "mean_distance" ="")


network_parameter_list<-get_colony_files(source_folder_networkpara, selected_colonies)
network_parameter_collected<-list()

for (i in 1:length(selected_colonies)){
  network_parameter_col<-read.csv(network_parameter_list[[selected_colonies[[i]]]])
  network_parameter_col<-network_parameter_col[1,]
  network_parameter_collected<-rbind(network_parameter_collected, network_parameter_col)
}
###adjust treatment and replicates
#treatments a,b,ba

network_parameter_collected_with_treatment<-data.frame(treatment = rep(c("a", "b", "ba"), each = 8), replicate = rep(1:8, times = 3))

network_parameter_complete<-data_frame(network_parameter_collected_with_treatment, network_parameter_collected)
#treatment only b, a 
network_parameter_ab<-network_parameter_complete %>%
  filter(treatment != "ba")
########################
##TEST FOR NORMAL DISTRIBUTION
##null hypothesis to Shapiro Wilks:
# sample distribution is normal: if p significant p distribution is not normal

distribution_test <- list()

for (para in parameters) {
  # Compute Shapiro–Wilk test per treatment for this parameter
  res <- network_parameter_complete %>%
    group_by(treatment) %>%
    summarise(
      shapiro_p = shapiro.test(.data[[para]])$p.value,
      .groups = "drop"
    )
  
  # Store results in list using parameter name
  distribution_test[[para]] <- res
}

# Combine all results into one data frame
distribution_test_df <- bind_rows(
  lapply(names(distribution_test), function(p) {
    data.frame(
      parameter = p,
      treatment = distribution_test[[p]]$treatment,
      shapiro_p = distribution_test[[p]]$shapiro_p
    )
  })
)



#### not normally distributed: KRUSKAL_WALLIS
kruskal_test_results <- list()
dunn_tests <- list()

for (para in parameters) {
  formula <- as.formula(paste(para, "~ treatment"))
  
  res <- kruskal.test(formula, data = network_parameter_complete)
  p_val <- res$p.value
  
  kruskal_test_results[[para]] <- p_val
  
  # Dunn only if KW significant
  if (p_val < 0.05) {
    dunn_tests[[para]] <- network_parameter_complete %>%
      dunn_test(as.formula(paste(para, "~ treatment")),
                p.adjust.method = "BH")
  }
}

folder_colonylevel_qqplots<-"network_parameter_plots/exp1_inf/STATS/qqplots/"
anova_pvals <- list()
posthoc_tests <- list()


for (para in parameters) {
  formula <- as.formula(paste(para, "~ treatment"))
  
  model <- aov(formula, data = network_parameter_complete)
  
  # Save diagnostic plots
  png(paste0(folder_colonylevel_qqplots, para, "_anova_diagnostics.png"),
      width = 800, height = 400)
  par(mfrow = c(1, 2))
  plot(model, which = 1)
  plot(model, which = 2)
  dev.off()
  
  # ANOVA p-value
  p_anova <- summary(model)[[1]][["Pr(>F)"]][1]
  anova_pvals[[para]] <- p_anova
  
  # Tukey only if significant
  if (p_anova < 0.05) {
    posthoc_tests[[para]] <- TukeyHSD(model)
  }
}

###########
#plot all violin plots

comparisons <- list(
  c("b", "a"),
  c("ba", "a"),
  c("ba", "b")
)



for (param in parameters){
  p_val<-round(kruskal_test_results[[param]],2)
  #p_val<-round(kruskal_test_results[[param]]$p.value,2)
  #p_val<-round(anova_test_results[[param]],2)
  param_plot<-ggplot(network_parameter_complete, 
                     aes(x = treatment, y = .data[[param]], fill = treatment)) +
    geom_violin(trim = FALSE, alpha = 0.4) +  # “bell shape”
    geom_boxplot(width = 0.15, outlier.shape = NA, color = "black") +
    geom_jitter(width = 0.1, alpha = 0.7) +   # individual points
    geom_text(aes(label = replicate), vjust = -1, size = 3.5, color = "black") +
    theme_minimal() +
     stat_compare_means(
      comparisons = comparisons,
      method = "wilcox.test",   # or "t.test"
      label = "p.signif",
      hide.ns = TRUE
    )+
    # Apply manual colors
    scale_fill_manual(values = treatment_colors) +
    labs(
      title = paste0(param, " p = ", p_val ),
      x = "Treatment",
      y = param) 
  #  scale_fill_viridis_d(option = "D")
  
  # Save each plot as PNG
  ggsave(
    filename = paste0(folder_path_networkplots, "/", param, "_violin_boxplot.png"),
    plot = param_plot,
    width = 6,
    height = 5,
    dpi = 300
  )
}



########################
########################
#GLMM ON ANTLEVEL
########################




antlevel_parameters<-c("strength_mean", "burstiness", "waitingtime_mean", "centrality", "interactionlength_mean", "local_efficiency", "waitingtime_sd")

png(file="network_parameter_plots/exp1_12hs/hist/",
    width=600, height=350)
hist(df_glmm$strength, xlab = "strength [hs]", col = "purple", main="")
dev.off()


########################################
##df wide preparation for all parameters
########################################
###
read_param_wide <- function(param) {
  
  df <- read.csv(
    paste0(
      "network_parameters/exp1_inf/30mins/ANTLEVEL_",
      param,
      "06022026.csv"
    )
  )
  
  df$ant <- expected_ants
  
  df_long <- df |>
    tidyr::pivot_longer(
      cols = -ant,
      names_to = "colony",
      values_to = param
    )
  
  df_long
}




df_list <- lapply(antlevel_parameters, read_param_wide)

df_wide <- Reduce(
  function(x, y) merge(x, y, by = c("ant", "colony")),
  df_list
)

df_wide <- df_wide |>
  dplyr::mutate(
    treatment = sub("\\d+.*$", "", colony),
    genotype = dplyr::case_when(
      treatment == "ba" ~ unlist(genotype_match[ant]),
      treatment == "a"  ~ "a",
      treatment == "b"  ~ "b"
    ),
    puremixed = ifelse(treatment == "ba", "m", "p")
  )

#plot histograms of the colonylevel parameters

for (param in antlevel_parameters){
  png(file=paste0("network_parameter_plots/exp1_inf/antlevel/hist/",param, ".png"),
      width=600, height=350)
  hist(df_wide[[param]], xlab = param, col = "purple", main= "")
  dev.off()
}

##################
#CENTRALITY
##################
#correct centrality for 1
eps <- 1e-6
df_wide$centrality_beta <- pmin(
  pmax(df_wide$centrality, eps),
  1 - eps )
 

model_c<- glmmTMB(
  centrality_beta ~ treatment + (1|colony),
  data = df_wide,
  #family = Gamma(link = log)# 'ensures is positive', right skewed (link = log)
  family = beta_family()
)
sim <- DHARMa::simulateResiduals(model_c)
plot(sim)
AIC(model_c)


model_c1 <- glmmTMB(
  centrality_beta ~ treatment + (1 | colony/ant) ,
  data = df_wide,
  #family = Gamma(link = log)# 'ensures is positive', right skewed (link = log)
  family = beta_family()
)
sim <- DHARMa::simulateResiduals(model_c1)
plot(sim)
AIC(model_c1)

model_c2 <- glmmTMB(
  centrality_beta ~ genotype*puremixed + (1 | colony/ant),
  data = df_wide,
  
  family = beta_family(link="probit")
  #family = beta_family(link = "cloglog")
)

sim <- DHARMa::simulateResiduals(model_c2)
plot(sim)
AIC(model_c2)

model_gt <- glmmTMB(
  centrality_beta ~ genotype*puremixed + (1 | colony/ant),
  data = df_wide,
  family = beta_family(link = "probit")
)

sim <- DHARMa::simulateResiduals(model_gt)
plot(sim)

model_gt2 <- glmmTMB(
  centrality_beta ~ genotype+puremixed + (1 | colony/ant),
  data = df_wide,
  #family = Gamma(link = log)# 'ensures is positive', right skewed (link = log)
  family = beta_family()
)

sim <- DHARMa::simulateResiduals(model_gt2)
plot(sim)

model_gt3 <- glmmTMB(
  centrality_beta ~ genotype + (1 | colony/ant),
  data = df_wide,
  #family = Gamma(link = log)# 'ensures is positive', right skewed (link = log)
  family = beta_family(link ="probit")
)

sim <- DHARMa::simulateResiduals(model_gt3)
plot(sim)

AIC(model_gt3)


model_gt4 <- glmmTMB(
  centrality_beta ~ puremixed + (1 | colony/ant),
  data = df_wide, 
  #family = Gamma(link = log)# 'ensures is positive', right skewed (link = log)
  family = beta_family()
)

sim <- DHARMa::simulateResiduals(model_gt4)
plot(sim)


AIC(model_gt4, model_gt3, model_gt2, model_gt)
################################
###################################
emm_trt <- emmeans(model_gt, ~ genotype*puremixed)
pairs(emm_trt, adjust = "tukey")

emmeans(model_gt, pairwise ~ genotype | puremixed)
emmeans(model_gt, pairwise ~ puremixed | genotype)




pairs <- emmeans(model_gt2, pairwise ~ genotype*puremixed, type = "response")
stat_test <- pairs$contrasts

############
#strength
############
model_s <- glmmTMB(
  strength_mean ~ genotype*puremixed + (1|colony/ant),
  data = df_wide,
  family = gaussian()# 'ensures is positive', right skewed (link = log)
  
)

sim <- DHARMa::simulateResiduals(model_s)
plot(sim)

emm_trt <- emmeans(model_s, ~ treatment)
pairs(emm_trt, adjust = "tukey")


model_s1 <- glmmTMB(
  strength_mean ~ genotype*puremixed + (1|colony),
  data = df_wide,
  family = gaussian()# 'ensures is positive', right skewed (link = log)
  
)

sim <- DHARMa::simulateResiduals(model_s1)
plot(sim)

model_s2 <- glmmTMB(
  strength_mean ~ genotype + (1|colony/ant),
  data = df_wide,
  family = gaussian()# 'ensures is positive', right skewed (link = log)
  
)

sim <- DHARMa::simulateResiduals(model_s2)
plot(sim)




model_s2 <- glmmTMB(
  strength_mean ~ genotype*puremixed + (1|colony/ant),
  data = df_wide,
  family = Gamma(link = "log")# 'ensures is positive', right skewed (link = log)
  
)

sim <- DHARMa::simulateResiduals(model_s2)
plot(sim)


model_s3 <- glmmTMB(
  strength_mean ~ treatment + (1|colony/ant),
  data = df_wide,
  family = Gamma(link = "log")# 'ensures is positive', right skewed (link = log)
  
)

sim <- DHARMa::simulateResiduals(model_s3)
plot(sim)


emm_trt <- emmeans(model_s3, ~ treatment)
pairs(emm_trt, adjust = "tukey")

############
#burstiness

model_b <- glmmTMB(
  burstiness ~ treatment+ (1|colony),
  data = df_wide,
  family = gaussian()# 'ensures is positive', right skewed (link = log)
  
)

sim <- DHARMa::simulateResiduals(model_b)
plot(sim)

AIC(model_b)

pairs <- emmeans(model_b2, pairwise ~ genotype*puremixed, type = "response")
stat_test <- pairs$contrasts


model_b2 <- glmmTMB(
  burstiness ~ genotype*puremixed + (1|colony/ant),
  data = df_wide,
  family = gaussian()# 'ensures is positive', right skewed (link = log)
  
)

sim <- DHARMa::simulateResiduals(model_b2)
plot(sim)

model_b3<- glmmTMB(
  burstiness ~ genotype + (1|colony/ant),
  data = df_wide,
  family = gaussian()# 'ensures is positive', right skewed (link = log)
  
)

sim <- DHARMa::simulateResiduals(model_b3)
plot(sim)

BIC(model_b3)

model_b4<- glmmTMB(
  burstiness ~ genotype + (1|colony/ant),
  data = df_wide,
  family = gaussian()# 'ensures is positive', right skewed (link = log)
  
)

sim <- DHARMa::simulateResiduals(model_b4)
plot(sim)

BIC(model_b4)


###############
#BURSTINESS
###############

model <- glmmTMB(
  strength_mean ~ genotype*puremixed + (1|colony/ant),
  data = df_wide,
  #family = Gamma(link = log)# 'ensures is positive', right skewed (link = log)
  family = gaussian()
)

############################
#WAITINGTIME
############################

model_w<- glmmTMB(
  waitingtime_mean ~ genotype + (1|colony),
  data = df_wide,
  family = Gamma(link = log)# 'ensures is positive', right skewed (link = log)
  
)

sim <- DHARMa::simulateResiduals(model_w)
plot(sim)

BIC(model_w)

##################################
model_w<- glmmTMB(
  waitingtime_mean ~ genotype*puremixed + (1|colony/ant),
  data = df_wide,
  family = Gamma(link = log)# 'ensures is positive', right skewed (link = log)
  
)

sim <- DHARMa::simulateResiduals(model_w)
plot(sim)

BIC(model_w)

pairs <- emmeans(model_w, pairwise ~ genotype*puremixed, type = "response")
stat_test <- pairs$contrasts
#############
model_w<- glmmTMB(
  waitingtime_mean ~ genotype*puremixed + (1|colony),
  data = df_wide,
  family = Gamma(link = log)# 'ensures is positive', right skewed (link = log)
  
)

sim <- DHARMa::simulateResiduals(model_w)
plot(sim)

BIC(model_w)
##############
model_w<- glmmTMB(
  waitingtime_mean ~ genotype+puremixed + (1|colony/ant),
  data = df_wide,
  family = Gamma(link = log)# 'ensures is positive', right skewed (link = log)
  
)

sim <- DHARMa::simulateResiduals(model_w)
plot(sim)

BIC(model_w)

#############
model_w<- glmmTMB(
  waitingtime_mean ~ puremixed + (1|colony),
  data = df_wide,
  family = Gamma(link = log)# 'ensures is positive', right skewed (link = log)
  
)

sim <- DHARMa::simulateResiduals(model_w)
plot(sim)

BIC(model_w)

############
#waitingtime_sd
###############
model_w<- glmmTMB(
  waitingtime_sd ~ puremixed*genotype + (1|colony/ant),
  data = df_wide,
  family = Gamma(link = log)# 'ensures is positive', right skewed (link = log)
  
)

sim <- DHARMa::simulateResiduals(model_w)
plot(sim)

AIC(model_w)

pairs <- emmeans(model_w, pairwise ~ genotype*puremixed, type = "response")
stat_test <- pairs$contrasts


model_w1<- glmmTMB(
  waitingtime_sd ~ puremixed*genotype + (1|colony),
  data = df_wide,
  family = Gamma(link = log)# 'ensures is positive', right skewed (link = log)
  
)

sim <- DHARMa::simulateResiduals(model_w1)
plot(sim)

AIC(model_w1)


##

model_w2<- glmmTMB(
  waitingtime_sd ~ genotype + (1|colony/ant),
  data = df_wide,
  family = Gamma(link = log)# 'ensures is positive', right skewed (link = log)
  
)

sim <- DHARMa::simulateResiduals(model_w2)
plot(sim)
AIC(model_w2)

######
#local efficiency
######


model_e<- glmmTMB(
  local_efficiency ~ genotype*puremixed + (1|colony/ant),
  data = df_wide,
  family = Gamma(link = "log")# 'ensures is positive', right skewed (link = log)
  
)

sim <- DHARMa::simulateResiduals(model_e)
plot(sim)
AIC(model_e)

pairs <- emmeans(model_e, pairwise ~ genotype*puremixed, type = "response")
stat_test <- pairs$contrasts




#############################
####plot antlevel data for all parameters antlevel in a loop
#############################
library(ggplot2)
library(viridisLite)
library(viridis)

treatment_colors<-c(
  "a" = "#05e0fc",
  "b" = "#fc0536",
  "ba" = "#db05fc"
)

output_folder <- "network_parameter_plots/exp1_inf/antlevel/"

for (param in antlevel_parameters){
  param_plot <- ggplot(df_wide, aes(x = colony, y = .data[[param]], fill = treatment)) +
    
    geom_jitter(aes(color = genotype), width = 0.1, alpha = 0.7, size = 2) +
    geom_boxplot(width = 0.3, outlier.shape = NA, color = "black", alpha = 0.4) +
    theme_minimal() +
    scale_fill_manual(values = c("a" = "#05e0fc", "b" = "#fc0536", "ba" = "#db05fc")) +
    scale_color_manual(values = c("a" = "#05e0fc", "b" = "#fc0536", "ba" = "#db05fc")) +
    labs(
      title = " ",
      x = "colony",
      y = param
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # rotate x labels for readability
  # Save each plot as PNG
  ggsave(
    filename = paste0(output_folder, "/", "ANTLEVEL_" , param, ".png"),
    plot = param_plot,
    width = 6,
    height = 5,
    dpi = 300
  )
  
}


#####################
#trash
for (para in parameters) {
  # Compute Shapiro–Wilk test per treatment for this parameter
  qqplot <- network_parameter_complete %>%
    group_by(treatment) %>%
    summarise(
      shapiro_p = shapiro.test(.data[[para]])$p.value,
      .groups = "drop"
    )
  
  # Store results in list using parameter name
  distribution_test[[para]] <- res
}

for (para in parameters) {
  # Compute Shapiro–Wilk test per treatment for this parameter
  res <- network_parameter_complete %>%
    group_by(treatment) %>%
    summarise(
      shapiro_p = shapiro.test(.data[[para]])$p.value,
      .groups = "drop"
    )
  
  # Store results in list using parameter name
  distribution_test[[para]] <- res
}
