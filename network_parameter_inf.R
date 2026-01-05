library(igraph)
library(sna)
library(tsna)
#library(ndtv)
library(reticulate)
library(igraph)
library(intergraph)
library(networkLite)
library(network)
library(networkDynamic)
library(dplyr)
library(glmmTMB)
library(tidyr)
library(ggplot2)
library(emmeans)
library(lme4)

#PULL DATA: present ants, edgelists, treatment, framerate and selected colonies for analysis
################################
present_ants_list<-list( "a16-1" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "a16-2" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "a16-3" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "a16-4" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "a16-5" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "a16-6" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "a16-7" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "a16-8" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "b16-1" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "b16-2" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "b16-3" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "b16-4" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "b16-5" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "b16-6" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "b16-7" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "b16-8" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "ba16-1" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "ba16-2" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "ba16-3" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "ba16-4" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "ba16-5" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "ba16-6" = c("BB", "BG", "GB", "GG", "GO", "OB", "OG", "OP", "PB", "PG", "PO", "PP"),
                         "ba16-7" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
                         "ba16-8" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP")
)

network_obj_list<-list("a16-1" = "edge_lists/exp2/a16-1_edgelist_02012026.rds",
                       "a16-2" = "edge_lists/exp2/a16-2_edgelist_02012026.rds",
                       "a16-3" = "edge_lists/exp2/a16-3_edgelist_02012026.rds",
                       "a16-4" = "edge_lists/exp2/a16-4_edgelist_02012026.rds",
                       "a16-5" = "edge_lists/exp2/a16-5_edgelist_05012026.rds",
                       "a16-6" = "edge_lists/exp2/a16-6_edgelist_02012026.rds",
                       "a16-7" = "edge_lists/exp2/a16-7_edgelist_02012026.rds",
                       "a16-8" = "edge_lists/exp2/a16-8_edgelist_02012026.rds",
                       "b16-1" = "edge_lists/exp2/b16-1_edgelist_02012026.rds",
                       #"b16-2" = "edge_lists/try/b16-2_edgelist_03102025.rds",
                       "b16-3" = "edge_lists/exp2/b16-3_edgelist_02012026.rds",
                       "b16-4" = "edge_lists/exp2/b16-4_edgelist_02012026.rds",
                       "b16-5" = "edge_lists/exp2/b16-5_edgelist_02012026.rds",
                       "b16-6" = "edge_lists/exp2/b16-6_edgelist_05012026.rds",
                       "b16-7" = "edge_lists/exp2/b16-7_edgelist_02012026.rds",
                       "b16-8" = "edge_lists/exp2/b16-8_edgelist_02012026.rds",
                       "ba16-1" = "edge_lists/exp2/ba16-1_edgelist_02012026.rds",
                       #"ba16-2" = "edge_lists/exp2/ba",
                       "ba16-3" = "edge_lists/exp2/ba16-3_edgelist_05012026.rds",
                       "ba16-4" = "edge_lists/exp2/ba16-4_edgelist_02012026.rds",
                       "ba16-5" = "edge_lists/exp2/ba16-5_edgelist_02012026.rds",
                       "ba16-6" = "edge_lists/exp2/ba16-6_edgelist_05012026.rds",
                       "ba16-7" = "edge_lists/exp2/ba16-7_edgelist_05012026.rds",
                       "ba16-8" = "edge_lists/exp2/ba16-8_edgelist_02012026.rds"
)


selected_colonies<-list("a16-1", "a16-2", "a16-3", "a16-4", "a16-5", "a16-6",
                        "b16-1", 
                        #"b16-2",
                        "b16-3", "b16-4","b16-5",
                        "b16-6",
                        "ba16-1",
                        #"ba16-2", 
                        "ba16-3", "ba16-4", "ba16-5", 
                        "ba16-6",
                        "a16-7","a16-8","b16-7","b16-8" ,"ba16-7","ba16-8"
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
                "b16-5"= 10,
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
###################
#FUNCTIONS
###################

#density
get_density<-function(network, present_ants){
  apply(network, 3, function(net_frame){
    
    net_frame_obj<- graph_from_adjacency_matrix(net_frame, weighted = NULL, add.colnames = NULL) %>% set_vertex_attr("name", value = present_ants)
    edge_density(net_frame_obj, loops = FALSE)
    
  })
}

#matching files in directory
get_colony_files<-function(folder,selected_colonies){
  
  data_list<-list()
  files<-list.files(folder,full.names = TRUE)
  for (i in 1:length(selected_colonies)){
    colony<-selected_colonies[[i]]
    pattern <- paste0("^", colony, "_")
    matching_files <- files[grepl(pattern, basename(files))]
    data_list[[colony]]<- matching_files
  }
  
  return(data_list)
}

aggregate_from_edgelist <- function(network_obj) {
  # Filter durations
  net<- network_obj
  # All unique ants from both head and tail
  ants <- sort(unique(c(net$head, net$tail)))
  
  # Aggregate
  agg <- with(net, tapply(duration, list(head, tail), sum))
  
  # Create a full square matrix with all ants
  agg_full <- matrix(0, nrow = length(ants), ncol = length(ants),
                     dimnames = list(ants, ants))
  
  # Fill in overlapping entries
  common_rows <- intersect(rownames(agg), ants)
  common_cols <- intersect(colnames(agg), ants)
  
  agg_full[common_rows, common_cols] <- agg[common_rows, common_cols]
  
  # Replace NAs with 0 just in case
  agg_full[is.na(agg_full)] <- 0
  
  return(agg_full)
}


##################
#PARAMETERS
##################
#functional
interaction_limit_min<-2
interaction_limit_sec<-1
#plots

treatment_colors<-c(
  "a" = "#05e0fc",
  "b" = "#fc0536",
  "ba" = "#db05fc"
)

# overview parameters
parameters <-c("total_number_interactions", "total_sum_of_interactions", "strength_mean", "strength_sd", "interaction_length_mean", "interaction_length_sd", "waiting_time_mean", "waiting_time_sd", "burstiness", "density_aggr", "global_eff", "assortativity")
#single ant level parameters
#parameters <-c("waiting_time_ant_mean", "waiting_time_ant_sd", "burstiness_ant_sd", "burstiness_ant_mean", "inter_perant_mean", "inter_perant_sd" )


#########
#organisational

folder_path_networkpara <-"network_parameters/infection/"
output_folder <- "network_parameter_plots/infection/"
timestamp <- format(Sys.time(), "%d%m%Y")

#get network parameters
#######################################
#antlevel parameter collection
#np <- import("numpy")
#network<-np$load(adjmatrix_list[[selected_colonies[[i]]]])
strength_antlevel_mean<-data.frame(row.names = expected_ants)
strength_antlevel_sd<-data.frame(row.names = expected_ants)
local_efficiency<-data.frame(row.names = expected_ants)
waiting_time_ant_mean_df<-data.frame(row.names = expected_ants)
waiting_time_ant_sd_df<-data.frame(row.names = expected_ants)
interactionlength_ant_mean<-data.frame(row.names = expected_ants)
burstiness_ant_df<-data.frame(row.names = expected_ants)

for (i in 1:length(selected_colonies)){
  
  colony_name<-selected_colonies[[i]]
  framerate_col<-framerate[[colony_name]]
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
  ##############
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
  #   scale_y_log10(limits = c(1, NA)) + 
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
  # ggsave(paste0("network_parameter_plots/infection/interaction_histograms/hist_interactions_",colony_name,".png"),
  #        plot = p_hist,
  #        width = 6,
  #        height = 5,
  #        dpi = 300)
  # 
  # 
  #########
  present_ants<-unique(c(network_obj_5$tail,network_obj_5$head))
  network_obj_5$tail<-match(network_obj_5$tail, present_ants)
  network_obj_5$head<-match(network_obj_5$head, present_ants)
  
  #aggregate network with new interation restrictions
  aggregated_network <- aggregate_from_edgelist(network_obj_5)
  
  aggregated_network<-aggregate_from_edgelist(network_obj_5)
  matrix_dim<- max (dim(aggregated_network))
  aggregated_network[,matrix_dim]<-t(aggregated_network[matrix_dim,])
  aggregated_network[matrix_dim,]<-0
  
  network_aggregated_symm<-aggregated_network
  network_aggregated_symm[lower.tri(network_aggregated_symm)]<-t(aggregated_network)[lower.tri(aggregated_network)]
  network_aggregated_symm<-round(network_aggregated_symm/(36000),2)
  #into igraph object
  graph_aggregated<-graph_from_adjacency_matrix(network_aggregated_symm, weighted = TRUE, mode = "undirected", add.colnames = NULL)
  
  #net_density <- get_density(network_obj, present_ants)
  #density_temporal_mean<-mean(net_density)
  #density_temporal_sd<-sd(net_density)
  
  ##################################################
  #AGGREGATED MEASURES
  number_ants<-length(present_ants)
  total_number_interactions<-length(network_obj_5$duration)/framerate_col
  #results_networkparameters$total_number_interactions<-total_number_interactions
  total_sum_of_interactions<-sum(network_obj_5$duration)/framerate_col
  #frame change, 10 frames and 5 frames
  duration<-12*60*60*framerate_col
  density_aggr<-total_sum_of_interactions/(((number_ants*(number_ants-1))/2)*duration)
  
  #DEGREE HETEROGENITY MEASURE
  strength_mean<-mean(strength(graph_aggregated))
  strength_sd<-sd(strength(graph_aggregated))
  strength_collection<-c(strength(graph_aggregated))
  #padding length in case ants are missing
  length(strength_collection) <- length(expected_ants)
  strength_antlevel_mean[[colony_name]]<-strength_collection
  V(graph_aggregated)$strength <- strength(graph_aggregated, weights = E(graph_aggregated)$weight)
  ass_ants<-assortativity(graph_aggregated, types1   = V(graph_aggregated)$strength,
                          directed = FALSE
  )
  
  ####Local and global Efficiency
  global_eff<-global_efficiency(graph_aggregated, weights = E(graph_aggregated)$weight)
  
  local_eff<-local_efficiency(
    graph_aggregated,
    vids = V(graph_aggregated),
    weights = E(graph_aggregated)$weight,
    mode = c("all")
  )
  length(local_eff) <- length(expected_ants)
  local_efficiency[[colony_name]]<-local_eff
  
  #TEMPORAL MEASURES
  #interaction durations
  interaction_length_mean<-mean(network_obj_5$duration)
  interaction_length_sd<-sd(network_obj_5$duration)
  
  #waiting times
  waiting_times<-data.frame(wt=diff(sort(network_obj_5$onset)) )
  wait <- waiting_times %>%
    filter(wt > 0)
  waiting_times <- waiting_times[waiting_times$wt > 0, , drop = FALSE]
  waiting_times_5<-data.frame(wt=diff(sort(network_obj_5$onset)))
  waiting_times_5<-waiting_times_5[waiting_times_5$wt>0,,drop = FALSE]/framerate_col
  max_wt<- max(wait$wt)
  max_wt_5<- max(waiting_times_5$wt)
  ##
  waiting_time_mean<-mean(unlist(waiting_times_5))
  waiting_time_sd<-sd(unlist(waiting_times_5))
  #burstyness
  burstiness<-(waiting_time_sd-waiting_time_mean)/(waiting_time_sd+waiting_time_mean)
  
  
  ######################
  #on ant level
  interactionlength_means<-c()
  waiting_time_ant_mean<-c()
  waiting_time_ant_sd<-c()
  burstiness_ant<-c()
  
  for (i in 1:length(present_ants)){
    df_ant <- network_obj_5[network_obj_5$head == i , ]
    interactionlength_means[i]<-mean(df_ant$duration)
    
    waiting_times<-data.frame(wt=diff(sort(df_ant$onset)) )
    wait <- waiting_times %>%
      filter(wt > 0)
    waiting_times <- waiting_times[waiting_times$wt > 0, , drop = FALSE]
    waiting_times_5<-data.frame(wt=diff(sort(df_ant$onset)))
    waiting_times_5<-waiting_times_5[waiting_times_5$wt>0,,drop = FALSE]
    
    waiting_time_ant_mean[i]<-mean(unlist(waiting_times_5))
    waiting_time_ant_sd[i]<-sd(unlist(waiting_times_5))
    #burstyness
    burstiness_ant[i]<-(waiting_time_ant_sd[i]-waiting_time_ant_mean[i])/(waiting_time_ant_sd[i]+waiting_time_ant_mean[i])
  }
  
  length(interactionlength_means) <- length(expected_ants)
  interactionlength_ant_mean[[colony_name]]<-interactionlength_means
  length(waiting_time_ant_mean) <- length(expected_ants)
  waiting_time_ant_mean_df[[colony_name]]<-waiting_time_ant_mean
  length(waiting_time_ant_sd) <- length(expected_ants)
  waiting_time_ant_sd_df[[colony_name]]<-waiting_time_ant_sd
  length(burstiness_ant) <- length(expected_ants)
  burstiness_ant_df[[colony_name]]<-burstiness_ant
  
  #########################single ant level
  # interactionlength_means<-c()
  # waiting_time_ant_mean<-c()
  # waiting_time_ant_sd<-c()
  # burstiness_ant<-c()
  
  # network_obj$head <- as.character(network_obj$head)
  # network_obj$tail <- as.character(network_obj$tail)
  # for (i in length(present_ants)){
  #   df_ant <- network_obj_5[network_obj_5$head == i , ]
  #   interactionlength_means[i]<-mean(df_ant$duration)
  #   
  #   waiting_times<-data.frame(wt=diff(sort(df_ant$onset)) )
  #   wait <- waiting_times %>%
  #     filter(wt > 0)
  #   waiting_times <- waiting_times[waiting_times$wt > 0, , drop = FALSE]
  #   waiting_times_5<-data.frame(wt=diff(sort(df_ant$onset)))
  #   waiting_times_5<-waiting_times_5[waiting_times_5$wt>0,,drop = FALSE]
  #   
  #   waiting_time_ant_mean[i]<-mean(unlist(waiting_times_5))
  #   waiting_time_ant_sd[i]<-sd(unlist(waiting_times_5))
  #   #burstyness
  #   burstiness_ant[i]<-(waiting_time_ant_sd[i]-waiting_time_ant_mean[i])/(waiting_time_ant_sd[i]+waiting_time_ant_mean[i])
  #   
  # }
  # 
  # ##cross ant comparisons for colonies
  # inter_perant_mean<-mean(interactionlength_means)
  # inter_perant_sd<-sd(interactionlength_means)
  # burstiness_ant_sd<-sd(burstiness_ant)
  # burstiness_ant_mean<-mean(burstiness_ant)
  # waiting_time_colony_mean<-mean(waiting_time_ant_mean)
  # waiting_time_colony_sd<-sd(waiting_time_ant_mean)
  # 
  
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
                                        
                                        
                                        burstiness = burstiness
                                        # waiting_time_ant_mean = waiting_time_ant_mean,
                                        # waiting_time_ant_sd = waiting_time_ant_sd,
                                        # burstiness_ant = burstiness_ant,
                                        # inter_perant_mean=inter_perant_mean,
                                        # inter_perant_sd=inter_perant_sd,
                                        # burstiness_ant_sd=burstiness_ant_sd,
                                        # burstiness_ant_mean=burstiness_ant_mean,
                                        # waiting_time_colony_mean=waiting_time_colony_mean,
                                        # waiting_time_colony_sd=waiting_time_colony_sd
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

###############
###stats and plots
###############

network_parameter_list<-get_colony_files(folder_path_networkpara, selected_colonies)

###########
#prepare network parameters
network_parameter_collected<-list()

for (i in 1:length(selected_colonies)){
  network_parameter_col<-read.csv(network_parameter_list[[selected_colonies[[i]]]])
  network_parameter_col<-network_parameter_col[1,]
  network_parameter_collected<-rbind(network_parameter_collected, network_parameter_col)
}
###adjust treatment and replicates
#treatments a,b,ba
#network_parameter_collected_with_treatment<-data.frame(treatment = rep(c("a", "b", "ba"), each = 8), replicate = rep(1:8, times = 3))
#treatment only b, a 
network_parameter_collected_with_treatment<-data.frame(treatment = rep(c("a", "b", "ba"), each = 8), replicate = rep(1:8, times = 3))

##remove missing replicates
network_parameter_collected_with_treatment <- 
  network_parameter_collected_with_treatment %>%
  filter(!(treatment == "b" & replicate == 2))

network_parameter_collected_with_treatment <- 
  network_parameter_collected_with_treatment %>%
 filter(!(treatment == "ba" & replicate == 2))

network_parameter_complete<-data_frame(network_parameter_collected_with_treatment, network_parameter_collected)

##################
##################
library(ggplot2)
library(viridisLite)
library(viridis)
################
##TEST FOR NORMAL DISTRIBUTION

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
kruskal_test_results<-list()

for (para in parameters){
  formula <- as.formula(paste(para, "~ treatment"))
  
  res<-kruskal.test(formula, data = network_parameter_complete)
  kruskal_test_results[[para]]<-res
}

#######################
#plot violin plots for network parameters

for (param in parameters){
  p_val<-round(kruskal_test_results[[param]]$p.value,2)
  param_plot<-ggplot(network_parameter_complete, 
                     aes(x = treatment, y = .data[[param]], fill = treatment)) +
    geom_violin(trim = FALSE, alpha = 0.4) +  # “bell shape”
    geom_boxplot(width = 0.15, outlier.shape = NA, color = "black") +
    geom_jitter(width = 0.1, alpha = 0.7) +   # individual points
    geom_text(aes(label = replicate), vjust = -1, size = 3.5, color = "black") +
    theme_minimal() +
    # Apply manual colors
    scale_fill_manual(values = treatment_colors) +
    labs(
      title = paste0(param, " p = ", p_val ),
      x = "Treatment",
      y = param) 
  #  scale_fill_viridis_d(option = "D")
  
  # Save each plot as PNG
  ggsave(
    filename = paste0(output_folder, "/", param, "_violin_boxplot1.png"),
    plot = param_plot,
    width = 6,
    height = 5,
    dpi = 300
  )
}



##################################
##all normal distr}buted
anova_test_results <- list()

for (para in parameters) {
  # Build formula dynamically
  formula <- as.formula(paste(para, "~ treatment"))
  
  # Run ANOVA
  res <- aov(formula, data = network_parameter_complete)
  
  # Store full ANOVA summary for that parameter
  #anova_test_results[[para]] <- summary(res)
  anova_test_results[[para]]<-summary(res)[[1]][["Pr(>F)"]][1]
}






###################################
##glmm 

#read as df
colony_activity_singleant<-read.csv(file = "network_parameters/infection/ANTLEVEL_strength_mean05012026.csv")
colony_activity_singleant<-read.csv(file = "network_parameters/infection/ANTLEVEL_burstiness05012026.csv")
colony_activity_singleant<-read.csv(file = "network_parameters/infection/ANTLEVEL_waitingtime_mean05012026.csv")
colony_activity_singleant<-read.csv(file = "network_parameters/infection/ANTLEVEL_local_efficiency05012026.csv")
colony_activity_singleant<-read.csv(file = "network_parameters/infection/ANTLEVEL_interactionlength_mean05012026.csv")
colony_activity_singleant$ant<-expected_ants

#sort data for glmm

colony_activity_singleant <- colony_activity_singleant %>%
  mutate(across(-ant, as.numeric)) 

colony_activity_singleant_glmm<-
  colony_activity_singleant <- colony_activity_singleant %>%
  mutate(across(-ant, as.numeric)) 


df_glmm <- colony_activity_singleant %>%
  pivot_longer(
    cols = -ant,          # all columns except the ant ID
    names_to = "colony",
    values_to = "rsdm"
  ) %>%
  
  mutate(
    treatment = sub("\\d+.*$", "", colony)   # extract 'a', 'b', or 'ba'
  )



df_glmm <- df_glmm %>%
  mutate(
    genotype = case_when(
      treatment == "ba" ~ unlist(genotype_match[ant]),  # lookup in genotype_match
      treatment == "a" ~ "a",                            # all ants in pure a colony
      treatment == "b" ~ "b"                             # all ants in pure b colony
    )
  )

df_glmm <- df_glmm %>%
  mutate(
    puremixed = case_when(
      treatment == "ba" ~ "m",  # lookup in genotype_match
      treatment == "a" ~ "p",                            # all ants in pure a colony
      treatment == "b" ~ "p"                             # all ants in pure b colony
    )
  )

png(file="network_parameter_plots/infection/hist_waitingtime.png",
    width=600, height=350)
hist(df_glmm$rsdm, xlab = "strength [hs]", col = "purple", main="")
dev.off()

#############
#apply glmm

model <- glmmTMB(
  rsdm ~ genotype * puremixed + (1 | colony),
  data = df_glmm,
 # family = Gamma(link = log)# 'ensures is positive', right skewed (link = log)
 family = gaussian()
)

modeltreatment <- glmmTMB(
  rsdm ~ treatment + (1 | colony),
  data = df_glmm,
  family = Gamma(link = log)# 'ensures is positive', right skewed (link = log)
)




pairs <- emmeans(model, pairwise ~ genotype, type = "response")
stat_test <- pairs$contrasts

lmemodel<-lmer()
##################################################
##################################################
#plot all colonies and individual values
library(ggplot2)
library(viridisLite)
library(viridis)

treatment_colors<-c(
  "a" = "#05e0fc",
  "b" = "#fc0536",
  "ba" = "#db05fc"
)

output_folder <- "network_parameter_plots/infection/"


param_plot <- ggplot(df_glmm, aes(x = colony, y = rsdm, fill = treatment)) +
  
  geom_jitter(aes(color = genotype), width = 0.1, alpha = 0.7, size = 2) +
  geom_boxplot(width = 0.3, outlier.shape = NA, color = "black", alpha = 0.4) +
  theme_minimal() +
  scale_fill_manual(values = c("a" = "#05e0fc", "b" = "#fc0536", "ba" = "#db05fc")) +
  scale_color_manual(values = c("a" = "#05e0fc", "b" = "#fc0536", "ba" = "#db05fc")) +
  labs(
    title = " ",
    x = "colony",
    y = "interactionlength"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # rotate x labels for readability
# Save each plot as PNG
ggsave(
  filename = paste0(output_folder, "/", "ANTLEVEL_interactionlength.png"),
  plot = param_plot,
  width = 6,
  height = 5,
  dpi = 300
)

##############################



