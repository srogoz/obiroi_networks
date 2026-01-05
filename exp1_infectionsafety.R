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

#list of adj-matrix files
adjmatrix_list<-list(
  "a16-1" = "adj_matrix/exp2/a16-1_proxMatrInterpo_20251229.npy",
  "a16-2" = "adj_matrix/exp2/a16-2_proxMatrInterpo_20251229.npy",
  "a16-3" = "adj_matrix/exp2/a16-3_proxMatrInterpo_20251229.npy",
  "a16-4" = "adj_matrix/exp2/a16-4_proxMatrInterpo_20251229.npy",
  "a16-5" = "adj_matrix/exp2/a16-5_proxMatrInterpo_20260105.npy",
  "a16-6" = "adj_matrix/exp2/a16-6_proxMatrInterpo_20251229.npy",
  "a16-7" = "adj_matrix/exp2/a16-7_proxMatrInterpo_20251231.npy",
  "a16-8" = "adj_matrix/exp2/a16-8_proxMatrInterpo_20251231.npy",
  "b16-1" = "adj_matrix/exp2/b16-1_proxMatrInterpo_20251229.npy",
  #"b16-2" = "adj_matrix/b16-y",
  "b16-3" = "adj_matrix/exp2/b16-3_proxMatrInterpo_20251229.npy",
  "b16-4" = "adj_matrix/exp2/b16-4_proxMatrInterpo_20251229.npy",
  "b16-5" = "adj_matrix/exp2/b16-5_proxMatrInterpo_20251229.npy",
  "b16-6" = "adj_matrix/exp2/b16-6_proxMatrInterpo_20260105.npy",
  "b16-7" = "adj_matrix/exp2/b16-7_proxMatrInterpo_20251231.npy",
  "b16-8" = "adj_matrix/exp2/b16-8_proxMatrInterpo_20251231.npy",
  "ba16-1" = "adj_matrix/exp2/ba16-1_proxMatrInterpo_20251229.npy",
  #"ba16-2" = "adj_matrix/",
  "ba16-3" = "adj_matrix/exp2/ba16-3_proxMatrInterpo_20260105.npy",
  "ba16-4" = "adj_matrix/exp2/ba16-4_proxMatrInterpo_20251229.npy",
  "ba16-5" = "adj_matrix/exp2/ba16-5_proxMatrInterpo_20251229.npy",
  "ba16-6" = "adj_matrix/exp2/ba16-6_proxMatrInterpo_20260105.npy",
  "ba16-7" = "adj_matrix/exp2/ba16-7_proxMatrInterpo_20260105.npy",
  "ba16-8" = "adj_matrix/exp2/ba16-8_proxMatrInterpo_20251231.npy"
)


present_ants_list<-list("a16-1" = "present_ants/exp2/a16-1_present_ants_20251229.npy",
                        "a16-2" = "present_ants/exp2/a16-2_present_ants_20251229.npy",
                        "a16-3" = "present_ants/exp2/a16-3_present_ants_20251229.npy",
                        "a16-4" = "present_ants/exp2/a16-4_present_ants_20251229.npy",
                        "a16-5" = "present_ants/exp2/a16-5_present_ants_20260105.npy",
                        "a16-6" = "present_ants/exp2/a16-6_present_ants_20251229.npy",
                        "a16-7" = "present_ants/exp2/a16-7_present_ants_20251231.npy",
                        "a16-8" = "present_ants/exp2/a16-8_present_ants_20251231.npy",
                        "b16-1" = "present_ants/exp2/b16-1_present_ants_20251229.npy",
                       # "b16-2" = "present_ants/b16-2_present_ants_20250930.npy",
                        "b16-3" = "present_ants/exp2/b16-3_present_ants_20251229.npy",
                        "b16-4" = "present_ants/exp2/b16-4_present_ants_20251229.npy",
                        "b16-5" = "present_ants/exp2/b16-5_present_ants_20251229.npy",
                        "b16-6" = "present_ants/exp2/b16-6_present_ants_20260105.npy",
                        "b16-7" = "present_ants/exp2/b16-7_present_ants_20251231.npy",
                        "b16-8" = "present_ants/exp2/b16-8_present_ants_20251231.npy",
                        "ba16-1" = "present_ants/exp2/ba16-1_present_ants_20251229.npy",
                        #"ba16-2" = "present_ants/ba16-2_present_ants_20250926.npy",
                        "ba16-3" = "present_ants/exp2/ba16-3_present_ants_20260105.npy",
                        "ba16-4" = "present_ants/exp2/ba16-4_present_ants_20251229.npy",
                        "ba16-5" = "present_ants/exp2/ba16-5_present_ants_20251229.npy",
                        "ba16-6" = "present_ants/exp2/ba16-6_present_ants_20260105.npy",
                        "ba16-7" = "present_ants/exp2/ba16-7_present_ants_20260105.npy",
                        "ba16-8" = "present_ants/exp2/ba16-8_present_ants_20251231.npy"
)

selected_colonies<-list(
  #"a16-1", "a16-2", "a16-3", "a16-4",
                        "a16-5",
                        #"a16-6",
                       #"a16-7",
                       #"a16-8",
                       #"b16-1",
                       #"b16-2", 
                       #"b16-3", "b16-4", "b16-5",
                       "b16-6", 
                       #"b16-7","b16-8",
                        #"ba16-1",
                       #"ba16-2",
                       "ba16-3", 
                       #"ba16-4", "ba16-5", 
                       "ba16-6","ba16-7"
                       #"ba16-8"
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
#create edgelists
# Load .npy file function
np <- import("numpy")

#edgelist function
interaction_dist <- function(net, present_ants) {
  edgelist <- data.frame(
    head = character(),
    tail = character(),
    onset = integer(),
    terminus = integer(),
    duration = integer(),
    stringsAsFactors = FALSE
  )
  
  n_ants <- length(present_ants)
  
  for (i in seq_len(n_ants)) {
    for (j in seq_len(n_ants)) {
      if (i < j) {  # skip self-interactions and dont run over bounds twice since undirected matrix
        # Get interaction vector over time
        interaction_vec <- net[i, j, ]
        
        # Run-length encoding to detect continuous stretches
        runs <- rle(interaction_vec)
        cumulative_time <- cumsum(runs$lengths)
        
        for (k in seq_along(runs$values)) {
          if (runs$values[k] == 1) {
            onset <- ifelse(k == 1, 1, cumulative_time[k - 1] + 1)
            terminus <- cumulative_time[k]
            duration <- runs$lengths[k]
            
            edgelist <- rbind(
              edgelist,
              data.frame(
                head = present_ants[i],
                tail = present_ants[j],
                onset = onset,
                terminus = terminus,
                duration = duration,
                stringsAsFactors = FALSE
              )
            )
          }
        }
      }
    }
  }
  
  return(edgelist)
}

for (i in 1:length(selected_colonies)){
  network<-np$load(adjmatrix_list[[selected_colonies[[i]]]])
  print(dim(network)[3])
}
for (i in 1:length(selected_colonies)){
  
  
  colony_name<-selected_colonies[[i]]
  network<-np$load(adjmatrix_list[[selected_colonies[[i]]]])
  #network <- np$load("adj_matrix/b16-2_proxMatrInterpo_20250930.npy")
  #aggregated colony networks reduced to 12hs in python and 
  present_ants <- np$load(present_ants_list[[selected_colonies[[i]]]], allow_pickle= TRUE)
  colnames(network)<-present_ants
  row.names(network)<-present_ants
  ###########
  #create edgelist from interaction array
  # transform into network object, undirected
  #reduce  all networks to 12hs: 432000 frames, starting from the end
  
  if (framerate[[selected_colonies[[i]]]]==10){
    network_cut<-network[,,1:47062]}
  else
   { network_cut<-network[,,1:23531]}
  
 
  
  
network_obj <-interaction_dist(network_cut, present_ants)

#save edgelist
timetag<-format(Sys.time(), "%d%m%Y")
saveRDS(network_obj,
        paste0("edge_lists/exp2/",colony_name, "_edgelist_", timetag, ".rds"))

print(paste0(colony_name, "done"))

}

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

folder_path_networkpara <-"network_parameters/antlevel_1_truncated_2min//"
output_folder <- "network_parameter_plots/antlevel_1_truncated_2min/"
timestamp <- format(Sys.time(), "%d%m%Y")

#get network parameters
#######################################
#antlevel parameter collection
strength_antlevel_mean<-data.frame(row.names = expected_ants)
strength_antlevel_sd<-data.frame(row.names = expected_ants)
local_efficiency<-data.frame(row.names = expected_ants)
waiting_time_ant_mean_df<-data.frame(row.names = expected_ants)
waiting_time_ant_sd_df<-data.frame(row.names = expected_ants)
interactionlength_ant_mean<-data.frame(row.names = expected_ants)
burstiness_ant_df<-data.frame(row.names = expected_ants)

for (i in 1:length(selected_colonies)){
  
  colony_name<-selected_colonies[[1]]
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
  
  p_hist <- ggplot(network_obj, aes(duration)) +
    geom_histogram(
      aes(fill = after_stat(
        ifelse(
          x < lower_limit_interaction, "below",
          ifelse(x > higher_limit_interaction, "above", "inside")
        )
      )),
      bins = 30,
      color = "grey"
    )  +
    scale_fill_manual(
      values = c(
        "below" = "grey",
        "inside" = "magenta",
        "above" = "seagreen2"
      ),
      labels = c(
        below  = paste0("< 1s (removed)\n(", lower_cut_framepercentage ," /",lower_cut_interactionumberpercentage, " )"),
        inside = paste0("unaffected"),
        above  = paste0("< 2 min (truncated)\n(", upper_cut_framepercentage, " /" ,upper_cut_interactionumberpercentage, " )")
      ),
      name = "Thresholds \n(% of frames, \ninteractions affected)"
    ) +
    scale_y_log10() +   # ✅ log-scale on y-axis
    geom_text(
      stat = "bin",
      bins = 40,
      aes(label = after_stat(count)),
      vjust = -0.5,
      size = 3
    ) +
    theme_minimal() +
    labs(fill = paste0("Bin midpoint > ", lower_limit_interaction),
         y = "Count (log scale)")+
    ggtitle(paste0("Interactionlength histogram for ", colony_name))+
    annotate(
      "text",
      x = Inf, y = Inf,
      hjust = 1.1, vjust = 2,
      label = paste0(" total frames :" ,sum_discarded_frames ," %", "\ntotal interactions: " ,affected_interactions ," %" ),
      size = 4
    ) +
    theme(
      # smaller labels
      legend.key.height = unit(2.2, "lines"),        # more vertical spacing
      # more horizontal spacing
    )
  
  ggsave(paste0("network_parameter_plots/antlevel_1_truncated_2min/interaction_histograms/hist_interactions_",colony_name,".png"), 
         plot = p_hist,
         width = 6,
         height = 5,
         dpi = 300)
  
  
  #########
  present_ants<-unique(c(network_obj_5$tail,network_obj_5$head))
  network_obj_5$tail<-match(network_obj_5$tail, present_ants)
  network_obj_5$head<-match(network_obj_5$head, present_ants)
  
  #aggregate network with new interation restrictions
  aggregated_network <- aggregate_from_edgelist(network_obj_5)
  
  aggregated_network<-aggregate_from_edgelist(network_obj_5)
  aggregated_network[,16]<-t(aggregated_network[16,])
  aggregated_network[16,]<-0
  
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
  number_ants<-length(present_ants_list[[colony_name]])
  total_number_interactions<-length(network_obj_5$duration)/framerate_col
  #results_networkparameters$total_number_interactions<-total_number_interactions
  total_sum_of_interactions<-sum(network_obj_5$duration)/framerate_col
  #frame change, 10 frames and 5 frames
  duration<-12*60*60*framerate_col
  density_aggr<-total_sum_of_interactions/(((number_ants*(number_ants-1))/2)*duration)
  
  #DEGREE HETEROGENITY MEASURE
  strength_mean<-mean(strength(graph_aggregated))
  strength_sd<-sd(strength(graph_aggregated))
  strength_antlevel_mean[[colony_name]]<-strength(graph_aggregated)
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
  
  interactionlength_ant_mean[[colony_name]]<-interactionlength_means
  waiting_time_ant_mean_df[[colony_name]]<-waiting_time_ant_mean
  waiting_time_ant_sd_df[[colony_name]]<-waiting_time_ant_sd
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

