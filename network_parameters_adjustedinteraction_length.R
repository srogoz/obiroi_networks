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
                         #"b16-5" = c("BB", "BG", "BO", "BP", "GB", "GG", "GO", "GP", "OB", "OG", "OO", "OP", "PB", "PG", "PO", "PP"),
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

network_obj_list<-list("a16-1" = "edge_lists/try/a16-1_edgelist_03102025.rds",
                       "a16-2" = "edge_lists/try/a16-2_edgelist_03102025.rds",
                       "a16-3" = "edge_lists/try/a16-3_edgelist_03102025.rds",
                       "a16-4" = "edge_lists/try/a16-4_edgelist_03102025.rds",
                       "a16-5" = "edge_lists/try/a16-5_edgelist_03102025.rds",
                       "a16-6" = "edge_lists/try/a16-6_edgelist_03102025.rds",
                       "a16-7" = "edge_lists/try/a16-7_edgelist_06102025.rds",
                       "a16-8" = "edge_lists/try/a16-8_edgelist_06102025.rds",
                       "b16-1" = "edge_lists/try/b16-1_edgelist_03102025.rds",
                       "b16-2" = "edge_lists/try/b16-2_edgelist_03102025.rds",
                       "b16-3" = "edge_lists/try/b16-3_edgelist_03102025.rds",
                       "b16-4" = "edge_lists/try/b16-4_edgelist_03102025.rds",
                       #"b16-5" = "",
                       "b16-6" = "edge_lists/try/b16-6_edgelist_03102025.rds",
                       "b16-7" = "edge_lists/try/b16-7_edgelist_06102025.rds",
                       "b16-8" = "edge_lists/try/b16-8_edgelist_06102025.rds",
                       "ba16-1" = "edge_lists/try/ba16-1_edgelist_03102025.rds",
                       "ba16-2" = "edge_lists/try/ba16-2_edgelist_03102025.rds",
                       "ba16-3" = "edge_lists/try/ba16-3_edgelist_03102025.rds",
                       "ba16-4" = "edge_lists/try/ba16-4_edgelist_03102025.rds",
                       "ba16-5" = "edge_lists/try/ba16-5_edgelist_03102025.rds",
                       "ba16-6" = "edge_lists/try/ba16-6_edgelist_03102025.rds",
                       "ba16-7" = "edge_lists/try/ba16-7_edgelist_06102025.rds",
                       "ba16-8" = "edge_lists/try/ba16-8_edgelist_06102025.rds"
)


selected_colonies<-list("a16-1", "a16-2", "a16-3", "a16-4", "a16-5", "a16-6",
                        "b16-1", "b16-2", "b16-3", "b16-4", "b16-6",
                        "ba16-1", "ba16-2", "ba16-3", "ba16-4", "ba16-5", "ba16-6",
                        "a16-7","a16-8","b16-7","b16-8","ba16-7","ba16-8"
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
parameters <-c("total_number_interactions", "total_sum_of_interactions", "strength_mean", "strength_sd", "interaction_length_mean", "interaction_length_sd", "waiting_time_mean", "waiting_time_sd", "burstiness", "density_aggr")
#single ant level parameters
#parameters <-c("waiting_time_ant_mean", "waiting_time_ant_sd", "burstiness_ant_sd", "burstiness_ant_mean", "inter_perant_mean", "inter_perant_sd" )


#########
#organisational

folder_path_networkpara <-"network_parameters/interaction_1_30s/"
output_folder <- "network_parameter_plots/interaction_1_30s/"
timestamp <- format(Sys.time(), "%d%m%Y")

#get network parameters
#######################################

for (i in 1:length(selected_colonies)){
  
   colony_name<-selected_colonies[[i]]
  framerate_col<-framerate[[colony_name]]
  lower_limit_interaction<-interaction_limit_sec*framerate_col
  higher_limit_interaction<-interaction_limit_min*60*framerate_col
    
  
  #prepare network obj, interaction limits
  network_obj <-readRDS(network_obj_list[[colony_name]])
  
  network_obj_5 <-network_obj[(network_obj$duration>lower_limit_interaction & network_obj$duration<higher_limit_interaction),]
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

network_parameter_collected_with_treatment<-data.frame(treatment = rep(c("a", "b", "ba"), each = 8), replicate = rep(1:8, times = 3))
network_parameter_collected_with_treatment <- 
  network_parameter_collected_with_treatment %>%
  filter(!(treatment == "b" & replicate == 5))

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
    filename = paste0(output_folder, "/", param, "_violin_boxplot.png"),
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










