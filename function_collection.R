#matching files in directory to treatment_names
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

############aggregate network from edgelist
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

###########density
get_density<-function(network, present_ants){
  apply(network, 3, function(net_frame){
    
    net_frame_obj<- graph_from_adjacency_matrix(net_frame, weighted = NULL, add.colnames = NULL) %>% set_vertex_attr("name", value = present_ants)
    edge_density(net_frame_obj, loops = FALSE)
    
  })
}

#########create edgelists
#edgelist function
make_edgelist <- function(net, present_ants) {
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

################

##create edgelist as basis for networkDynamic and network object
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