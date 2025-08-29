#library(igraph)
library(sna)
library(tsna)
#library(ndtv)
library(reticulate)

##check network file
file.exists("processed/networks/b16-1_proxMatrInterpo_20250806.npy")
# Load .npy file
np <- import("numpy")

network <- np$load("processed/networks/b16-1_proxMatrInterpo_20250806.npy")
present_ants <- np$load("processed/present_ants/b16-1_present_ants_20250806.npy", allow_pickle= TRUE)

###########
#create edgelist from interaction array

#############
# transform into network object, undirected

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

network_obj <-interaction_dist(network, present_ants)
network_obj_5 <-network_obj[network_obj$duration>5,]
##################
# interaction_length distribution 

library(ggplot2)
# change duration into mins from frames
network_obj$duration_min <- network_obj$duration/600
hist <- ggplot(network_obj, aes(x = duration_min))+
  geom_histogram(binwidth = 0.1, fill ="purple", color ="black" , alpha = 0.7) +
  labs( x ="interaction length [mins]", 
        y = "frequency",
        title = "distribution of interaction lengths"
  ) + 
  coord_cartesian(ylim = c(0, 1000), xlim= c(0,2.5)) + 
  #scale_y_log10(limits = c(1, 100))+
  theme_minimal()

##
library (igraph)

##OUTPUT: clustering coefficient global and per ant and for each ant and each frame
get_ccoef_vec <- function(network, present_ants) {
  apply(network, 3, function(net_frame) {
    #look at network per frame
    rownames(net_frame) <- present_ants
    colnames(net_frame) <- present_ants
    
    net_frame_obj <- graph_from_adjacency_matrix(net_frame, weighted = NULL)
    global <- transitivity(net_frame_obj, type = "global")
    local <- transitivity(net_frame_obj, vid = present_ants, type = "local")
    
    c(global, local)
  })
}

############
#density per time and sum of all

get_density<-function(network, present_ants){
  apply(network, 3, function(net_frame){
    
  net_frame_obj<- graph_from_adjacency_matrix(net_frame, weighted = NULL, add.colnames = NULL) %>% set_vertex_attr("name", value = present_ants)
  edge_density(net_frame_obj, loops = FALSE)
  
  })
}

net_density <- get_density(network, present_ants)
mean(net_density)
sd(net_density)
##################################################
library(ggplot2)
library(dplyr)
net_density_df <- data.frame(net_density)
net_density_df$index <- seq_along(net_density_df$net_density)

ggplot(filter(net_density_df, index >=3000 & index <= 10000), aes(y = net_density, x = index))+
 geom_line(color = "purple", size = 0.3)+
  labs( x ="time [frames]", 
        y = "density",
        title = "Density of connections over time for colony b16-1"
  ) +theme_minimal()
#####################################################
##edge formation 
#####################################################
########### 
library(network)
library(networkDynamic)
library(tsna)
library(dplyr)

####################
# into networkDynamic object: structure of edge_df needs to be: onset, terminus , tail, head // tail, head must be NUMERIC
edge_df <- data.frame(onset=network_obj$onset, terminus=network_obj$terminus, tail=network_obj$tail, head=network_obj$head)
#tail and head must be numeric
edge_df$tail<-match(edge_df$tail, present_ants)
edge_df$head<-match(edge_df$head, present_ants)

##CHECK NUMBER OF ONE FRAME INTERACTIONS
bad_edges <- filter(edge_df,(onset >= terminus))
onset_diff <- diff(bad_edges$onset)

# Identify which differences are < 5
close_onsets <- onset_diff < 5

# Count them
sum(close_onsets)
print(bad_edges)

##########check for duplicates
duplicated_edges <- edge_df_corrected[duplicated(edge_df_corrected[, c("tail","head","onset","terminus")]), ]
duplicated_edges
########## each edge appear mirrowed, due to full su
##############################
# correct one frame interactions
edge_df_corrected <- edge_df %>% filter(onset < terminus)
nd<- networkDynamic(edge.spells = as.matrix(edge_df_corrected)
                    ,create.TEAs = TRUE
                    )

library(dplyr)

# Check for true duplicates of the exact same spell
exact_dupes <- edge_df_corrected %>%
  group_by(tail, head, onset, terminus) %>%
  filter(n() > 1) %>%
  ungroup()

if (nrow(exact_dupes) > 0) {
  message("Found ", nrow(exact_dupes), " exact duplicate edge spells:")
  print(exact_dupes)
} else {
  message("exact duplicate edge spells found.")
}

# Check for overlapping spells for the same direction (not reversed)
overlap_dupes <- edge_df_corrected %>%
  group_by(tail, head) %>%
  arrange(onset, terminus) %>%
  mutate(next_onset = lead(onset),
         next_terminus = lead(terminus)) %>%
  filter(!is.na(next_onset) & onset <= next_onset & terminus > next_onset) %>%
  ungroup()

if (nrow(overlap_dupes) > 0) {
  message("Found overlapping spells for the same edge:")
  print(overlap_dupes)
} else {
  message("no overlapping spells for the same edge.")
}

############
edge_df_corrected_2 <- edge_df_corrected %>%
  distinct(tail, head, onset, terminus, .keep_all = TRUE)

######################### edge formation and dissolution

library(tsna)
##for one ant only
nd_singleant <- network.extract(nd, v = 1, neighborhood= "combined")
diso1<-tEdgeDissolution(nd, v=1, result.type = "count")
form1<-tEdgeFormation(nd, v=1,  result.type = "count")

plot(form,col='green', 
     main='edge formation and dissolution rates per timestep of base')
points(diso,col='red',type='l')