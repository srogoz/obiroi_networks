#open and plot graph same node size, no number/ numbers less aggregated and connected
library(intergraph) # convert between alternative network objects
library(sna)        # alternative network objects
library(RColorBrewer)
library(magrittr)
library(dplyr)
library(igraph)
library(GGally)
library(ggrepel)
library(ggplot2)
library(viridisLite)
library(viridis)
library(igraph)
library(ggnetwork)

#upper and lower interaction limit
interaction_limit_min<-2
interaction_limit_sec<-10
interval_time<-8 #12hs [minutes]
#get network parameters 
network_obj_list<-get_colony_files(source_folder_edgelist, selected_colonies)
#######################################
#select colony to plot
  colony_name<-selected_colonies[[17]]
  framerate_col<-framerate[[colony_name]]
  present_ants<-present_ants_list[[colony_name]]
  duration<-interval_time*60*framerate_col
  lower_limit_interaction<-interaction_limit_sec*framerate_col
  higher_limit_interaction<-interaction_limit_min*60*framerate_col
  
  ###################
  
  #prepare network obj, interaction limits
  network_obj <-readRDS(network_obj_list[[colony_name]])
  
  ##adjust intervallsize for other timelines
  
  
    #start and stop
    #j, add same length if later timeintervalls
    time_window_start<-1
    #j+1 add same length if later timeintervalls
    time_window_end<-duration
    
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
    
    
treatment<-"bB"
ig<-graph_aggregated
ig<-ig%>%set_vertex_attr("strength",value=igraph::strength(ig))
new_layout<-layout_with_fr(ig)  # try to adjust to the best layout for each type of appendage
#new_layout<-layout_with_kk(ig)
V(ig)$strength <- as.numeric(V(ig)$strength)

V(ig)$strength <- as.numeric(V(ig)$strength)
V(ig)$name <- expected_ants
#add anttype, genotype and age
V(ig)$anttype <- sapply(
  V(ig)$name,
  function(x) anttype_match[[treatment]][[x]]
)

E(ig)$weight <- as.numeric(E(ig)$weight)
E(ig)$scaled_weight <- scales::rescale(E(ig)$weight, to = c(0.2, 3))

#same size nodes
anttypes_colors<-list( "b" ="magenta3",
                       "B" ="seagreen2" )

p<-ggnet2(
  ig,
  size = 9,                    # all nodes same size
  color = "anttype",           # use vertex attribute
  palette = anttypes_colors,    # named color vector
  edge.size = "scaled_weight",        # edge thickness from weight
  label = FALSE   # no node names
  ) +
  guides(           #to remove legend
    color = "none",
    size = "none"
  )

p <- ggnet2(
  ig,
  mode = new_layout,
  label = FALSE,
 
  size = 2,
  edge.size = "scaled_weight",
  edge.alpha = 0.6,
  edge.label = "weight",
  edge.label.size = 1
) +
  geom_text_repel(label = V(ig)$vertex.names, size = 3) +
  scale_color_viridis_d(option = "D", name = "Node strength [hs]") +
  scale_size_continuous(name = "Node strength [hs]") +
  labs(title = "Aggregated interaction network for a16-4 over 12hs") +
  theme_minimal()

p<-ggnet2(ig,
          mode=new_layout,
          label=FALSE,
          color="strength",
          color.legend="strength, interactiontime[hs]",
          #color.palette=viridis_pal(option = "D")(50),
          size="strength", 
          size.min = 8,
          edge.size = "weight", 
          edge.label = "weight",
          edge.size.min = 0.1,
          edge.alpha=0.6,
          edge.label.size = 3) +
  labs(title = "Aggregated interaction network for a16-4 over 12hs")+
  scale_color_viridis_d(option = "D", name = "Node strength [hs]") +
  # scale_size_continuous(name = "Node strength [hs]") +
  guides(color= guide_legend("strength [hs]")) +
  geom_text_repel(label = V(net)$vertex.names, size = 3)

ggnet2(
  ig,
  mode = new_layout,
  label = FALSE,
  color = "strength",    # numeric
  size = "strength", 
  edge.size = "weight", 
  edge.label = "weight",
  edge.size.min = 0.1,
  edge.alpha = 0.6,
  edge.label.size = 3
) +
  geom_text_repel(label = V(ig)$vertex.names) +
  scale_color_viridis_c(option = "D", name = "Node strength [hs]") +  # continuous bar
  scale_size_continuous(name = "Node strength [hs]") +
  labs(title = "Aggregated interaction network for a16-4 over 12hs")


svg(filename=paste("minimal_network_plots/",networks[i],"minimal network plot with colored modules.svg"),width=10,height=7)
print(p)
dev.off()