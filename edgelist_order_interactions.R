#order edgelist by interactions with infected individual, so get a proxy for network measures
interval_time<-5
number_intervalls<-1
#unpper and lower interaction limit
interaction_limit_min<-2
interaction_limit_sec<-1
colony_name<-selected_colonies[[1]]
framerate_col<-framerate[[colony_name]]
present_ants<-present_ants_list[[colony_name]]
duration<-interval_time*60*framerate_col
lower_limit_interaction<-interaction_limit_sec*framerate_col
higher_limit_interaction<-interaction_limit_min*60*framerate_col

###################
edgelist_source<-"edge_lists/exp1_inf/30mins/"
network_obj_list<-get_colony_files(edgelist_source,selected_colonies )

#prepare network obj, interaction limits
network_obj <-readRDS(network_obj_list[[colony_name]])

##adjust intervallsize for other timelines

interval_size_1<-interval_time*60*framerate_col


  #start and stop
  #j, add same length if later timeintervalls
  time_window_start<-1
  #j+1 add same length if later timeintervalls
  time_window_end<-interval_size
  
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
  
  exposed_ant<-exposed_ants[[colony_name]]
  network_obj_infectious_interactions<-network_obj_5%>%filter(head == exposed_ant | tail == exposed_ant )
  #get subset of groomers interacting in that time
  interacting_ants<-unique(c(network_obj_infectious_interactions$head, network_obj_infectious_interactions$tail))
  network_obj_complete_subset_interacting_ants<-network_obj_5%>%filter(head %in% interacting_ants | tail %in% interacting_ants)
  
sorted_interactions_of_interest<-network_obj_infectious_interactions %>%
  arrange(onset) 



duration_per_partner <- network_obj_infectious_interactions %>%
  mutate(partner = ifelse(head == exposed_ant, tail, head)) %>%
  group_by(partner) %>%
  summarise(total_duration = sum(duration, na.rm = TRUE))%>%arrange(total_duration)
#unique interaction partners
first_responders <- network_obj_infectious_interactions %>%
  mutate(partner = ifelse(head == exposed_ant, tail, head)) %>%
  arrange(onset) %>%
  distinct(partner)

first_responders$order<-seq(1:length(first_responders$partner))
first_responder_order <- setNames(first_responders$order,
                                  first_responders$partner)
#plot interaction networks for time intervalls and with size and color by firstresponders
#cretae aggregated network with matching interaction-time filters
agg_network<-aggregate_from_edgelist(sorted_interactions_of_interest)
agg_network_complete_subset<-aggregate_from_edgelist(network_obj_complete_subset_interacting_ants)

network_list<-list(agg_network, agg_network_complete_subset)
finished_network_list<-list()

for (i in 1:length(network_list)){
  agg_network_iteration<-network_list[[i]]
network_aggregated_symm<-agg_network_iteration
network_aggregated_symm[lower.tri(network_aggregated_symm)]<-t(agg_network_iteration)[lower.tri(agg_network_iteration)]
#network_aggregated symm transformed into hours as unit, correcting for the framerate difference
network_aggregated_symm<-round(network_aggregated_symm/((framerate_col)),2) 
#into igraph object
graph_aggregated<-graph_from_adjacency_matrix(network_aggregated_symm, weighted = TRUE, mode = "undirected", add.colnames = NULL)

V(graph_aggregated)$strength <- strength(graph_aggregated, weights = E(graph_aggregated)$weight)
V(graph_aggregated)$first_responders<-first_responder_order[V(graph_aggregated)$name]

finished_network_list[[i]]<-graph_aggregated
}


#plot networks
library(GGally)
library(ggnetwork)
library(viridis)
library(patchwork)


colors_first_responder <- viridis(length(first_responders$partner))
color_map <- setNames(colors_first_responder, first_responders$partner)


 

graph_aggregated<-finished_network_list[[1]]
layout <- layout_with_drl(graph_aggregated)
layout<-layout*0.5

net <- ggnetwork(graph_aggregated, layout= layout)
incomplete<-ggplot(net, aes(x = x, y = y, xend = xend, yend = yend)) +
geom_edges(aes(linewidth = weight), color = "grey70", alpha = 0.6) +
geom_nodes(aes(size = strength,
                 color = first_responders)) +
  geom_nodetext(aes(label = name),
                 size = 3) +
scale_color_viridis_c(
    name = "Responder order",
    na.value = "grey85",
    option = "plasma") +
scale_size_continuous(range = c(9,20), guide = "none") +
theme_blank()


graph_aggregated<-finished_network_list[[2]]
layout <- layout_with_drl(graph_aggregated)
layout<-layout*0.5

net <- ggnetwork(graph_aggregated, layout= layout)
complete<-ggplot(net, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(aes(linewidth= weight),color = "grey70", alpha = 0.4) +
  geom_nodes(aes(size = strength,
                 color = first_responders)) +
  geom_nodetext(aes(label = name),
                size = 3) +
  scale_color_viridis_c(
    name = "Responder order",
    na.value = "grey85",
    option = "plasma") +
  scale_size_continuous(range = c(9,20), guide = "none") +
  theme_blank()

grid_plot<-(incomplete|complete) +
  plot_layout( axes = "collect" , guides = "collect")+
  theme(legend.position = "right")+plot_annotation(tag_levels = list("A", "B"))

ggsave(
  filename = paste0("network_parameter_plots/exp1_12hs/5_mins/time_mean/single_ants/ants_interaction_networks" , colony, ".png"),
  plot = grid_plot,
  width = 20,
  height = 5,
  dpi = 300
)