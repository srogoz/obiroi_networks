library(EpiModel)
library(igraph)
library(intergraph)
library(networkLite)
library(network)
library(networkDynamic)
library(dplyr)
library(tidyr)
#plot prevalence
library(ggplot2)


#########load parameters and data

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

##################
#PARAMETERS
##################
#functional
interaction_limit_min<-2
interaction_limit_sec<-1
interval_time<-90 # in minutes
duration<- 50000 #simulation-duration



#plots
treatment_colors<-c(
  "a" = "#05e0fc",
  "b" = "#fc0536",
  "ba" = "#db05fc"
)

# overview parameters
parameters <-c("total_number_interactions", "total_sum_of_interactions", "strength_mean", "strength_sd", "interaction_length_mean", "interaction_length_sd", "waiting_time_mean", "waiting_time_sd", "burstiness", "density_aggr", "global_eff", "assortativity", "started_but_not_ended")
#single ant level parameters
#parameters <-c("waiting_time_ant_mean", "waiting_time_ant_sd", "burstiness_ant_sd", "burstiness_ant_mean", "inter_perant_mean", "inter_perant_sd" )


#########
#organisational

folder_path_networkpara <-"simulations/"
output_folder <- "simulations"
timestamp <- format(Sys.time(), "%d%m%Y")
############SI-simulations-function
spread_SI <-function(edge_list, spreading_prob, seed_ant, duration){
  
  #infection process, SI
  #s=0
  #i=1
  #parameters
  
  present_ants<-unique(c(edge_list$head, edge_list$tail))
  
  
  
  
  infection_hist<-matrix(0,ncol=length(present_ants), nrow=duration)
  colnames(infection_hist)<-present_ants
  #SI, no recovery infect seed ant
  
  infection_hist[1:duration, seed_ant]<-1
  
  
  #order interactions by time
  timeordered <- edge_list[order(edge_list$onset),]
  number_interactions<-nrow(timeordered)
  
  
  for (i in seq_len(number_interactions)){
    time_step = timeordered$onset[i]
    
    #respect simulation duration
    if (time_step > duration) next
    
    #interaction-partners
    partner1<-edge_list$tail[i]
    partner2<-edge_list$head[i]
    
    #collect infected ants
    infected_ants<-present_ants[infection_hist[time_step, ] == 1]
    
    
    #if both infected, skip
    
    if (partner2 %in% infected_ants && partner1 %in% infected_ants){
      next
    }
    
    # transmission from 1 to 2 with prob
    if (partner1 %in% infected_ants && !(partner2 %in% infected_ants)){
      if (runif(1) < spreading_prob) {
        # once infected, stays infected for all future times
        infection_hist[time_step:duration, partner2 ] <-1
      }
    }
    #transmission from 2 to 1
    if (partner2 %in% infected_ants && !(partner1 %in% infected_ants)){
      if(runif(1) < spreading_prob){
        infection_hist[time_step:duration, partner1 ] <-1
      }
    }
    #end simulation if everyone is infected  
    if (length(infected_ants) == length(present_ants)) break
  }
  
  prevalence<-rowSums(infection_hist)
  return(list(infection_hist = infection_hist, 
              prevalence = prevalence))
  
}

###########
#loop over all colonies
for (i in 1:length(selected_colonies)){
  
  
  colony_name<-selected_colonies[[i]]
  cat("Running simulation for colony:", colony_name, "\n")
  framerate_col<-framerate[[colony_name]]
  lower_limit_interaction<-interaction_limit_sec*framerate_col
  higher_limit_interaction<-interaction_limit_min*60*framerate_col
  
  #intervall_list<-[]
  ###################
  
  #prepare network obj, interaction limits
  network_obj <-readRDS(network_obj_list[[colony_name]])
  
  #pick time slice of network object
  
  number_intervalls<-round((12*60*60*framerate_col)/(interval_time*60*framerate_col),0)
  interval_size<-interval_time*60*framerate_col
  intervall_timepoints <- seq(
    from = 0,
    by   = interval_size,
    length.out = number_intervalls + 1
  )
  
  time_window_start<-intervall_timepoints[1]
  time_window_end<-intervall_timepoints[1+1]
  
  #limit to timeintervall
  network_obj_5<- network_obj[
    network_obj$onset < time_window_end &
      network_obj$terminus  < time_window_end,]
  
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
  
  
  #transform into network dynamic object

  edge_df <- data.frame(onset=network_obj_5$onset, terminus=network_obj_5$terminus, tail=network_obj_5$tail, head=network_obj_5$head)
  #tail and head must be numeric
 
  #edge_df$tail<-match(edge_df$tail, present_ants)
  #edge_df$head<-match(edge_df$head, present_ants)
  #edge_df$onset<-as.integer(edge_df$onset)
  
  ##CHECK NUMBER OF ONE FRAME INTERACTIONS
  bad_edges <- filter(edge_df,(onset >= terminus))
  onset_diff <- diff(bad_edges$onset)
  ##transform into undirected network object
  edge_df_corrected <- edge_df %>% filter(onset < terminus)
  
 
##########write a function for simulating spread on a network
#duration, whole length of simulation
#spreading prob: 1
#seed_ant: infected individual
#edgelist structure: onset, terminus, head, tail

edge_list<-edge_df_corrected



#result <-spread_SI(edge_list=edge_df_corrected, spreading_prob=1, seed_ant = "PP", duration=50000)
##############
###for all ants in one colony
##############

seedvariation_prevalence <- matrix(0, nrow = duration, ncol = length(present_ants))
colnames(seedvariation_prevalence) <- present_ants

# Now loop over seeds
for (ant in present_ants) {
  cat("Running simulation for seed:", ant, "\n")
  
  result <- spread_SI(edge_list = edge_df_corrected, 
                      spreading_prob = 1, 
                      seed_ant = ant, 
                      duration = duration)
  
  # Fill this ant's column with its prevalence trajectory
  seedvariation_prevalence[, ant] <- result$prevalence
}

#################

#convert to dataframe
prevalence_df<-data.frame(time = 1:length(result$prevalence), seedvariation_prevalence)

prevalence_long <- prevalence_df %>%
  pivot_longer(cols = -time, 
               names_to = "seed_ant", 
               values_to = "prevalence")
#all ants
prevalence_plot<-ggplot(prevalence_long, aes(x= time, y= prevalence, color= seed_ant))+
geom_line(alpha = 0.6, size = 0.8) +
  labs(x = "time", y = "prevalence", 
       title = paste0("SI- model: Prevalence by Seed Ant\n colony: ", colony_name)) +
  theme_minimal() +
  theme()

ggsave(filename = paste0("simulations/prevalence_" , 
                         colony_name , ".png"), plot = prevalence_plot)
}




#one ant
prevalence_plot<-ggplot(prevalence_long, aes(x= time, y= prevalence, color = seed_ant))+
  geom_point(size = 0.5, color= "purple") +
  labs(x = "Time step", y = "Number infected ants", 
       title = "SI Prevalence over time") +
  theme_minimal()



####################
#networkDynamic object TRASHE
#map ants to numbers
all_ants <- sort(unique(c(edge_df_corrected$tail, edge_df_corrected$head)))

# create a mapping
ant_map <- data.frame(
  ant_name = all_ants,
  ant_id   = 1:length(all_ants)
)
ant_map

edge_df_corrected$tail <- ant_map$ant_id[match(edge_df_corrected$tail, ant_map$ant_name)]
edge_df_corrected$head <- ant_map$ant_id[match(edge_df_corrected$head, ant_map$ant_name)]
#################################################
#################################################
base_net = data.frame(tail = edge_df_corrected$tail, head = edge_df_corrected$head)
base_net_unique<-unique(base_net)
base_net_unique <- edge_df_corrected %>% select(tail, head) %>% unique() %>%
  arrange(tail, head)
#network_obj_undirected <- network(base_net_unique, directed = FALSE, multiple= FALSE, matrix.type = "edge.list")
network_obj_undirected <- network(base_net_unique, directed = FALSE, multiple= FALSE, matrix.type = "edge.list")

##into networkDynamic object 
net_dyn <- networkDynamic(
  base.net    = network_obj_undirected,
  edge.spells = edge_df_corrected[, c("onset", "terminus", "tail", "head")],
  create.TEA = TRUE)
