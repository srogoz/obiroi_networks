#antlevel data plotting
################################
library(tidyr)

antlevel_files<-c("network_parameters/exp1_12hs/5_mins/1_15/ANTLEVEL_allintervalls_strength_mean26022026.csv",
                  "network_parameters/exp1_12hs/5_mins/16_30/ANTLEVEL_allintervalls_strength_mean26022026.csv",
                  "network_parameters/exp1_12hs/5_mins/106_120/ANTLEVEL_allintervalls_strength_mean27022026.csv",
                    "network_parameters/exp1_12hs/5_mins/120_135/ANTLEVEL_allintervalls_strength_mean27022026.csv",
                  "network_parameters/exp1_inf/5mins/ANTLEVEL_allintervalls_strength_mean27022026.csv" #infection
)

antlevel_data_unmerged<-lapply(antlevel_files, read.csv)

#timepoint marker
names(antlevel_data_unmerged[[1]])<-paste0(names(antlevel_data_unmerged[[1]]), "_0-1.25h")
names(antlevel_data_unmerged[[2]])<-paste0(names(antlevel_data_unmerged[[2]]), "_1.25-2.5h")
names(antlevel_data_unmerged[[3]])<-paste0(names(antlevel_data_unmerged[[3]]), "_8.75-10h")
names(antlevel_data_unmerged[[4]])<-paste0(names(antlevel_data_unmerged[[4]]), "_10-11.25h")
names(antlevel_data_unmerged[[5]])<-paste0(names(antlevel_data_unmerged[[5]]), "_fungal exp")




#adding ant level
antlevel_data_unmerged<-lapply(antlevel_data_unmerged, function(df){
  df$ant<-expected_ants
  return(df)
})


antlevel_data_merged<-bind_rows(antlevel_data_unmerged)


antlevel_data<- antlevel_data_merged%>%
  pivot_longer(
    cols = -ant,
    names_to = "colony_name",
    values_to = "strength_mean"
  )%>%
  separate(colony_name,
           into = c("colony_name", "time_interval", "globaltime"),
           sep = "_",
           remove = TRUE)%>%
  mutate(time_interval = as.numeric(time_interval))%>%
  mutate(
    treatment = sub("\\d+.*$", "", colony_name)   # extract 'a', 'b', or 'ba'
  )%>%
  mutate(
    genotype = case_when(
      treatment == "ba" ~ unlist(genotype_match[ant]),  # lookup in genotype_match
      treatment == "a" ~ "a",                            # all ants in pure a colony
      treatment == "b" ~ "b"                             # all ants in pure b colony
    )
  )%>%
  mutate(
    puremixed = case_when(
      treatment == "ba" ~ "m",  # lookup in genotype_match
      treatment == "a" ~ "p",                            # all ants in pure a colony
      treatment == "b" ~ "p"                             # all ants in pure b colony
    )
  )

#filter for early timepoint in infection
antlevel_data_30min<-antlevel_data%>%filter(time_interval<= 6)

#mean over all 30 mins per ant
antlevel_data_mean<-antlevel_data_30min%>%group_by(colony_name, globaltime,ant, puremixed, genotype) %>%
  summarise(
    across(where(is.numeric), mean, na.rm = TRUE),
    .groups = "drop"
  )%>%select(-time_interval)

#substract deltas for dingle ant leve
antlevel_1<-antlevel_data_mean%>%filter(globaltime=="0-1.25h")
antlevel_2<-antlevel_data_mean%>%filter(globaltime=="1.25-2.5h")
antlevel_3<-antlevel_data_mean%>%filter(globaltime=="10-11.25h")
antlevel_4<-antlevel_data_mean%>%filter(globaltime=="8.75-10h")

antlevel_inf<-antlevel_data_mean%>%filter(globaltime=="fungal exp")


antdelta_1 <- antlevel_1 %>%
  mutate(
    across(
      where(is.numeric),
      ~ . - antlevel_inf[[cur_column()]]))

antdelta_2 <- antlevel_2 %>%
  mutate(
    across(
      where(is.numeric),
      ~ . - antlevel_inf[[cur_column()]]))

antdelta_3 <- antlevel_3 %>%
  mutate(
    across(
      where(is.numeric),
      ~ . - antlevel_inf[[cur_column()]]))

antdelta_4 <- antlevel_4 %>%
  mutate(
    across(
      where(is.numeric),
      ~ . - antlevel_inf[[cur_column()]]))


#dim= colony 24 x ants 16 x datapoints 4 ~1536
antlevel_delta <- bind_rows(antdelta_4, antdelta_3, antdelta_2, antdelta_1)

##plot  deltas for ants in different colonies
#mark exposed ant and genotype
infection_ant_colors<-c("exposed" = "springgreen", "nestmate" = "blue")
param<-"strength_mean"
for (colony in selected_colonies_dot){
  
single_colony<-antlevel_delta|>filter(colony_name == colony)
single_colony$color <- ifelse(
  single_colony$ant %in% exposed_ants_plot[[colony]],
  "exposed",
  "nestmate"
)

boxplot_delta<-ggplot(single_colony,
                      aes(x = ant,
                          y = .data[[param]],
                          color = color,
                          fill  = genotype,
                          group = ant
                      )) +
  geom_boxplot( alpha = 0.3)+
  geom_point( size = 2)+
  # geom_signif(
  #   comparisons = comparisons,
  #   annotations = labels,
  #   y_position = c( 0.03, 0.04, 0.05)
  # )+
  scale_color_manual(
    values = infection_ant_colors
  ) +
  scale_fill_manual(
    values = anttypes_colors
  ) +
  labs(
    title = colony,
    x = "ants",
    y = paste0(param,"_delta" ),
    color = "colony:",
    fill  = "genotype"
  ) +
  theme_minimal()

ggsave(
  filename = paste0("network_parameter_plots/exp1_12hs/5_mins/time_mean/single_ants/delta_exposed_marked/" , colony, ".png"),
  plot = boxplot_delta,
  width = 10,
  height = 5,
  dpi = 300
)
}

##plot direct!! parameter values for ants in different colonies
#mark exposed ant and genotype
infection_ant_colors<-c("exposed" = "springgreen", "nestmate" = "blue")
param<-"strength_mean"


for (colony in selected_colonies_dot){
  
  single_colony<-antlevel_data_mean|>filter(colony_name == colony)
  #position of infected ant
  exposed_ant<-exposed_ants_plot[colony]
  position_exposed<-which(levels(factor(single_colony$ant)) == exposed_ant)
  single_colony$color <- ifelse(
    single_colony$ant %in% exposed_ants_plot[[colony]],
    "exposed",
    "nestmate"
  )
  single_colony$infection <- ifelse(
    single_colony$globaltime == "fungal exp",
    "after",
    "before"
  )
  
  if (single_colony$puremixed[1]=="m"){
  boxplot_delta<-ggplot(single_colony,
                        aes(x = ant,
                            y = .data[[param]],
                            color = infection,
                            group = ant
                        )) +
    geom_rect(
      xmin = -Inf,
      xmax = 8.5,     # halfway point (depends on number of ants)
      ymin = -Inf,
      ymax = Inf,
      fill = "lightpink",
      alpha = 0.2
    ) +
    geom_rect(
      xmin = 8.5,
      xmax = Inf,     # halfway point (depends on number of ants)
      ymin = -Inf,
      ymax = Inf,
      fill = "cyan",
      alpha = 0.2
    ) +
    geom_rect(
      xmin = (position_exposed-0.5),
      xmax = (position_exposed+0.5),    
      ymin = -Inf,
      ymax = Inf,
      fill = "springgreen",
      alpha = 0.2
    )+
    geom_boxplot( alpha = 0.3)+
    geom_point( size = 2)+
    # geom_signif(
    #   comparisons = comparisons,
    #   annotations = labels,
    #   y_position = c( 0.03, 0.04, 0.05)
    # )+
    scale_color_manual(
      values = infection_status
    ) +
    # scale_fill_manual(
    #   values = anttypes_colors
    # ) +
    
    labs(
      title = colony,
      x = "ants",
      y = param,
      color = "colony:",
      #fill  = "genotype"
    ) +
    theme_minimal()
  
 
  ggsave(
    filename = paste0("network_parameter_plots/exp1_12hs/5_mins/time_mean/single_ants/strength_pre_post/" , colony, ".png"),
    plot = boxplot_delta,
    width = 10,
    height = 5,
    dpi = 300
  )
  }
  else{
    boxplot_delta<-ggplot(single_colony,
                          aes(x = ant,
                              y = .data[[param]],
                              color = infection,
                              group = ant
                          )) +
     
      geom_rect(
        xmin = (position_exposed-0.5),
        xmax = (position_exposed+0.5),    
        ymin = -Inf,
        ymax = Inf,
        fill = "springgreen",
        alpha = 0.2
      )+
      geom_boxplot( alpha = 0.3)+
      geom_point( size = 2)+
      # geom_signif(
      #   comparisons = comparisons,
      #   annotations = labels,
      #   y_position = c( 0.03, 0.04, 0.05)
      # )+
      scale_color_manual(
        values = infection_status
      ) +
      # scale_fill_manual(
      #   values = anttypes_colors
      # ) +
      
      labs(
        title = colony,
        x = "ants",
        y =param,
        color = "colony:",
        #fill  = "genotype"
      ) +
      theme_minimal()
    
    
    ggsave(
      filename = paste0("network_parameter_plots/exp1_12hs/5_mins/time_mean/single_ants/strength_pre_post/" , colony, ".png"),
      plot = boxplot_delta,
      width = 10,
      height = 5,
      dpi = 300
    )
  }
}



#########################
#plot values for all ants
#########################
library(ggrepel)
global_times_colors_a<-c("0-1.25h" = "cyan", 
                         "1.25-2.5h" = "slateblue2",
                         "8.75-10h" = "blue",
                         "10-11.25h" = "dodgerblue1",
                         "fungal exp"= "springgreen1")


output_folder <- "network_parameter_plots/antlevel_1_truncated_2min/"

ants_a <- antlevel_data%>%filter(treatment == "b")
halflist_a<-c("b16.1", "b16.2", "b16.3", "b16.4")
ants_a <- antlevel_data%>%filter(colony %in% halflist_a)
ants_a<-ants_a%>%filter(time_interval !=15)

ants_a <- ants_a %>%
  arrange(colony, ant) %>%
  mutate(colony_ant = interaction(ant,colony, sep = "_"))
ants_a<-ants_a%>%mutate(infected = colony_ant %in% exposed_ants_string)

ggplot(ants_a,
       aes(x = colony_ant,
           y = strength_mean,
           color = globaltime)) +
  geom_point(
    aes(fill = globaltime),
    shape = 21,             # circle with outline
    color = ifelse(ants_a$infected, "red","white"),  # outline
    size = 3,
    stroke = 1.5,           # thickness of outline
    alpha = 0.7
  ) +
  theme_minimal() +
  scale_color_manual(values = global_times_colors_a) +
  labs(
    x = "Ant (within colony)",
    y = "strength_mean [s]"
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggplot(ants_a, aes(x = colony, y = strength_mean, fill = treatment)) +
  
  geom_jitter(aes(color = genotype), width = 0.1, alpha = 0.7, size = 2) +
  geom_boxplot(width = 0.3, outlier.shape = NA, color = "black", alpha = 0.4) +
  theme_minimal() +
  scale_fill_manual(values = global_times_colors_a) +
  scale_color_manual(values = global_times_colors_a) +
  labs(
    title = " ",
    x = "colony",
    y = "strength_mean [s]"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # rotate x labels for readability

ggplot(ants_a, aes(x = colony, y = strength_mean, color= globaltime, group = ant)) +
  geom_point(
    aes(color = globaltime),
    position = position_jitterdodge(
      jitter.width = 0.1,
      dodge.width = 1
    ),
    size = 1,
    alpha = 0.7
  )+
  theme_minimal() +
  scale_fill_manual(values = global_times_colors_a) +
  scale_color_manual(values = global_times_colors_a) +
  labs(
    x = "Colony",
    y = "strength_mean [s]"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# Save each plot as PNG
ggsave(
  filename = paste0(output_folder, "/", "ANTLEVEL_meaninteractionlength.png"),
  plot = param_plot,
  width = 6,
  height = 5,
  dpi = 300
)