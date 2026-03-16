#separate into different timepoints
prep_delta_1<-antlevel_data_mean%>%filter(globaltime=="0-1.25h")
prep_delta_2<-antlevel_data_mean%>%filter(globaltime=="1.25-2.5h")
prep_delta_3<-antlevel_data_mean%>%filter(globaltime=="10-11.25h")
prep_delta_4<-antlevel_data_mean%>%filter(globaltime=="8.75-10h")

prep_delta_inf<-antlevel_data_mean%>%filter(globaltime=="fungal exp")

antlevel_delta_prep<-list()


  
  delta_1 <- prep_delta_1 %>%
    mutate(
      across(
        where(is.numeric),
        ~ . - prep_delta_inf[[cur_column()]]
      )
    )
  
  delta_2 <- prep_delta_2 %>%
    mutate(
      across(
        where(is.numeric),
        ~ . - prep_delta_inf[[cur_column()]]
      )
    )
  delta_3 <- prep_delta_3 %>%
    mutate(
      across(
        where(is.numeric),
        ~ . - prep_delta_inf[[cur_column()]]
      )
    )
  delta_4 <- prep_delta_4 %>%
    mutate(
      across(
        where(is.numeric),
        ~ . - prep_delta_inf[[cur_column()]]
      )
    )
  
 antlevel_delta <- bind_rows(delta_4, delta_3, delta_2, delta_1)
 antlevel_delta_mean <- antlevel_delta_mean %>%
   mutate(
     infected = ifelse(
       ant == exposed_ants[colony_name],
       "infected",
       "nestmate"
     )
   )
 
 saveRDS(antlevel_data, "network_parameters/exp2_12hs/5_mins/antlevel_data_full.rds")
 saveRDS(antlevel_data_mean, "network_parameters/exp2_12hs/5_mins/antlevel_data_mean_full.rds")
 
 antlevel_delta_mean<-readRDS("network_parameters/exp2_12hs/5_mins/antlevel_data_mean_full.rds")
 #plot antlevel deltas per colony
 param= "strength_mean"
 for (colony_id in selected_colonies_dot){
   antlevel_delta_colony<-antlevel_delta%>%filter(colony_name == colony_id)
   
   param_name<-paste0(param, "_delta")
  #  boxplot_delta<-ggplot(antlevel_delta_colony,
  #                        aes(x = ant,
  #                            y = .data[[param]],
  #                            fill = ant)) +
  #    
  #    geom_boxplot(alpha = 0.3) +
  #    
  #    geom_jitter(width = 0.2,
  #                size = 2,
  #                color = "black") +
  #    
  #    geom_hline(yintercept = 0,
  #               linetype = "dashed",
  #               color = "grey40") +
  #    
  #    labs(
  #      x = "ant",
  #      y = param_name,
  #      fill = "ant"
  #    ) +
  #    
  #    theme_minimal()
   boxplot_delta<-ggplot(antlevel_delta_colony,
          aes(x = ant,
              y = .data[[param]])) +
     
     geom_boxplot(
       aes(color = infected),
       fill = "grey85",
       alpha = 0.4
     ) +
     
     geom_jitter(width = 0.2,
                 size = 2,
                 color = "black") +
     
     geom_hline(yintercept = 0,
                linetype = "dashed",
                color = "grey40") +
     
     scale_color_manual(
       values = c(
         "infected" = "red",
         "nestmate" = "black"
       )
     ) +
     
     labs(
       x = "ant",
       y = param_name,
       color = "infection status"
     ) +
     
     theme_minimal()
   #plot hists
   ggsave(
     filename = paste0("network_parameter_plots/exp2_12hs/5mins/antlevel_mean/" , param,"_",colony_id, ".png"),
     plot = boxplot_delta,
     width = 6,
     height = 5,
     dpi = 300
   )
 }

 ## average the change of anttypes    
 
 #plot only infected
 infected_delta<-antlevel_delta_mean%>%filter(infected == "infected")
 #include genotype of infected ant
 infected_delta$infected_genotpe<-
 
 for (param in antlevel_parameters){
  
   # comparisons <- list(c("b","B"), c("b","bA"), c("b","bB"), c("B","bA"), c("B","bB"), c("bA","bB"))
   # labels <- contrasts_delta[[param]]$label
   param_name<-paste0(param, "_delta")
   boxplot_delta<-ggplot(infected_delta,
                         aes(x = treatment,
                             y = .data[[param]],
                             color = anttype,
                             fill  = treatment,
                             group = treatment
                         )) +
     geom_boxplot(aes(fill = treatment), alpha = 0.3)+
     
     geom_jitter( size = 2)+
     stat_summary(
       fun = mean,
       geom = "point",
       color = "seagreen"
       
     )+
     # geom_signif(
     #   comparisons = comparisons,
     #   annotations = labels,
     #   y_position = c(0,7000,0,8000,0,0)
     # )+
     scale_color_manual(
       values = anttypes_colors
     ) +
     scale_fill_manual(
       values = treatment_colors
     ) +
     labs(
       x = "treatments",
       y = param_name,
       color = "colony:",
       fill  = "treatment"
     ) +
     theme_minimal()
   
   ggsave(
     filename = paste0("network_parameter_plots/exp2_12hs/5mins/antlevel_mean/infected_ants_delta/" , param, ".png"),
     plot = boxplot_delta,
     width = 10,
     height = 5,
     dpi = 300
   )
   
 }
 
 #seperate boxplots for reactions
 for (param in antlevel_parameters){
   
   # comparisons <- list(c("b","B"), c("b","bA"), c("b","bB"), c("B","bA"), c("B","bB"), c("bA","bB"))
   # labels <- contrasts_delta[[param]]$label
   param_name<-paste0(param, "_delta")
   boxplot_delta<-ggplot(infected_delta,
                         aes(x = treatment,
                             y = .data[[param]],
                             fill = anttype)) +
     
     geom_boxplot(
       position = position_dodge(width = 0.8),
       alpha = 0.4
     ) +
     
     geom_jitter(
       aes(color = anttype),
       position = position_jitterdodge(
         jitter.width = 0.15,
         dodge.width = 0.8
       ),
       size = 2
     ) +
     
     stat_summary(
       aes(group = anttype),
       fun = mean,
       geom = "point",
       position = position_dodge(width = 0.8),
       color = "seagreen",
       size = 3
     ) +
     
     scale_color_manual(values = anttypes_colors) +
     scale_fill_manual(values = anttypes_colors) +
     
     labs(
       x = "treatment",
       y = param_name,
       color = "ant type \n infected ant",
       fill  = "ant type \n infected ant"
     ) +
     
     theme_minimal()
   ggsave(
     filename = paste0("network_parameter_plots/exp2_12hs/5mins/antlevel_mean/infected_ants_delta/seperate_boxplots/" , param, ".png"),
     plot = boxplot_delta,
     width = 10,
     height = 5,
     dpi = 300
   )
   
 }
 