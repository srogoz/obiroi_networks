###############
#plot interaction histogramms and save values of replaced values
#for reloading function 
#count how many interaction and percentages of frames ommited
interaction_control <- setNames(
  data.frame(matrix(NA, nrow = length(selected_colonies), ncol = 6)),
  c("sum_interactions", "sum_frames", "frames_upper", 
    "frames_lower", "interactions_upper", "interactions_lower")
)

rownames(interaction_control) <- selected_colonies


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
  print(paste0(colony_name, "done"))
}

write.csv(interaction_control,
          file =  "network_parameter/antlevel_1_truncated_2min/DATA_interaction_control_1s_2min.csv",
          row.names = FALSE)

mean(interaction_control$sum_frames)
sd(interaction_control$sum_frames)
mean(interaction_control$sum_interactions)
sd(interaction_control$sum_interactions)
