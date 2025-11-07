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

selected_colonies<-c("a16-1", "a16-2", "a16-3", "a16-4", "a16-5", "a16-6",
                        "b16-1", "b16-2", "b16-3", "b16-4", "b16-6",
                        "ba16-1", "ba16-2", "ba16-3", "ba16-4", "ba16-5", "ba16-6",
                        "a16-7","a16-8","b16-7","b16-8","ba16-7","ba16-8"
                        )


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
                      "ba16-8" = "3")
###############
#read colony_activity

# Load .npy file function
np <- import("numpy")

colony_activity_path<-"colony_activity/activity_mean_mm_20251028.npy"
colony_activity_mean<-as.numeric(np$load(colony_activity_path, allow_pickle = TRUE))

######format

###adjust treatment and replicates

treatment_replicate<-data.frame(treatment = rep(c("a", "b", "ba"), each = 8), replicate = rep(1:8, times = 3))
treatment_replicate <- 
  treatment_replicate %>%
  filter(!(treatment == "b" & replicate == 5))

activity_complete<-data_frame(selected_colonies, colony_activity_mean, treatment_replicate)


###########################
#########test for normality

library(ggplot2)
library(viridisLite)
library(viridis)

treatment_colors<-c(
  "a" = "#05e0fc",
  "b" = "#fc0536",
  "ba" = "#db05fc"
)

# Compute Shapiro–Wilk test per treatment for this parameter
res <- activity_complete %>%
  group_by(treatment) %>%
  summarise(
    shapiro_p = shapiro.test(colony_activity_mean)$p.value,
    .groups = "drop"
  )

#######not normal
##Kruskal-Wallis

  formula <- as.formula(paste("colony_activity_mean", "~ treatment"))
  
  res<-kruskal.test(formula, data = activity_complete)

#################
#dunns test or pairwise wilcoxon

library(FSA)
library(rstatix)
library(ggpubr)    

##################
stat_test <- activity_complete %>%
  dunn_test(colony_activity_mean ~ treatment, p.adjust.method = "holm")
############
stat_test<-pairwise.wilcox.test(
  x = activity_complete$colony_activity_mean,
  g = activity_complete$treatment,
  p.adjust.method = "BH"  # Benjamini–Hochberg correction
)
  
############
#activity-plot
  output_folder <- "colony_activity"
  
    
    param_plot<-ggplot(activity_complete, 
                       aes(x = treatment, y = colony_activity_mean, fill = treatment)) +
      geom_violin(trim = FALSE, alpha = 0.4) +  # “bell shape”
      geom_boxplot(width = 0.15, outlier.shape = NA, color = "black") +
      geom_jitter(width = 0.1, alpha = 0.7) +   # individual points
      geom_text(aes(label = replicate), vjust = -1, size = 3.5, color = "black") +
      stat_pvalue_manual(stat_test, label = "p.adj.signif", 
                         y.position = c(1.1, 1.2, 1.3) * max(activity_complete$colony_activity_mean))+
      theme_minimal() +
      # Apply manual colors
      scale_fill_manual(values = treatment_colors) +
      labs(
        title = paste0("colony-activity  ", "p = ", round(res$p.value,2)),
        x = "Treatment",
        y = "r.s.d.m.[mm]") 
    #  scale_fill_viridis_d(option = "D")
    
    # Save each plot as PNG
    ggsave(
      filename = paste0(output_folder, "/", "colony_activity_violin_boxplot_dunns.png"),
      plot = param_plot,
      width = 6,
      height = 5,
      dpi = 300
    )
  

