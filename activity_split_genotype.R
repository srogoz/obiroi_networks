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

mixed_colonies<-c("ba16-1", "ba16-2", "ba16-3", "ba16-4", "ba16-5", "ba16-6","ba16-7","ba16-8")

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
###############
#read colony_activity

# Load .npy file function
np <- import("numpy")

colony_activity_path<-"colony_activity/activity_colonies_mm_20251028.npy"
colony_activity_split_path<-"colony_activity/split_means_20251029.npy"
colony_activity_mean<-as.numeric(np$load(colony_activity_path, allow_pickle = TRUE))
colony_activity_split_mean<-as.numeric(np$load(colony_activity_split_path, allow_pickle = TRUE))

######format

col_act_split_form<-data.frame("colonies"=mixed_colonies, "a"=colony_activity_split_mean[1:8], "b"=colony_activity_split_mean[9:16] )
###adjust treatment and replicates

treatment_replicate<-data.frame(treatment = rep(c("a", "b", "ba"), each = 8), replicate = rep(1:8, times = 3))
treatment_replicate <- 
  treatment_replicate %>%
  filter(!(treatment == "b" & replicate == 5))

treatment_replicate_split<-data.frame(treatment = rep(c("a", "b"), each = 4), replicate = rep(1:8, times = 3))
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


##############################
##############################

#%% using GLMM

library(glmmTMB)
library(tidyr)
library(dplyr)
#read as df
colony_activity_singleant<-read.csv(file = "colony_activity/activity_colonies_mm_20251104.csv")

#sort data for glmm
colony_activity_singleant_glmm<-
colony_activity_singleant <- colony_activity_singleant %>%
  mutate(across(-ant, as.numeric)) 


df_glmm <- colony_activity_singleant %>%
  pivot_longer(
    cols = -X,          # all columns except the ant ID
    names_to = "colony",
    values_to = "rsdm"
  ) %>%
  rename(ant = X) %>%  
  mutate(
    treatment = sub("\\d+.*$", "", colony)   # extract 'a', 'b', or 'ba'
  )

df_glmm <- df_glmm %>%
  mutate(
    genotype = case_when(
      treatment == "ba" ~ unlist(genotype_match[ant]),  # lookup in genotype_match
      treatment == "a" ~ "a",                            # all ants in pure a colony
      treatment == "b" ~ "b"                             # all ants in pure b colony
    )
  )

df_glmm <- df_glmm %>%
  mutate(
    puremixed = case_when(
      treatment == "ba" ~ "m",  # lookup in genotype_match
      treatment == "a" ~ "p",                            # all ants in pure a colony
      treatment == "b" ~ "p"                             # all ants in pure b colony
    )
  )

#############
#apply glmm

model <- glmmTMB(
  rsdm ~ genotype * puremixed + (1 | colony),
  data = df_glmm,
  family = Gamma()# 'ensures is positive', right skewed (link = log)
)

##################################################
##################################################
#plot all colonies and individual values
library(ggplot2)
library(viridisLite)
library(viridis)

treatment_colors<-c(
  "a" = "#05e0fc",
  "b" = "#fc0536",
  "ba" = "#db05fc"
)

output_folder <- "colony_activity"


param_plot <- ggplot(df_glmm, aes(x = colony, y = rsdm, fill = treatment)) +
 
  geom_jitter(aes(color = genotype), width = 0.1, alpha = 0.7, size = 2) +
  geom_boxplot(width = 0.3, outlier.shape = NA, color = "black", alpha = 0.4) +
  theme_minimal() +
  scale_fill_manual(values = c("a" = "#05e0fc", "b" = "#fc0536", "ba" = "#db05fc")) +
  scale_color_manual(values = c("a" = "#05e0fc", "b" = "#fc0536", "ba" = "#db05fc")) +
  labs(
    title = "colony activity (rsdm over 12hs per ant)",
    x = "colony",
    y = "rsdm [mm]"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # rotate x labels for readability
# Save each plot as PNG
ggsave(
  filename = paste0(output_folder, "/", "rsdmperant.png"),
  plot = param_plot,
  width = 6,
  height = 5,
  dpi = 300
)
