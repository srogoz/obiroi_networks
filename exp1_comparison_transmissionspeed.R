#try to correlate total sum of all interactions with speed of transmission

prevalence_data_merged<-readRDS("simulations/exp1_inf/10s_5min_5min/data/prevalenceALL_10s_normed.rds")
prevalence_data_merged<-prevalence_data_merged%>%mutate( infection = if_else(
  seed_ant == exposed_ants[colony_name],
  "exposed",
  "nestmate"
))

# only exposed, looks the same
prevalence_data_exposed<-prevalence_data_merged%>%filter(infection == "exposed")

prevalence_mean_treatment_exposed<-prevalence_data_exposed|>
  group_by(treatment, time) |>
  summarise(
    mean_preva = mean(prevalence_normed, na.rm = TRUE),
    n             = n(),
    sd_val   = sd(prevalence_normed, na.rm = TRUE),
    se_val   = sd_val / sqrt(n),
    .groups = "drop"
  )

#plot mean over treatments
ggplot(prevalence_mean_treatment_exposed, 
       aes(x = time, y = mean_preva, color = treatment,
           fill  = treatment)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(
      ymin = mean_preva - se_val,
      ymax = mean_preva + se_val
    ),
    alpha = 0.25,
    color = NA
  ) +
  scale_color_manual(
    values = treatment_colors
  ) +
  scale_fill_manual(
    values = treatment_colors
  ) +
  xlim(0,15000)+
  labs(x = "time [frame = 0.1s]", y = " mean prevalence",
       color = "treatment", fill = "treatment",
       title = paste0("LONG INTERACTIONS (10s-3min) only exposed ants, \n p = ")) +
  theme_minimal()

#######################
#all ants
#slope over colony mean 
transmission_slopes<-readRDS("simulations/exp1_inf/10s_5min_5min/data/slopes.rds")

#network parameters right before exposure
network_parameters_before<-read.csv("network_parameters/exp1_beforefungus/5_min/121_135/10s_3min_20minx1/ALL_network_parameters_allintervalls_23062026.csv")

comparison_colonies<-cbind(network_parameters_before, transmision_slopes)

#testing if the transmission speed is affected by treatment and total sum of interactions

gpval<-summary(emmeans(slope_fit_B_glmm, pairwise ~ infection, type = "response")$contrasts)$p.value


model_sit<-glmmTMB(
   slope ~ treatment*total_sum_of_interactions+ (1|colony_name) ,
  data = comparison_colonies,
  family = beta_family()
  #family = Gamma(link = "log")
)

sim <- DHARMa::simulateResiduals(model_sit)
plot(sim)

drop1(model_sit)
gpval<-summary(emmeans(slope_fit_B_glmm, pairwise ~ infection, type = "response")$contrasts)$p.value
glmm_b[[param]]<-gpval
emmeans(model_gt, pairwise ~ genotype | puremixed)