mean_overtimes<-before_after_infection_merged%>%group_by(colony_name, globaltime,infection) %>%
  summarise(
    across(where(is.numeric), mean, na.rm = TRUE),
    .groups = "drop"
  )
#mean over before and after infection
prep_delta_1<-before_after_infection_mean%>%filter(globaltime=="0-1.25h")
prep_delta_2<-before_after_infection_mean%>%filter(globaltime=="1.25-2.5h")
prep_delta_3<-before_after_infection_mean%>%filter(globaltime=="10-11.25h")
prep_delta_4<-before_after_infection_mean%>%filter(globaltime=="8.75-10h")

prep_delta_inf<-before_after_infection_mean%>%filter(globaltime=="fungal exp")


delta_1 <- prep_delta_1 %>%
  mutate(
    across(
      where(is.numeric),
      ~ . - prep_delta_inf[[cur_column()]]))

delta_2 <- prep_delta_2 %>%
  mutate(
    across(
      where(is.numeric),
      ~ . - prep_delta_inf[[cur_column()]]))

delta_3 <- prep_delta_3 %>%
  mutate(
    across(
      where(is.numeric),
      ~ . - prep_delta_inf[[cur_column()]]))

delta_4 <- prep_delta_4 %>%
  mutate(
    across(
      where(is.numeric),
      ~ . - prep_delta_inf[[cur_column()]]))



before_after_infection_delta <- bind_rows(delta_4, delta_3, delta_2, delta_1)%>%
  select(-c(number_ants, time_interval, infection))

#correct global efficiency for exp1
before_after_infection_delta <- before_after_infection_delta%>%mutate(global_eff_corr = global_eff/10000)

#exp1
strength_interval<-c(800, 900, 1000)
assortativity<-c( 0.05, 0.1, 0.15)
global_eff_corr<-c( 0.05, 0.1, 0.15)
density_aggr<-c( 0.03, 0.04, 0.05)

#exp2


for (param in parameters){
  png(file= paste0("network_parameter_plots/exp2_12hs/5mins/time_mean/deltas/", param, ".png"),
      width=600, height=350)
  hist(before_after_infection_delta[[param]], xlab = param , col = "purple", main="")
  dev.off()
}
##plot deltas in dependence of treatment
#exp1
parameters_delta<-c("strength_mean", "strength_sd", "density_aggr","assortativity" , "global_eff_corr", "interaction_length_mean")
#exp2
parameters_delta<-c("strength_mean", "strength_sd", "density_aggr","assortativity" , "global_eff_corr", "interaction_length_mean", "global_clustering", "total_number_of_interactons", "total_sum_of_interactions")

for (param in parameters){
 param<-"strength_mean"
  comparisons <- list(c("b","B"), c("b","bA"), c("b","bB"), c("B","bA"), c("B","bB"), c("bA","bB"))
  labels <- contrasts_delta[[param]]$label
  param_name<-paste0(param, "_delta")
  boxplot_delta<-ggplot(before_after_infection_delta,
                        aes(x = treatment,
                            y = .data[[param]],
                            color = colony_name,
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
    geom_signif(
      comparisons = comparisons,
      annotations = labels,
     y_position = c( 1000, 1050, 1100, 1150, 1200, 1250)
    )+
    scale_color_manual(
      values = colony_colors
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
    filename = paste0("network_parameter_plots/exp2_12hs/5mins/time_mean/deltas/" , param, ".png"),
    plot = boxplot_delta,
    width = 10,
    height = 5,
    dpi = 300
  )
  
}
##glmm for parameters for exp2
##save in file for plots
contrasts_delta<-list()

#strength
hist(before_after_infection_delta$strength_mean)
model_delta_s<-lm(strength_mean~treatment, 
                  data = before_after_infection_delta,
                  
)

sim <- DHARMa::simulateResiduals(model_delta_s)
plot(sim)

drop1(model_delta_s, test = "Chisq")
contrasts_delta[["strength_mean"]]<-as.data.frame(summary(emmeans(model_delta_s, pairwise ~ treatment, type = "response"))$contrasts)
#add significance labels
contrasts_delta[["strength_mean"]]$label<- ifelse(contrasts_delta[["strength_mean"]]$p.value < 0.001, "***",
                                                  ifelse(contrasts_delta[["strength_mean"]]$p.value < 0.01, "**",
                                                         ifelse(contrasts_delta[["strength_mean"]]$p.value < 0.05, "*", "ns")))


#strength
hist(before_after_infection_delta$strength_mean)
hist((before_after_infection_delta%>%filter(treatment == "bB"))$strength_mean)
model_delta_s<-lm(strength_mean~treatment, 
                  data = before_after_infection_delta
)

sim <- DHARMa::simulateResiduals(model_delta_s)
plot(sim)

model_delta_s_glm<-glmmTMB(strength_mean~treatment, 
                           data = before_after_infection_delta,
                           family = gaussian()
)

sim <- DHARMa::simulateResiduals(model_delta_s_glm)
plot(sim)
drop1(model_delta_s_glm, test = "Chisq")
summary(emmeans(model_delta_s_glm, pairwise ~ treatment, type = "response"))$contrasts

#subsets on global time
subset_means<- before_after_infection_delta%>%filter(globaltime=="8.75-10h") #1.25-2.5h, "10-11.25h" ,"8.75-10h"



model_delta_s_subset<-glmmTMB(strength_mean~treatment, 
                              data = subset_means,
                              family = gaussian()
)

sim <- DHARMa::simulateResiduals(model_delta_s_subset)
plot(sim)

drop1(model_delta_s_subset, test = "Chisq")
summary(emmeans(model_delta_s_subset, pairwise ~ treatment, type = "response"))$contrasts

drop1(model_delta_s, test = "Chisq")
contrasts_delta[["strength_mean"]]<-as.data.frame(summary(emmeans(model_delta_s, pairwise ~ treatment, type = "response"))$contrasts)
#add significance labels
contrasts_delta[["strength_mean"]]$label<- ifelse(contrasts_delta[["strength_mean"]]$p.value < 0.001, "***",
                                                  ifelse(contrasts_delta[["strength_mean"]]$p.value < 0.01, "**",
                                                         ifelse(contrasts_delta[["strength_mean"]]$p.value < 0.05, "*", "ns")))
#density_aggr
hist(before_after_infection_delta$density_aggr)
model_delta_d<-glmmTMB(density_aggr~treatment, 
                  data = before_after_infection_delta,
                  family = gaussian()
)

sim <- DHARMa::simulateResiduals(model_delta_d)
plot(sim)

drop1(model_delta_d, test = "Chisq")
emmeans(model_delta_d, pairwise ~ treatment, type = "response")
contrasts_delta[["density_aggr"]]<-as.data.frame(summary(emmeans(model_delta_d, pairwise ~ treatment, type = "response"))$contrasts)
#add significance labels
contrasts_delta[["density_aggr"]]$label<- ifelse(contrasts_delta[["density_aggr"]]$p.value < 0.001, "***",
                                                 ifelse(contrasts_delta[["density_aggr"]]$p.value < 0.01, "**",
                                                        ifelse(contrasts_delta[["density_aggr"]]$p.value < 0.05, "*", "ns")))


#assortativity
hist(before_after_infection_delta$assortativity)
model_delta_a<-glmmTMB(assortativity~treatment, 
                  data = before_after_infection_delta,
                  #dispformula = ~1,
                  family = gaussian()
)

sim <- DHARMa::simulateResiduals(model_delta_a)
plot(sim)

drop1(model_delta_a, test = "Chisq")
emmeans(model_delta_a, pairwise ~ treatment, type = "response")
contrasts_delta[["assortativity"]]<-as.data.frame(summary(emmeans(model_delta_a, pairwise ~ treatment, type = "response"))$contrasts)
#add significance labels
contrasts_delta[["assortativity"]]$label<- ifelse(contrasts_delta[["assortativity"]]$p.value < 0.001, "***",
                                                  ifelse(contrasts_delta[["assortativity"]]$p.value < 0.01, "**",
                                                         ifelse(contrasts_delta[["assortativity"]]$p.value < 0.05, "*", "ns")))


###strength_sd
hist(before_after_infection_delta$strength_sd)
model_delta_ssd<-glmmTMB(strength_sd~treatment, 
                    data = before_after_infection_delta,
                    #dispformula = ~1,
                    family = gaussian()
)

sim <- DHARMa::simulateResiduals(model_delta_ssd)
plot(sim)

drop1(model_delta_ssd, test = "Chisq")
emmeans(model_delta_ssd, pairwise ~ treatment, type = "response")
contrasts_delta[["strength_sd"]]<-as.data.frame(summary(emmeans(model_delta_ssd, pairwise ~ treatment, type = "response"))$contrasts)
#add significance labels
contrasts_delta[["strength_sd"]]$label<- ifelse(contrasts_delta[["strength_sd"]]$p.value < 0.001, "***",
                                                ifelse(contrasts_delta[["strength_sd"]]$p.value < 0.01, "**",
                                                       ifelse(contrasts_delta[["strength_sd"]]$p.value < 0.05, "*", "ns")))

#burstiness
hist(before_after_infection_delta$burstiness)
model_delta_b<-glmmTMB(burstiness~treatment, 
                  data = before_after_infection_delta,
                  #dispformula = ~1,
                  family = gaussian()
)

sim <- DHARMa::simulateResiduals(model_delta_b)
plot(sim)

drop1(model_delta_b, test = "Chisq")

#global effectiveness
hist(before_after_infection_delta$global_eff)
model_delta_ge<-glmmTMB(global_eff~treatment, 
                   data = before_after_infection_delta,
                   #dispformula = ~1,
                   family = gaussian()
)

sim <- DHARMa::simulateResiduals(model_delta_ge)
plot(sim)

drop1(model_delta_ge, test = "Chisq")
emmeans(model_delta_ge, pairwise ~ treatment, type = "response")
contrasts_delta[["global_eff_corr"]]<-as.data.frame(summary(emmeans(model_delta_ge, pairwise ~ treatment, type = "response"))$contrasts)
#add significance labels
contrasts_delta[["global_eff_corr"]]$label<- ifelse(contrasts_delta[["global_eff_corr"]]$p.value < 0.001, "***",
                                                    ifelse(contrasts_delta[["global_eff_corr"]]$p.value < 0.01, "**",
                                                           ifelse(contrasts_delta[["global_eff_corr"]]$p.value < 0.05, "*", "ns")))

##waitingtime_mean
hist(before_after_infection_delta$waiting_time_mean)
model_delta_wt<-glmmTMB(waiting_time_mean~treatment, 
                   data = before_after_infection_delta,
                   dispformula = ~1,
                   family = t_family()
)

sim <- DHARMa::simulateResiduals(model_delta_wt)
plot(sim)

drop1(model_delta_wt, test = "Chisq")
#interactionlength mean
hist(before_after_infection_delta$interaction_length_mean)
model_delta_im<-glmmTMB(interaction_length_mean~treatment, 
                   data = before_after_infection_delta,
                   #dispformula = ~1,
                   family = gaussian()
)

sim <- DHARMa::simulateResiduals(model_delta_im)
plot(sim)

drop1(model_delta_im, test = "Chisq")
emmeans(model_delta_im, pairwise ~ treatment, type = "response")
contrasts_delta[["interaction_length_mean"]]<-as.data.frame(summary(emmeans(model_delta_im, pairwise ~ treatment, type = "response"))$contrasts)
#add significance labels
contrasts_delta[["interaction_length_mean"]]$label<- ifelse(contrasts_delta[["interaction_length_mean"]]$p.value < 0.001, "***",
                                                            ifelse(contrasts_delta[["interaction_length_mean"]]$p.value < 0.01, "**",
                                                                   ifelse(contrasts_delta[["interaction_length_mean"]]$p.value < 0.05, "*", "ns")))

#global clustering
hist(before_after_infection_delta$clustering_global)
model_delta_gc<-glmmTMB(clustering_global~treatment, 
                    data = before_after_infection_delta,
                    dispformula = ~treatment,
                    family = gaussian()
)

sim <- DHARMa::simulateResiduals(model_delta_gc)
plot(sim)

drop1(model_delta_gc, test = "Chisq")
emmeans(model_delta_gc, pairwise ~ treatment, type = "response")
contrasts_delta[["clustering_global"]]<-as.data.frame(summary(emmeans(model_delta_gc, pairwise ~ treatment, type = "response"))$contrasts)
#add significance labels
contrasts_delta[["clustering_global"]]$label<- ifelse(contrasts_delta[["clustering_global"]]$p.value < 0.001, "***",
                                                            ifelse(contrasts_delta[["clustering_global"]]$p.value < 0.01, "**",
                                                                   ifelse(contrasts_delta[["clustering_global"]]$p.value < 0.05, "*", "ns")))


#mean distance
#interactionlength mean
hist(before_after_infection_delta$mean_distance)
model_delta_md<-glmmTMB(mean_distance~treatment, 
                        data = before_after_infection_delta,
                        dispformula = ~treatment,
                        family = gaussian()
)

sim <- DHARMa::simulateResiduals(model_delta_md)
plot(sim)

drop1(model_delta_md, test = "Chisq")
#total number interactions
hist(before_after_infection_delta$mean_distance)
model_delta_md<-glmmTMB(mean_distance~treatment, 
                        data = before_after_infection_delta,
                        dispformula = ~treatment,
                        family = gaussian()
)

sim <- DHARMa::simulateResiduals(model_delta_md)
plot(sim)

drop1(model_delta_md, test = "Chisq")
emmeans(model_delta_im, pairwise ~ treatment, type = "response")
contrasts_delta[["interaction_length_mean"]]<-as.data.frame(summary(emmeans(model_delta_im, pairwise ~ treatment, type = "response"))$contrasts)
#add significance labels
contrasts_delta[["interaction_length_mean"]]$label<- ifelse(contrasts_delta[["interaction_length_mean"]]$p.value < 0.001, "***",
                                                            ifelse(contrasts_delta[["interaction_length_mean"]]$p.value < 0.01, "**",
                                                                   ifelse(contrasts_delta[["interaction_length_mean"]]$p.value < 0.05, "*", "ns")))

#total interactionlsum

hist(before_after_infection_delta$total_sum_of_interactions)
model_delta_ts<-glmmTMB(total_sum_of_interactions~treatment, 
                        data = before_after_infection_delta,
                        #dispformula = ~treatment,
                        family = gaussian()
)

sim <- DHARMa::simulateResiduals(model_delta_ts)
plot(sim)

drop1(model_delta_ts, test = "Chisq")
emmeans(model_delta_ts, pairwise ~ treatment, type = "response")
contrasts_delta[["total_sum_of_interactions"]]<-as.data.frame(summary(emmeans(model_delta_ts, pairwise ~ treatment, type = "response"))$contrasts)
#add significance labels
contrasts_delta[["total_sum_of_interactions"]]$label<- ifelse(contrasts_delta[["total_sum_of_interactions"]]$p.value < 0.001, "***",
                                                            ifelse(contrasts_delta[["total_sum_of_interactions"]]$p.value < 0.01, "**",
                                                                   ifelse(contrasts_delta[["total_sum_of_interactions"]]$p.value < 0.05, "*", "ns")))


#total number of interactions
hist(before_after_infection_delta$total_number_interactions)
model_delta_tn<-glmmTMB(total_number_interactions~treatment, 
                        data = before_after_infection_delta,
                        #dispformula = ~treatment,
                        family = gaussian()
)

sim <- DHARMa::simulateResiduals(model_delta_tn)
plot(sim)

drop1(model_delta_tn, test = "Chisq")
emmeans(model_delta_tn, pairwise ~ treatment, type = "response")
contrasts_delta[["total_number_interactions"]]<-as.data.frame(summary(emmeans(model_delta_tn, pairwise ~ treatment, type = "response"))$contrasts)
#add significance labels
contrasts_delta[["total_number_interactions"]]$label<- ifelse(contrasts_delta[["total_number_interactions"]]$p.value < 0.001, "***",
                                                              ifelse(contrasts_delta[["total_number_interactions"]]$p.value < 0.01, "**",
                                                                     ifelse(contrasts_delta[["total_number_interactions"]]$p.value < 0.05, "*", "ns")))

##glmm for parameters for exp1 
##save in file for plots
contrasts_delta<-list()


###################################
#parameters where deltas are being looked at:

#strength
hist(before_after_infection_delta$strength_mean)
hist((before_after_infection_delta%>%filter(treatment == "bB"))$strength_mean)
model_delta_s<-lm(strength_mean~treatment, 
                       data = before_after_infection_delta
                       )

sim <- DHARMa::simulateResiduals(model_delta_s)
plot(sim)

model_delta_s_glm<-glmmTMB(strength_mean~treatment, 
                  data = before_after_infection_delta,
                  family = gaussian()
)

sim <- DHARMa::simulateResiduals(model_delta_s_glm)
plot(sim)
drop1(model_delta_s_glm, test = "Chisq")
summary(emmeans(model_delta_s_glm, pairwise ~ treatment, type = "response"))$contrasts

#subsets on global time
subset_means<- before_after_infection_delta%>%filter(globaltime=="8.75-10h") #1.25-2.5h, "10-11.25h" ,"8.75-10h"



model_delta_s_subset<-glmmTMB(strength_mean~treatment, 
                           data = subset_means,
                           family = gaussian()
)

sim <- DHARMa::simulateResiduals(model_delta_s_subset)
plot(sim)

drop1(model_delta_s_subset, test = "Chisq")
summary(emmeans(model_delta_s_subset, pairwise ~ treatment, type = "response"))$contrasts

drop1(model_delta_s, test = "Chisq")
contrasts_delta[["strength_mean"]]<-as.data.frame(summary(emmeans(model_delta_s, pairwise ~ treatment, type = "response"))$contrasts)
#add significance labels
contrasts_delta[["strength_mean"]]$label<- ifelse(contrasts_delta[["strength_mean"]]$p.value < 0.001, "***",
                                                  ifelse(contrasts_delta[["strength_mean"]]$p.value < 0.01, "**",
                                                         ifelse(contrasts_delta[["strength_mean"]]$p.value < 0.05, "*", "ns")))
#density_aggr
hist(before_after_infection_delta$density_aggr)
model_delta_d<-lm(density_aggr~treatment, 
                  data = before_after_infection_delta,
                  #dispformula = ~1,
                  #family = gaussian()
)

sim <- DHARMa::simulateResiduals(model_delta_d)
plot(sim)

drop1(model_delta_d, test = "Chisq")
emmeans(model_delta_d, pairwise ~ treatment, type = "response")
contrasts_delta[["density_aggr"]]<-as.data.frame(summary(emmeans(model_delta_d, pairwise ~ treatment, type = "response"))$contrasts)
#add significance labels
contrasts_delta[["density_aggr"]]$label<- ifelse(contrasts_delta[["density_aggr"]]$p.value < 0.001, "***",
                                                  ifelse(contrasts_delta[["density_aggr"]]$p.value < 0.01, "**",
                                                         ifelse(contrasts_delta[["density_aggr"]]$p.value < 0.05, "*", "ns")))


#assortativity
hist(before_after_infection_delta$assortativity)
model_delta_a<-lm(assortativity~treatment, 
                  data = before_after_infection_delta
                  #dispformula = ~1,
                  #family = gaussian()
)

sim <- DHARMa::simulateResiduals(model_delta_a)
plot(sim)

drop1(model_delta_a, test = "Chisq")
emmeans(model_delta_a, pairwise ~ treatment, type = "response")
contrasts_delta[["assortativity"]]<-as.data.frame(summary(emmeans(model_delta_d, pairwise ~ treatment, type = "response"))$contrasts)
#add significance labels
contrasts_delta[["assortativity"]]$label<- ifelse(contrasts_delta[["assortativity"]]$p.value < 0.001, "***",
                                                 ifelse(contrasts_delta[["assortativity"]]$p.value < 0.01, "**",
                                                        ifelse(contrasts_delta[["assortativity"]]$p.value < 0.05, "*", "ns")))


###strength_sd
hist(before_after_infection_delta$strength_sd)
model_delta_ssd<-lm(strength_sd~treatment, 
                  data = before_after_infection_delta,
                  #dispformula = ~1,
                  #family = gaussian()
)

sim <- DHARMa::simulateResiduals(model_delta_ssd)
plot(sim)

drop1(model_delta_ssd, test = "Chisq")
emmeans(model_delta_ssd, pairwise ~ treatment, type = "response")
contrasts_delta[["strength_sd"]]<-as.data.frame(summary(emmeans(model_delta_sdd, pairwise ~ treatment, type = "response"))$contrasts)
#add significance labels
contrasts_delta[["strength_sd"]]$label<- ifelse(contrasts_delta[["strength_sd"]]$p.value < 0.001, "***",
                                                 ifelse(contrasts_delta[["strength_sd"]]$p.value < 0.01, "**",
                                                        ifelse(contrasts_delta[["strength_sd"]]$p.value < 0.05, "*", "ns")))

#burstiness
hist(before_after_infection_delta$burstiness)
model_delta_b<-lm(burstiness~treatment, 
                    data = before_after_infection_delta,
                    #dispformula = ~1,
                    #family = gaussian()
)

sim <- DHARMa::simulateResiduals(model_delta_b)
plot(sim)

drop1(model_delta_b, test = "Chisq")

#global effectiveness
hist(before_after_infection_delta$global_eff_corr)
model_delta_ge<-lm(global_eff_corr~treatment, 
                    data = before_after_infection_delta,
                    #dispformula = ~1,
                    #family = gaussian()
)

sim <- DHARMa::simulateResiduals(model_delta_ge)
plot(sim)

drop1(model_delta_ge, test = "Chisq")
emmeans(model_delta_ge, pairwise ~ treatment, type = "response")
contrasts_delta[["global_eff_corr"]]<-as.data.frame(summary(emmeans(model_delta_ge, pairwise ~ treatment, type = "response"))$contrasts)
#add significance labels
contrasts_delta[["global_eff_corr"]]$label<- ifelse(contrasts_delta[["global_eff_corr"]]$p.value < 0.001, "***",
                                                ifelse(contrasts_delta[["global_eff_corr"]]$p.value < 0.01, "**",
                                                       ifelse(contrasts_delta[["global_eff_corr"]]$p.value < 0.05, "*", "ns")))

##waitingtime_mean
hist(before_after_infection_delta$waiting_time_mean)
model_delta_wt<-lm(waiting_time_mean~treatment, 
                   data = before_after_infection_delta,
                   #dispformula = ~1,
                   #family = gaussian()
)

sim <- DHARMa::simulateResiduals(model_delta_wt)
plot(sim)

drop1(model_delta_wt, test = "Chisq")
#interactionlength mean
hist(before_after_infection_delta$interaction_length_mean)
model_delta_im<-lm(interaction_length_mean~treatment, 
                   data = before_after_infection_delta,
                   #dispformula = ~1,
                   #family = gaussian()
)

sim <- DHARMa::simulateResiduals(model_delta_im)
plot(sim)

drop1(model_delta_im, test = "Chisq")
emmeans(model_delta_im, pairwise ~ treatment, type = "response")
contrasts_delta[["interaction_length_mean"]]<-as.data.frame(summary(emmeans(model_delta_im, pairwise ~ treatment, type = "response"))$contrasts)
#add significance labels
contrasts_delta[["interaction_length_mean"]]$label<- ifelse(contrasts_delta[["interaction_length_mean"]]$p.value < 0.001, "***",
                                                    ifelse(contrasts_delta[["interaction_length_mean"]]$p.value < 0.01, "**",
                                                           ifelse(contrasts_delta[["interaction_length_mean"]]$p.value < 0.05, "*", "ns")))

#global clustering
hist(before_after_infection_delta$clustering_global)
model_delta_gc<-glm(log(clustering_global)~treatment, 
                   data = before_after_infection_delta,
                  #dispformula = ~1,
                   family = gaussian()
)

sim <- DHARMa::simulateResiduals(model_delta_gc)
plot(sim)

drop1(model_delta_gc, test = "Chisq")
emmeans(model_delta_gc, pairwise ~ treatment, type = "response")
contrasts_delta[["interaction_length_mean"]]<-as.data.frame(summary(emmeans(model_delta_im, pairwise ~ treatment, type = "response"))$contrasts)
#add significance labels
contrasts_delta[["interaction_length_mean"]]$label<- ifelse(contrasts_delta[["interaction_length_mean"]]$p.value < 0.001, "***",
                                                            ifelse(contrasts_delta[["interaction_length_mean"]]$p.value < 0.01, "**",
                                                                   ifelse(contrasts_delta[["interaction_length_mean"]]$p.value < 0.05, "*", "ns")))

