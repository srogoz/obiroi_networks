mean_overtimes<-before_after_infection_merged%>%group_by(colony_name, globaltime,infection) %>%
  summarise(
    across(where(is.numeric), mean, na.rm = TRUE),
    .groups = "drop"
  )
#mean over before and after infection

prep_delta_1<-before_after_infection_mean%>%filter(globaltime=="0-1.25h")
prep_delta_2<-before_after_infection_mean%>%filter(globaltime=="1.25-2.5h")
prep_delta_3<-before_after_infection_mean%>%filter(globaltime=="10h-11.25h")
prep_delta_4<-before_after_infection_mean%>%filter(globaltime=="8.75h-10h")

prep_delta_inf<-before_after_infection_mean%>%filter(globaltime=="fungus exp")

param<-"strength_mean"
param_name<-paste0("delta_", param)
before_after_infection_delta<-data_frame()
before_after_infection_delta[[param_name]]<-prep_delta_1[[param]]-prep_delta_inf[[param]]

prep_delta_2<-prep_delta_2%>%mutate(delta = mean_value- prep_delta_inf$mean_value)
prep_delta_2<-prep_delta_2%>%mutate(delta_norm = (mean_value- prep_delta_inf$mean_value)/mean_value)

prep_delta_3<-prep_delta_3%>%mutate(delta = mean_value- prep_delta_inf$mean_value)
prep_delta_3<-prep_delta_3%>%mutate(delta_norm = (mean_value- prep_delta_inf$mean_value)/mean_value)

prep_delta_4<-prep_delta_4%>%mutate(delta = mean_value- prep_delta_inf$mean_value)
prep_delta_4<-prep_delta_4%>%mutate(delta_norm = (mean_value- prep_delta_inf$mean_value)/mean_value)

#
prep_delta_all<-bind_rows(prep_delta_1, prep_delta_2, prep_delta_3, prep_delta_4)

prep_delta_all<-prep_delta_all%>%mutate(treatment = sub("\\d+.*$", "", colony_name))

delta_strength<-prep_delta_all%>%filter(parameter=="strength_mean"
)


hist(delta_strength$delta_norm)

##glmm for strength
###################################
model_delta_s<-glmmTMB(delta_norm~treatment + (1|colony_name), 
                       data = delta_strength,
                       dispformula = ~1,
                       family = gaussian())

sim <- DHARMa::simulateResiduals(model_delta_s)
plot(sim)

pairs <- emmeans(model_delta_s, pairwise ~ treatment, type = "response")
stat_test <- pairs$contrasts
stat_test

########
model_delta_s<-glmmTMB(delta~treatment + (1|colony_name), 
                       data = delta_strength,
                       dispformula = ~1,
                       family = gaussian())

sim <- DHARMa::simulateResiduals(model_delta_s)
plot(sim)

pairs <- emmeans(model_delta_s, pairwise ~ treatment, type = "response")
stat_test <- pairs$contrasts
stat_test

#for density
hist(delta_density$delta)
delta_density<-delta%>%filter(parameter=="density_aggr")
model_delta_d<-glmmTMB(delta_norm~treatment, 
                       data = delta_density,
                       family = gaussian())

sim <- DHARMa::simulateResiduals(model_delta_d)
plot(sim)

pairs <- emmeans(model_delta_d, pairwise ~ treatment, type = "response")
stat_test <- pairs$contrasts
stat_test