# plot log-log to get the powerlaw distribution for a>2 for mean and sd to exist and be finite
# for the whole colony
#filter out zeros for waiting times
library(dplyr)
waiting_times<-data.frame(wt=diff(sort(network_obj$onset)) )
wait <- waiting_times %>%
  filter(wt > 0)
waiting_times <- waiting_times[waiting_times$wt > 0, , drop = FALSE]
waiting_times_5<-data.frame(wt=diff(sort(network_obj_5$onset)))
waiting_times_5<-waiting_times_5[waiting_times_5$wt>0,,drop = FALSE]
#############################
####
B_urst<-(sd(waiting_times_5)-mean(waiting_times_5))/(sd(waiting_times_5)+mean(waiting_times_5))
############
#plot waiting times with right distribution

max_wt<- max(wait$wt)
max_wt_5<- max(waiting_times_5$wt)
library(ggplot2)
hist <- ggplot(wait, aes(x = wt))+
  geom_histogram(binwidth = 0.001, fill ="purple", color ="black" , alpha = 0.7) +
  scale_x_log10()+
  scale_y_log10()+
  labs( x =" log(waiting time [frame])", 
        y = "log (frequency)",
        title = "log-log distribution of waiting times"
  ) + 
  #coord_cartesian(ylim = c(0, 3000), xlim= c(0,max_wt)) + 
  #scale_y_log10(limits = c(1, 100))+
  theme_minimal()
hist

##################
library(poweRlaw)
wait_list<-wait$wt

#create discrete powerlaw object
pl_obj <- displ$new(wait_list)
#get x_min and gamme estimates
est_pl <- estimate_xmin(pl_obj)

est_pl$xmin
est_pl$pars
#####################
#check fit of distributions
library(fitdistrplus)
x<-waiting_times_5$wt
fit_ln<-fitdist(x, "lnorm")   # log-normal
fit_wb  <- fitdist(x, "weibull")

gof <- gofstat(list(fit_ln, fit_wb))
print(gof)
###############bootstrapping estimation of powerlaw to check goodness of powerlaw fit

bs_pl <- bootstrap_p(pl_obj, no_of_sims=1000, threads=8, seed = 123)

df_bs_pl <- bs_pl$bootstraps
ggplot(data=df_bs_pl, aes(pars)) + geom_histogram() + labs(x="gamma", y="frequency") + theme_bw()

ggplot(data=df_bs_pl, aes(xmin)) + geom_histogram() + labs(x="K_min", y="frequency") + theme_bw()
#################

data.s <- unique(wait)

d_est <- data.frame(K_min=sort(data.s)[1:(length(data.s)-2)], gamma=rep(0,length(data.s)-2), D=rep(0,length(data.s)-2))
gamma_D.min <- d_est[which.min(d_est$D), 2]

ggplot(data=df_bs_pl, aes(x=xmin, y=pars)) + labs(x="K_min", y="gamma") + theme_bw() + 
  geom_point(shape=21, colour="black", fill="red", size=0.5, stroke=2, 
             position = position_jitter(), alpha=0.6) +
  geom_vline(xintercept=K.min_D.min, colour="blue") +
  geom_hline(yintercept=gamma_D.min, colour="blue") +
  annotate("text", x=K.min_D.min, y=min(df_bs_pl$pars), label=K.min_D.min, col="blue") +
  annotate("text", x=min(df_bs_pl$xmin), y=gamma_D.min, label=round(gamma_D.min, digits=2), col="blue")
