############################################################################   
############################################################################
##                                                                        ##
##  R-code: Meta-regression to estimate the effect of traits on           ##
##          the demographic sensitivity of populations to precipiration   ##
##                                                                        ##
############################################################################
############################################################################

rm (list = ls ( ) ) 

# Load packages
library(brms)
library(data.table)
library(DHARMa)
library(mgcv)
library(ggplot2)
library(tidyverse)
library(tictoc)
library(purrr)
library(ggdist)
library(ape)
library(phangorn)
library(patchwork)


setwd("")

# Import dataset
data_meta_regression <- read.csv("final_dataset.txt", sep="")

data_meta_regression <- data_meta_regression[!is.na(data_meta_regression$coop_num),]

########################
#### Load phylogeny ####
########################

MyTree <- read.tree("mcctree_mctavish_v1.6.nex") 

tree_phylo <- MyTree$`TREE1=`

setdiff(data_meta_regression$species, tree_phylo$tip.label) # check that all species are in the phylogeny

tree_phylo <- keep.tip(tree_phylo, data_meta_regression$species) # drop tip from phylo: remove the species not present in the dataset

phylo_mcct <- vcv.phylo(tree_phylo) # transform into a variance covariance matrix

data_meta_regression$phylo <- data_meta_regression$species

#########################
#### meta-regression ####
#########################

meta_regression_cov <- brm(
  abs_slope_rainfall ~ 1 + 
    scale(logmass) +
    scale(hand_wing_index) +
    scale(sample_size) + 
    scale(abs_latitude) +
    mo(coop_num) +
    mo(migration) +
    (1 | species) +
    (1 | gr(phylo, cov = phylo_mcct))
  ,
  data = data_meta_regression,
  data2 = list(phylo_mcct = phylo_mcct),
  family = Gamma(link = "log"),
  prior = c(
    set_prior("normal(0,1)", class = "b"),
    set_prior("normal(0,1)", class = "Intercept"),
    set_prior("exponential(2)", class = "sd")),
  cores = 3, 
  chains = 3, 
  warmup = 1000, 
  iter = 2000, 
  thin = 1, 
  seed = 123,
  control = list(adapt_delta = 0.9),
  backend = "cmdstanr"
)

summary(meta_regression_cov, prob=0.89)

pp_check(meta_regression_cov, ndraws = 100) + xlim(c(0,1)) + theme_bw()

posterior <- as.matrix(meta_regression_cov)
dat_plot <- as.data.frame(posterior)


########################
#### Plot body mass ####
########################

x2.sim = seq(min(scale(data_meta_regression$logmass)[,1]),
             max(scale(data_meta_regression$logmass)[,1]),
             by =  0.1) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- exp(dat_plot$b_Intercept + dat_plot$b_scalelogmass * (x2.sim[i]))
}

# calculate quantiles of predictions
bayes.c.eff.mean <- exp(median(dat_plot$b_Intercept) + median(dat_plot$b_scalelogmass) * (x2.sim)) 
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.05)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.95)))
bayes.c.eff.lower.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.15)))
bayes.c.eff.upper.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.85)))
bayes.c.eff.lower.bis2 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis2 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.75)))
bayes.c.eff.lower.bis3 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.35)))
bayes.c.eff.upper.bis3 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.65)))
bayes.c.eff.lower.bis4 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.45)))
bayes.c.eff.upper.bis4 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.55)))
plot.dat <- data.frame(x2.sim, bayes.c.eff.mean,
                       bayes.c.eff.lower, bayes.c.eff.upper,
                       bayes.c.eff.lower.bis, bayes.c.eff.upper.bis,
                       bayes.c.eff.lower.bis2, bayes.c.eff.upper.bis2,
                       bayes.c.eff.lower.bis3, bayes.c.eff.upper.bis3,
                       bayes.c.eff.lower.bis4, bayes.c.eff.upper.bis4)

inv_fun_mass <- function(x){(x*sd(data_meta_regression$logmass))+mean(data_meta_regression$logmass)}


plot_mass <- ggplot(plot.dat, aes(x = inv_fun_mass(x2.sim), y = bayes.c.eff.mean)) +
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "cornflowerblue", alpha = 0.15)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "cornflowerblue", alpha = 0.15)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "cornflowerblue", alpha = 0.15)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "cornflowerblue", alpha = 0.15)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "cornflowerblue", alpha = 0.15)+
  geom_line(color = "cornflowerblue", linewidth = 1.8, alpha=0.7)+
  xlab("Log mass")+
  ylab(expression(paste("Absolute rainfall sensitivity (log ", lambda, " change per 100mm)")))+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.title.y=element_text(size=12))
plot_mass


posterior_mass <- as.data.frame(dat_plot$b_scalelogmass)

plot_posterior_mass <- ggplot(posterior_mass, aes(x=`dat_plot$b_scalelogmass`))+
      stat_slab(aes(fill = after_stat(level), alpha = after_stat(level)), .width = c(.1,.3,.5,.7,.9, 1)) +
      scale_fill_manual(values = c("cornflowerblue", "cornflowerblue", "cornflowerblue", "cornflowerblue", "cornflowerblue", "cornflowerblue"))+
      geom_vline(aes(xintercept=0), color="grey10", linewidth=1)+
      xlab("Posterior estimate")+
      theme_bw()+
      theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())+
      theme(legend.position = "none")+
      theme(axis.line.x = element_line(colour = "black"),
            axis.text.x = element_text(size=12),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            panel.background = element_blank(), panel.border = element_blank())
plot_posterior_mass 

plot_mass_final <- plot_mass + inset_element(plot_posterior_mass, 0.5, 0.75, 0.95, 0.95)
plot_mass_final


########################
#### Plot sociality ####
########################

# no family
int.sim  <-  exp(dat_plot$b_Intercept) 
# family
int.sim2 <-  exp(dat_plot$b_Intercept + 
  dat_plot$bsp_mocoop_num * (dat_plot$`simo_mocoop_num1[1]` / 0.5))
# cooperative family
int.sim3 <-  exp(dat_plot$b_Intercept + 
  dat_plot$bsp_mocoop_num * (dat_plot$`simo_mocoop_num1[1]` / 0.5) +
  dat_plot$bsp_mocoop_num * (dat_plot$`simo_mocoop_num1[2]` / 0.5))


# calculate quantiles of predictions
# no family
bayes.c.eff.mean <- exp(median(dat_plot$b_Intercept))
bayes.c.eff.lower <- quantile(int.sim, probs = c(0.05))
bayes.c.eff.upper <- quantile(int.sim, probs = c(0.95))
bayes.c.eff.lower.bis <- quantile(int.sim, probs = c(0.15))
bayes.c.eff.upper.bis <- quantile(int.sim, probs = c(0.85))
bayes.c.eff.lower.bis2 <- quantile(int.sim, probs = c(0.25))
bayes.c.eff.upper.bis2 <- quantile(int.sim, probs = c(0.75))
bayes.c.eff.lower.bis3 <- quantile(int.sim, probs = c(0.35))
bayes.c.eff.upper.bis3 <- quantile(int.sim, probs = c(0.65))
bayes.c.eff.lower.bis4 <- quantile(int.sim, probs = c(0.45))
bayes.c.eff.upper.bis4 <- quantile(int.sim, probs = c(0.55))
plot.dat <- data.frame(bayes.c.eff.mean,
                       bayes.c.eff.lower, bayes.c.eff.upper,
                       bayes.c.eff.lower.bis, bayes.c.eff.upper.bis,
                       bayes.c.eff.lower.bis2, bayes.c.eff.upper.bis2,
                       bayes.c.eff.lower.bis3, bayes.c.eff.upper.bis3,
                       bayes.c.eff.lower.bis4, bayes.c.eff.upper.bis4)

# family
bayes.c.eff.mean <- exp(median(dat_plot$b_Intercept) + 
  median(dat_plot$bsp_mocoop_num) * (median(dat_plot$`simo_mocoop_num1[1]`) / 0.5))
bayes.c.eff.lower <- quantile(int.sim2, probs = c(0.05))
bayes.c.eff.upper <- quantile(int.sim2, probs = c(0.95))
bayes.c.eff.lower.bis <- quantile(int.sim2, probs = c(0.15))
bayes.c.eff.upper.bis <- quantile(int.sim2, probs = c(0.85))
bayes.c.eff.lower.bis2 <- quantile(int.sim2, probs = c(0.25))
bayes.c.eff.upper.bis2 <- quantile(int.sim2, probs = c(0.75))
bayes.c.eff.lower.bis3 <- quantile(int.sim2, probs = c(0.35))
bayes.c.eff.upper.bis3 <- quantile(int.sim2, probs = c(0.65))
bayes.c.eff.lower.bis4 <- quantile(int.sim2, probs = c(0.45))
bayes.c.eff.upper.bis4 <- quantile(int.sim2, probs = c(0.55))
plot.dat.new <- data.frame(bayes.c.eff.mean,
                           bayes.c.eff.lower, bayes.c.eff.upper,
                           bayes.c.eff.lower.bis, bayes.c.eff.upper.bis,
                           bayes.c.eff.lower.bis2, bayes.c.eff.upper.bis2,
                           bayes.c.eff.lower.bis3, bayes.c.eff.upper.bis3,
                           bayes.c.eff.lower.bis4, bayes.c.eff.upper.bis4)

# cooperative family
bayes.c.eff.mean <- exp(median(dat_plot$b_Intercept) + 
  median(dat_plot$bsp_mocoop_num) * (median(dat_plot$`simo_mocoop_num1[1]`) / 0.5) +
  median(dat_plot$bsp_mocoop_num) * (median(dat_plot$`simo_mocoop_num1[2]`) / 0.5))
bayes.c.eff.lower <- quantile(int.sim3, probs = c(0.05))
bayes.c.eff.upper <- quantile(int.sim3, probs = c(0.95))
bayes.c.eff.lower.bis <- quantile(int.sim3, probs = c(0.15))
bayes.c.eff.upper.bis <- quantile(int.sim3, probs = c(0.85))
bayes.c.eff.lower.bis2 <- quantile(int.sim3, probs = c(0.25))
bayes.c.eff.upper.bis2 <- quantile(int.sim3, probs = c(0.75))
bayes.c.eff.lower.bis3 <- quantile(int.sim3, probs = c(0.35))
bayes.c.eff.upper.bis3 <- quantile(int.sim3, probs = c(0.65))
bayes.c.eff.lower.bis4 <- quantile(int.sim3, probs = c(0.45))
bayes.c.eff.upper.bis4 <- quantile(int.sim3, probs = c(0.55))
plot.dat.new.2 <- data.frame(bayes.c.eff.mean,
                             bayes.c.eff.lower, bayes.c.eff.upper,
                             bayes.c.eff.lower.bis, bayes.c.eff.upper.bis,
                             bayes.c.eff.lower.bis2, bayes.c.eff.upper.bis2,
                             bayes.c.eff.lower.bis3, bayes.c.eff.upper.bis3,
                             bayes.c.eff.lower.bis4, bayes.c.eff.upper.bis4)


plot.dat <- rbind(plot.dat, plot.dat.new, plot.dat.new.2)
sp <- as.data.frame(matrix(c("NFL", "FL", "COOP"), nrow = 3, ncol = 1))
plot.dat <- cbind(plot.dat, sp)



plot_coop <- ggplot(data_meta_regression, aes(x = coop_num, y = abs_sens)) +
  xlab("")+
  ylab(expression(paste("Absolute rainfall sensitivity (log ", lambda, " change per 100mm)")))+
  ylim(c(0.04, 0.21))+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.title.y=element_text(size=12))

plot_coop <- plot_coop + annotate("rect", xmin = 0.7, xmax = 1.3, ymin =  plot.dat$bayes.c.eff.lower.bis4[1], ymax = plot.dat$bayes.c.eff.upper.bis4[1], alpha=0.15, fill="cornflowerblue")
plot_coop <- plot_coop + annotate("rect", xmin = 1.7, xmax = 2.3, ymin =  plot.dat$bayes.c.eff.lower.bis4[2], ymax = plot.dat$bayes.c.eff.upper.bis4[2], alpha=0.15, fill="cornflowerblue")
plot_coop <- plot_coop + annotate("rect", xmin = 2.7, xmax = 3.3, ymin =  plot.dat$bayes.c.eff.lower.bis4[3], ymax = plot.dat$bayes.c.eff.upper.bis4[3], alpha=0.15, fill="cornflowerblue")

plot_coop <- plot_coop + annotate("rect", xmin = 0.7, xmax = 1.3, ymin =  plot.dat$bayes.c.eff.lower.bis3[1], ymax = plot.dat$bayes.c.eff.upper.bis3[1], alpha=0.15, fill="cornflowerblue")
plot_coop <- plot_coop + annotate("rect", xmin = 1.7, xmax = 2.3, ymin =  plot.dat$bayes.c.eff.lower.bis3[2], ymax = plot.dat$bayes.c.eff.upper.bis3[2], alpha=0.15, fill="cornflowerblue")
plot_coop <- plot_coop + annotate("rect", xmin = 2.7, xmax = 3.3, ymin =  plot.dat$bayes.c.eff.lower.bis3[3], ymax = plot.dat$bayes.c.eff.upper.bis3[3], alpha=0.15, fill="cornflowerblue")

plot_coop <- plot_coop + annotate("rect", xmin = 0.7, xmax = 1.3, ymin =  plot.dat$bayes.c.eff.lower.bis2[1], ymax = plot.dat$bayes.c.eff.upper.bis2[1], alpha=0.15, fill="cornflowerblue")
plot_coop <- plot_coop + annotate("rect", xmin = 1.7, xmax = 2.3, ymin =  plot.dat$bayes.c.eff.lower.bis2[2], ymax = plot.dat$bayes.c.eff.upper.bis2[2], alpha=0.15, fill="cornflowerblue")
plot_coop <- plot_coop + annotate("rect", xmin = 2.7, xmax = 3.3, ymin =  plot.dat$bayes.c.eff.lower.bis2[3], ymax = plot.dat$bayes.c.eff.upper.bis2[3], alpha=0.15, fill="cornflowerblue")

plot_coop <- plot_coop + annotate("rect", xmin = 0.7, xmax = 1.3, ymin =  plot.dat$bayes.c.eff.lower.bis[1], ymax = plot.dat$bayes.c.eff.upper.bis[1], alpha=0.15, fill="cornflowerblue")
plot_coop <- plot_coop + annotate("rect", xmin = 1.7, xmax = 2.3, ymin =  plot.dat$bayes.c.eff.lower.bis[2], ymax = plot.dat$bayes.c.eff.upper.bis[2], alpha=0.15, fill="cornflowerblue")
plot_coop <- plot_coop + annotate("rect", xmin = 2.7, xmax = 3.3, ymin =  plot.dat$bayes.c.eff.lower.bis[3], ymax = plot.dat$bayes.c.eff.upper.bis[3], alpha=0.15, fill="cornflowerblue")

plot_coop <- plot_coop + annotate("rect", xmin = 0.7, xmax = 1.3, ymin =  plot.dat$bayes.c.eff.lower[1], ymax = plot.dat$bayes.c.eff.upper[1], alpha=0.15, fill="cornflowerblue")
plot_coop <- plot_coop + annotate("rect", xmin = 1.7, xmax = 2.3, ymin =  plot.dat$bayes.c.eff.lower[2], ymax = plot.dat$bayes.c.eff.upper[2], alpha=0.15, fill="cornflowerblue")
plot_coop <- plot_coop + annotate("rect", xmin = 2.7, xmax = 3.3, ymin =  plot.dat$bayes.c.eff.lower[3], ymax = plot.dat$bayes.c.eff.upper[3], alpha=0.15, fill="cornflowerblue")

plot_coop <- plot_coop + annotate("segment", x = 0.7, xend = 1.3, y = plot.dat$bayes.c.eff.mean[1], yend = plot.dat$bayes.c.eff.mean[1], size=1.8, alpha=0.7, colour="cornflowerblue")
plot_coop <- plot_coop + annotate("segment", x = 1.7, xend = 2.3, y = plot.dat$bayes.c.eff.mean[2], yend = plot.dat$bayes.c.eff.mean[2], size=1.8, alpha=0.7, colour="cornflowerblue")
plot_coop <- plot_coop + annotate("segment", x = 2.7, xend = 3.3, y = plot.dat$bayes.c.eff.mean[3], yend = plot.dat$bayes.c.eff.mean[3], size=1.8, alpha=0.7, colour="cornflowerblue")

plot_coop <- plot_coop + scale_x_continuous(breaks = c(1,2,3), labels = c("Non-family living", "Family living", "Cooperative breeding"))
plot_coop




posterior_coop <- as.data.frame(dat_plot$bsp_mocoop_num)

plot_posterior_coop <- ggplot(posterior_coop, aes(x=`dat_plot$bsp_mocoop_num`))+
  stat_slab(aes(fill = after_stat(level), alpha = after_stat(level)), .width = c(.1,.3,.5,.7,.9, 1)) +
  scale_fill_manual(values = c("cornflowerblue", "cornflowerblue", "cornflowerblue", "cornflowerblue", "cornflowerblue", "cornflowerblue"))+
  geom_vline(aes(xintercept=0), color="grey40", linewidth=1)+
  xlab("Posterior estimate")+
  theme_bw()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  theme(legend.position = "none")+
  theme(axis.line.x = element_line(colour = "black"),
        axis.text.x = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank())
plot_posterior_coop


plot_coop_final <- plot_coop + inset_element(plot_posterior_coop, 0.5, 0.75, 0.95, 0.95)
plot_coop_final



#######################
#### Plot latitude ####
#######################

x2.sim = seq(min(scale(data_meta_regression$abs_latitude)[,1]),
             max(scale(data_meta_regression$abs_latitude)[,1]),
             by =  0.1) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- exp(dat_plot$b_Intercept + dat_plot$b_scaleabs_latitude * (x2.sim[i]))
}

# calculate quantiles of predictions
bayes.c.eff.mean <- exp(median(dat_plot$b_Intercept) + median(dat_plot$b_scaleabs_latitude) * (x2.sim)) 
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.05)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.95)))
bayes.c.eff.lower.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.15)))
bayes.c.eff.upper.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.85)))
bayes.c.eff.lower.bis2 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis2 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.75)))
bayes.c.eff.lower.bis3 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.35)))
bayes.c.eff.upper.bis3 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.65)))
bayes.c.eff.lower.bis4 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.45)))
bayes.c.eff.upper.bis4 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.55)))
plot.dat <- data.frame(x2.sim, bayes.c.eff.mean,
                       bayes.c.eff.lower, bayes.c.eff.upper,
                       bayes.c.eff.lower.bis, bayes.c.eff.upper.bis,
                       bayes.c.eff.lower.bis2, bayes.c.eff.upper.bis2,
                       bayes.c.eff.lower.bis3, bayes.c.eff.upper.bis3,
                       bayes.c.eff.lower.bis4, bayes.c.eff.upper.bis4)

inv_fun_latitude <- function(x){(x*sd(data_meta_regression$abs_latitude))+mean(data_meta_regression$abs_latitude)}


plot_latitude <- ggplot(plot.dat, aes(x = inv_fun_latitude(x2.sim), y = bayes.c.eff.mean)) +
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "cornflowerblue", alpha = 0.15)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "cornflowerblue", alpha = 0.15)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "cornflowerblue", alpha = 0.15)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "cornflowerblue", alpha = 0.15)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "cornflowerblue", alpha = 0.15)+
  geom_line(color = "cornflowerblue", linewidth = 1.8, alpha=0.7)+
  xlab("Absolute latitude")+
  ylab(expression(paste("Absolute rainfall sensitivity (log ", lambda, " change per 100mm)")))+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.title.y=element_text(size=12))
plot_latitude


posterior_latitude <- as.data.frame(dat_plot$b_scaleabs_latitude)

plot_posterior_latitude <- ggplot(posterior_latitude, aes(x=`dat_plot$b_scaleabs_latitude`))+
  stat_slab(aes(fill = after_stat(level), alpha = after_stat(level)), .width = c(.1,.3,.5,.7,.9, 1)) +
  scale_fill_manual(values = c("cornflowerblue", "cornflowerblue", "cornflowerblue", "cornflowerblue", "cornflowerblue", "cornflowerblue"))+
  geom_vline(aes(xintercept=0), color="grey10", linewidth=1)+
  xlab("Posterior estimate")+
  theme_bw()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  theme(legend.position = "none")+
  theme(axis.line.x = element_line(colour = "black"),
        axis.text.x = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank())
plot_posterior_latitude 

plot_latitude_final <- plot_latitude + inset_element(plot_posterior_latitude, 0.05, 0.75, 0.5, 0.95)
plot_latitude_final




##################
#### Plot HWI ####
##################

x2.sim = seq(min(scale(data_meta_regression$hand_wing_index)[,1]),
             max(scale(data_meta_regression$hand_wing_index)[,1]),
             by =  0.1) 

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- exp(dat_plot$b_Intercept + dat_plot$b_scalehand_wing_index * (x2.sim[i]))
}

# calculate quantiles of predictions
bayes.c.eff.mean <- exp(median(dat_plot$b_Intercept) + median(dat_plot$b_scalehand_wing_index) * (x2.sim)) 
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.05)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.95)))
bayes.c.eff.lower.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.15)))
bayes.c.eff.upper.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.85)))
bayes.c.eff.lower.bis2 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis2 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.75)))
bayes.c.eff.lower.bis3 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.35)))
bayes.c.eff.upper.bis3 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.65)))
bayes.c.eff.lower.bis4 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.45)))
bayes.c.eff.upper.bis4 <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.55)))
plot.dat <- data.frame(x2.sim, bayes.c.eff.mean,
                       bayes.c.eff.lower, bayes.c.eff.upper,
                       bayes.c.eff.lower.bis, bayes.c.eff.upper.bis,
                       bayes.c.eff.lower.bis2, bayes.c.eff.upper.bis2,
                       bayes.c.eff.lower.bis3, bayes.c.eff.upper.bis3,
                       bayes.c.eff.lower.bis4, bayes.c.eff.upper.bis4)

inv_fun_hwi <- function(x){(x*sd(data_meta_regression$hand_wing_index))+mean(data_meta_regression$hand_wing_index)}


plot_hwi <- ggplot(plot.dat, aes(x = inv_fun_hwi(x2.sim), y = bayes.c.eff.mean)) +
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis4, ymax = bayes.c.eff.upper.bis4), fill = "cornflowerblue", alpha = 0.15)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis3, ymax = bayes.c.eff.upper.bis3), fill = "cornflowerblue", alpha = 0.15)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis2, ymax = bayes.c.eff.upper.bis2), fill = "cornflowerblue", alpha = 0.15)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower.bis, ymax = bayes.c.eff.upper.bis), fill = "cornflowerblue", alpha = 0.15)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "cornflowerblue", alpha = 0.15)+
  geom_line(color = "cornflowerblue", linewidth = 1.8, alpha=0.7)+
  xlab("Hand-wing index")+
  ylab(expression(paste("Absolute rainfall sensitivity (log ", lambda, " change per 100mm)")))+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.title.y=element_text(size=12))
plot_hwi


posterior_hwi <- as.data.frame(dat_plot$b_scalehand_wing_index)

plot_posterior_hwi <- ggplot(posterior_hwi, aes(x=`dat_plot$b_scalehand_wing_index`))+
  stat_slab(aes(fill = after_stat(level), alpha = after_stat(level)), .width = c(.1,.3,.5,.7,.9, 1)) +
  scale_fill_manual(values = c("cornflowerblue", "cornflowerblue", "cornflowerblue", "cornflowerblue", "cornflowerblue", "cornflowerblue"))+
  geom_vline(aes(xintercept=0), color="grey10", linewidth=1)+
  xlab("Posterior estimate")+
  theme_bw()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  theme(legend.position = "none")+
  theme(axis.line.x = element_line(colour = "black"),
        axis.text.x = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank())
plot_posterior_hwi

plot_hwi_final <- plot_hwi + inset_element(plot_posterior_hwi, 0.05, 0.75, 0.5, 0.95)
plot_hwi_final




########################
#### Plot migration ####
########################

# resident
int.sim  <-  exp(dat_plot$b_Intercept) 
# partial migrant
int.sim2 <-  exp(dat_plot$b_Intercept + 
                   dat_plot$bsp_momigration * (dat_plot$`simo_momigration1[1]` / 0.5))
# full migrant
int.sim3 <-  exp(dat_plot$b_Intercept + 
                   dat_plot$bsp_momigration * (dat_plot$`simo_momigration1[1]` / 0.5) +
                   dat_plot$bsp_momigration * (dat_plot$`simo_momigration1[2]` / 0.5))


# calculate quantiles of predictions
# resident
bayes.c.eff.mean <- exp(median(dat_plot$b_Intercept))
bayes.c.eff.lower <- quantile(int.sim, probs = c(0.05))
bayes.c.eff.upper <- quantile(int.sim, probs = c(0.95))
bayes.c.eff.lower.bis <- quantile(int.sim, probs = c(0.15))
bayes.c.eff.upper.bis <- quantile(int.sim, probs = c(0.85))
bayes.c.eff.lower.bis2 <- quantile(int.sim, probs = c(0.25))
bayes.c.eff.upper.bis2 <- quantile(int.sim, probs = c(0.75))
bayes.c.eff.lower.bis3 <- quantile(int.sim, probs = c(0.35))
bayes.c.eff.upper.bis3 <- quantile(int.sim, probs = c(0.65))
bayes.c.eff.lower.bis4 <- quantile(int.sim, probs = c(0.45))
bayes.c.eff.upper.bis4 <- quantile(int.sim, probs = c(0.55))
plot.dat <- data.frame(bayes.c.eff.mean,
                       bayes.c.eff.lower, bayes.c.eff.upper,
                       bayes.c.eff.lower.bis, bayes.c.eff.upper.bis,
                       bayes.c.eff.lower.bis2, bayes.c.eff.upper.bis2,
                       bayes.c.eff.lower.bis3, bayes.c.eff.upper.bis3,
                       bayes.c.eff.lower.bis4, bayes.c.eff.upper.bis4)

# partial migrant
bayes.c.eff.mean <- exp(median(dat_plot$b_Intercept) + 
                          median(dat_plot$bsp_momigration) * (median(dat_plot$`simo_momigration1[1]`) / 0.5))
bayes.c.eff.lower <- quantile(int.sim2, probs = c(0.05))
bayes.c.eff.upper <- quantile(int.sim2, probs = c(0.95))
bayes.c.eff.lower.bis <- quantile(int.sim2, probs = c(0.15))
bayes.c.eff.upper.bis <- quantile(int.sim2, probs = c(0.85))
bayes.c.eff.lower.bis2 <- quantile(int.sim2, probs = c(0.25))
bayes.c.eff.upper.bis2 <- quantile(int.sim2, probs = c(0.75))
bayes.c.eff.lower.bis3 <- quantile(int.sim2, probs = c(0.35))
bayes.c.eff.upper.bis3 <- quantile(int.sim2, probs = c(0.65))
bayes.c.eff.lower.bis4 <- quantile(int.sim2, probs = c(0.45))
bayes.c.eff.upper.bis4 <- quantile(int.sim2, probs = c(0.55))
plot.dat.new <- data.frame(bayes.c.eff.mean,
                           bayes.c.eff.lower, bayes.c.eff.upper,
                           bayes.c.eff.lower.bis, bayes.c.eff.upper.bis,
                           bayes.c.eff.lower.bis2, bayes.c.eff.upper.bis2,
                           bayes.c.eff.lower.bis3, bayes.c.eff.upper.bis3,
                           bayes.c.eff.lower.bis4, bayes.c.eff.upper.bis4)

# full migrant
bayes.c.eff.mean <- exp(median(dat_plot$b_Intercept) + 
                          median(dat_plot$bsp_momigration) * (median(dat_plot$`simo_momigration1[1]`) / 0.5) +
                          median(dat_plot$bsp_momigration) * (median(dat_plot$`simo_momigration1[2]`) / 0.5))
bayes.c.eff.lower <- quantile(int.sim3, probs = c(0.05))
bayes.c.eff.upper <- quantile(int.sim3, probs = c(0.95))
bayes.c.eff.lower.bis <- quantile(int.sim3, probs = c(0.15))
bayes.c.eff.upper.bis <- quantile(int.sim3, probs = c(0.85))
bayes.c.eff.lower.bis2 <- quantile(int.sim3, probs = c(0.25))
bayes.c.eff.upper.bis2 <- quantile(int.sim3, probs = c(0.75))
bayes.c.eff.lower.bis3 <- quantile(int.sim3, probs = c(0.35))
bayes.c.eff.upper.bis3 <- quantile(int.sim3, probs = c(0.65))
bayes.c.eff.lower.bis4 <- quantile(int.sim3, probs = c(0.45))
bayes.c.eff.upper.bis4 <- quantile(int.sim3, probs = c(0.55))
plot.dat.new.2 <- data.frame(bayes.c.eff.mean,
                             bayes.c.eff.lower, bayes.c.eff.upper,
                             bayes.c.eff.lower.bis, bayes.c.eff.upper.bis,
                             bayes.c.eff.lower.bis2, bayes.c.eff.upper.bis2,
                             bayes.c.eff.lower.bis3, bayes.c.eff.upper.bis3,
                             bayes.c.eff.lower.bis4, bayes.c.eff.upper.bis4)


plot.dat <- rbind(plot.dat, plot.dat.new, plot.dat.new.2)
sp <- as.data.frame(matrix(c("R", "PM", "FM"), nrow = 3, ncol = 1))
plot.dat <- cbind(plot.dat, sp)



plot_mig <- ggplot(data_meta_regression, aes(x = migration, y = abs_sens)) +
  xlab("")+
  ylab(expression(paste("Absolute rainfall sensitivity (log ", lambda, " change per 100mm)")))+
  ylim(c(0.05, 0.2))+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.title.y=element_text(size=12))

plot_mig <- plot_mig + annotate("rect", xmin = 0.7, xmax = 1.3, ymin =  plot.dat$bayes.c.eff.lower.bis4[1], ymax = plot.dat$bayes.c.eff.upper.bis4[1], alpha=0.15, fill="cornflowerblue")
plot_mig <- plot_mig + annotate("rect", xmin = 1.7, xmax = 2.3, ymin =  plot.dat$bayes.c.eff.lower.bis4[2], ymax = plot.dat$bayes.c.eff.upper.bis4[2], alpha=0.15, fill="cornflowerblue")
plot_mig <- plot_mig + annotate("rect", xmin = 2.7, xmax = 3.3, ymin =  plot.dat$bayes.c.eff.lower.bis4[3], ymax = plot.dat$bayes.c.eff.upper.bis4[3], alpha=0.15, fill="cornflowerblue")

plot_mig <- plot_mig + annotate("rect", xmin = 0.7, xmax = 1.3, ymin =  plot.dat$bayes.c.eff.lower.bis3[1], ymax = plot.dat$bayes.c.eff.upper.bis3[1], alpha=0.15, fill="cornflowerblue")
plot_mig <- plot_mig + annotate("rect", xmin = 1.7, xmax = 2.3, ymin =  plot.dat$bayes.c.eff.lower.bis3[2], ymax = plot.dat$bayes.c.eff.upper.bis3[2], alpha=0.15, fill="cornflowerblue")
plot_mig <- plot_mig + annotate("rect", xmin = 2.7, xmax = 3.3, ymin =  plot.dat$bayes.c.eff.lower.bis3[3], ymax = plot.dat$bayes.c.eff.upper.bis3[3], alpha=0.15, fill="cornflowerblue")

plot_mig <- plot_mig + annotate("rect", xmin = 0.7, xmax = 1.3, ymin =  plot.dat$bayes.c.eff.lower.bis2[1], ymax = plot.dat$bayes.c.eff.upper.bis2[1], alpha=0.15, fill="cornflowerblue")
plot_mig <- plot_mig + annotate("rect", xmin = 1.7, xmax = 2.3, ymin =  plot.dat$bayes.c.eff.lower.bis2[2], ymax = plot.dat$bayes.c.eff.upper.bis2[2], alpha=0.15, fill="cornflowerblue")
plot_mig <- plot_mig + annotate("rect", xmin = 2.7, xmax = 3.3, ymin =  plot.dat$bayes.c.eff.lower.bis2[3], ymax = plot.dat$bayes.c.eff.upper.bis2[3], alpha=0.15, fill="cornflowerblue")

plot_mig <- plot_mig + annotate("rect", xmin = 0.7, xmax = 1.3, ymin =  plot.dat$bayes.c.eff.lower.bis[1], ymax = plot.dat$bayes.c.eff.upper.bis[1], alpha=0.15, fill="cornflowerblue")
plot_mig <- plot_mig + annotate("rect", xmin = 1.7, xmax = 2.3, ymin =  plot.dat$bayes.c.eff.lower.bis[2], ymax = plot.dat$bayes.c.eff.upper.bis[2], alpha=0.15, fill="cornflowerblue")
plot_mig <- plot_mig + annotate("rect", xmin = 2.7, xmax = 3.3, ymin =  plot.dat$bayes.c.eff.lower.bis[3], ymax = plot.dat$bayes.c.eff.upper.bis[3], alpha=0.15, fill="cornflowerblue")

plot_mig <- plot_mig + annotate("rect", xmin = 0.7, xmax = 1.3, ymin =  plot.dat$bayes.c.eff.lower[1], ymax = plot.dat$bayes.c.eff.upper[1], alpha=0.15, fill="cornflowerblue")
plot_mig <- plot_mig + annotate("rect", xmin = 1.7, xmax = 2.3, ymin =  plot.dat$bayes.c.eff.lower[2], ymax = plot.dat$bayes.c.eff.upper[2], alpha=0.15, fill="cornflowerblue")
plot_mig <- plot_mig + annotate("rect", xmin = 2.7, xmax = 3.3, ymin =  plot.dat$bayes.c.eff.lower[3], ymax = plot.dat$bayes.c.eff.upper[3], alpha=0.15, fill="cornflowerblue")

plot_mig <- plot_mig + annotate("segment", x = 0.7, xend = 1.3, y = plot.dat$bayes.c.eff.mean[1], yend = plot.dat$bayes.c.eff.mean[1], size=1.8, alpha=0.7, colour="cornflowerblue")
plot_mig <- plot_mig + annotate("segment", x = 1.7, xend = 2.3, y = plot.dat$bayes.c.eff.mean[2], yend = plot.dat$bayes.c.eff.mean[2], size=1.8, alpha=0.7, colour="cornflowerblue")
plot_mig <- plot_mig + annotate("segment", x = 2.7, xend = 3.3, y = plot.dat$bayes.c.eff.mean[3], yend = plot.dat$bayes.c.eff.mean[3], size=1.8, alpha=0.7, colour="cornflowerblue")

plot_mig <- plot_mig + scale_x_continuous(breaks = c(1,2,3), labels = c("Sedentary", "Partially migratory", "Migratory"))
plot_mig




posterior_mig <- as.data.frame(dat_plot$bsp_momigration)

plot_posterior_mig <- ggplot(posterior_mig, aes(x=`dat_plot$bsp_momigration`))+
  stat_slab(aes(fill = after_stat(level), alpha = after_stat(level)), .width = c(.1,.3,.5,.7,.9, 1)) +
  scale_fill_manual(values = c("cornflowerblue", "cornflowerblue", "cornflowerblue", "cornflowerblue", "cornflowerblue", "cornflowerblue"))+
  geom_vline(aes(xintercept=0), color="grey10", linewidth=1)+
  xlab("Posterior estimate")+
  theme_bw()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  theme(legend.position = "none")+
  theme(axis.line.x = element_line(colour = "black"),
        axis.text.x = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank())
plot_posterior_mig

plot_mig_final <- plot_mig + inset_element(plot_posterior_mig, 0.5, 0.75, 0.95, 0.95)
plot_mig_final




plot_coop_final /
  (plot_latitude_final + plot_mass_final) /
  (plot_hwi_final + plot_mig_final)


