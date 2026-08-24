############################################################################   
############################################################################
##                                                                        ##
##  R-code: Meta-regression to estimate the effect of traits on           ##
##          the demographic sensitivity of populations to temperature     ##
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

##################
#### 20 years ####
##################

model_20years <- brm(
  abs_slope_temperature ~ 1 + 
    scale(logmass) +
    scale(hand_wing_index) +
    scale(sample_size) + 
    scale(abs_latitude) +
    mo(coop_num) +
    mo(migration) +
    (1 | species) +
    (1 | gr(phylo, cov = phylo_mcct))
  ,
  data = data_meta_regression[data_meta_regression$sample_size>=20,],
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

summary(model_20years, prob=0.89)

pp_check(model_20years, ndraws = 100) + xlim(c(0,2)) + theme_bw()

posterior_20y <- as.matrix(model_20years)
dat_plot <- as.data.frame(posterior_20y)



#################################
#### Plot posterior 20 years ####
#################################

posterior_coop <- as.data.frame(dat_plot$bsp_mocoop_num)

plot_posterior_coop_20y <- ggplot(posterior_coop, aes(x=`dat_plot$bsp_mocoop_num`))+
  stat_slab(aes(fill = after_stat(level), alpha = after_stat(level)), .width = c(.1,.3,.5,.7,.9, 1)) +
  scale_fill_manual(values = c("tomato3", "tomato3", "tomato3", "tomato3", "tomato3", "tomato3"))+
  geom_vline(aes(xintercept=0), color="grey10", linewidth=1)+
  xlab("Posterior estimate")+
  xlim(c(-0.45, 0.25))+
  ggtitle("Time series > 20 years")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  theme(legend.position = "none")+
  theme(axis.line.x = element_line(colour = "black"),
        axis.text.x = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank())
plot_posterior_coop_20y




##################
#### 15 years ####
##################

model_15years <- brm(
  abs_slope_temperature ~ 1 + 
    scale(logmass) +
    scale(hand_wing_index) +
    scale(sample_size) + 
    scale(abs_latitude) +
    mo(coop_num) +
    mo(migration) +
    (1 | species) +
    (1 | gr(phylo, cov = phylo_mcct))
  ,
  data = data_meta_regression[data_meta_regression$sample_size>=15,],
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

summary(model_15years, prob=0.89)

pp_check(model_15years, ndraws = 100) + xlim(c(0,2)) + theme_bw()

posterior_15y <- as.matrix(model_15years)
dat_plot <- as.data.frame(posterior_15y)



#################################
#### Plot posterior 15 years ####
#################################

posterior_coop <- as.data.frame(dat_plot$bsp_mocoop_num)

plot_posterior_coop_15y <- ggplot(posterior_coop, aes(x=`dat_plot$bsp_mocoop_num`))+
  stat_slab(aes(fill = after_stat(level), alpha = after_stat(level)), .width = c(.1,.3,.5,.7,.9, 1)) +
  scale_fill_manual(values = c("tomato3", "tomato3", "tomato3", "tomato3", "tomato3", "tomato3"))+
  geom_vline(aes(xintercept=0), color="grey10", linewidth=1)+
  xlab("Posterior estimate")+
  ggtitle("Time series > 15 years")+
  xlim(c(-0.45, 0.25))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  theme(legend.position = "none")+
  theme(axis.line.x = element_line(colour = "black"),
        axis.text.x = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank())
plot_posterior_coop_15y




##################
#### 10 years ####
##################

model_10years <- brm(
  abs_slope_temperature ~ 1 + 
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

summary(model_10years, prob=0.89)

pp_check(model_10years, ndraws = 100) + xlim(c(0,2)) + theme_bw()

posterior_10y <- as.matrix(model_10years)
dat_plot <- as.data.frame(posterior_10y)



#################################
#### Plot posterior 10 years ####
#################################

posterior_coop <- as.data.frame(dat_plot$bsp_mocoop_num)

plot_posterior_coop_10y <- ggplot(posterior_coop, aes(x=`dat_plot$bsp_mocoop_num`))+
  stat_slab(aes(fill = after_stat(level), alpha = after_stat(level)), .width = c(.1,.3,.5,.7,.9, 1)) +
  scale_fill_manual(values = c("tomato3", "tomato3", "tomato3", "tomato3", "tomato3", "tomato3"))+
  geom_vline(aes(xintercept=0), color="grey10", linewidth=1)+
  xlab("Posterior estimate")+
  ggtitle("Time series > 10 years")+
  xlim(c(-0.45, 0.25))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  theme(legend.position = "none")+
  theme(axis.line.x = element_line(colour = "black"),
        axis.text.x = element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank())
plot_posterior_coop_10y


figure_s5 <- plot_posterior_coop_10y / plot_posterior_coop_15y / plot_posterior_coop_20y


