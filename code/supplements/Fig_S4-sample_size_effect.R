rm (list = ls ( ) ) 

# Load packages

library(ggplot2)
library(tidyverse)
library(ggdist)
library(patchwork)

setwd("")

# Import dataset
population_sensitivity_temp <- read.csv("population_sensitivity_temperature.txt", sep="")  # output of script 3
population_sensitivity_rain <- read.csv("population_sensitivity_rainfall.txt", sep="")  # output of script 4

population_sensitivity_temp$abs_sens <- abs(population_sensitivity_temp$mean)
population_sensitivity_rain$abs_sens <- abs(population_sensitivity_rain$mean)


plot_year_temp <- ggplot(population_sensitivity_temp, aes(x=sample_size, y=abs_sens))+
  geom_jitter(alpha = 0.2, width = 0.2, height = 0)+
  xlab("Number of years")+
  ylab(expression(paste("Absolute temperature sensitivity (log ", lambda, " change per degree C)")))+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.title.y=element_text(size=12))

plot_year_rain <- ggplot(population_sensitivity_rain, aes(x=sample_size, y=abs_sens))+
  geom_jitter(alpha = 0.2, width = 0.2, height = 0)+
  xlab("Number of years")+
  ylab(expression(paste("Absolute rainfall sensitivity (log ", lambda, " change per 100mm)")))+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.title.y=element_text(size=12))

figure_s4 <- plot_year_temp / plot_year_rain
