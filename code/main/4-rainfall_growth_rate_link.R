############################################################################   
############################################################################
##                                                                        ##
##  R-code: estimate the effect of the environment on lambda for each     ##
##          sample                                                        ##
##                                                                        ##
##  This file estimate the effect of rainfall on population growth        ##
##  rate for each samples in our dataset, while accounting for temporal   ##
##  autocorrelation and non-linear temporal effects                       ##
##                                                                        ##
############################################################################
############################################################################

#####################################################
#### clean working environment and load packages ####
#####################################################

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

######################
#### Load dataset ####
######################

setwd("")

# Import dataset
bird_data <- read.csv("growth_rate_climate.txt", sep="")  # output of script 2

# Sort dataset by time before present
bird_data <- bird_data[order(bird_data$year, decreasing = F),]

# delta precipitation in 100mm unit
bird_data$delta_precipitation_100 <- bird_data$delta_precipitation/100 

# Get the list of samples
# list_of_samples <- unique(data_for_first_analysis$name)
list_of_samples <- unique(bird_data$ID_block)

#remove samples with no climate data
samples_to_remove <- unique(bird_data$ID_block[is.na(bird_data$delta_precipitation)])
list_of_samples <- list_of_samples[!c(list_of_samples %in% samples_to_remove)]

#create empty lists to store all the model summary and posterior samples
store_posterior_effect_size <- list()
store_model_summary <- list()
store_warnings <- list()
list_pp_check <- list()
store_n_years <- list()

quiet_brm <- quietly(brm) # to be able to save warnings

tic()
## Loop over all samples, to perform a model for each sample
for (i in 1:length(list_of_samples)) {
  
  # get dataset for a given sample
  temporary_dataset <- subset(bird_data, ID_block==list_of_samples[i])
  rownames(temporary_dataset) <- NULL
  
  store_n_years[[i]] <-  as.data.frame(cbind(list_of_samples[i], dim(temporary_dataset)[1]))
  
  # Model, including spline for non-linear temporal effects, and autoregressive (ar) model for temporal autocorrelation
  model_spline_ar <- quiet_brm(
    pop_growth_rate ~ delta_precipitation_100 +
      s(year, bs="tp", k=5) +
      ar(time=year)
    ,
    data = temporary_dataset, cores = 3, chains = 3,
    warmup = 1000, iter = 2000, thin = 1, seed = 123,
    family = gaussian(),
    prior = c(
      set_prior("normal(0,0.5)", class = "ar"),
      set_prior("exponential(2)", class = "sds"),
      set_prior("normal(0,1)", class = "b", coef = "delta_precipitation_100"),
      set_prior("normal(0,1)", class = "Intercept")),
    control = list(adapt_delta = 0.999),
    backend = "cmdstanr"
  )
  
  store_warnings[[i]] <- model_spline_ar$messages
  
  list_pp_check[[i]] <- pp_check(model_spline_ar$result, ndraws = 100) + theme_bw()
  
  # save posterior samples for the effect of covariate on population growth rate
  post_model_spline <- as.data.frame(model_spline_ar$result)
  store_posterior_effect_size[[i]] <- as.data.frame(
    cbind(as.numeric(post_model_spline$b_delta_precipitation_100),
          rep(list_of_samples[i], 3000)))
  
  # save model summary
  store_model_summary[[i]] <- summary(model_spline_ar$result)
  
  print(i)
}
toc()

# 
all_effect_sizes <- rbindlist(store_posterior_effect_size)
all_effect_sizes$V1 <- as.numeric(all_effect_sizes$V1)

#scale effect sizes
all_effect_sizes$scaled_beta <- scale(all_effect_sizes$V1, center = F)

# summarise results for all samples
result_summary <- all_effect_sizes %>%
  group_by(V2) %>% 
  summarise(
    quantile_low = quantile(V1, 0.025),
    quantile_low_89 = quantile(V1, 0.055),
    quantile_low_50 = quantile(V1, 0.25),
    mean = mean(V1),
    quantile_high = quantile(V1, 0.975),
    quantile_high_89 = quantile(V1, 0.945),
    quantile_high_50 = quantile(V1, 0.75),
    quantile_low_scaled = quantile(scaled_beta, 0.025),
    mean_scaled = mean(scaled_beta),
    quantile_high_scaled = quantile(scaled_beta, 0.975))

# Estimate standard error for each effect size estimate
result_summary$se <- (result_summary$quantile_high - result_summary$quantile_low) / 3.92
result_summary$se_scaled <- (result_summary$quantile_high_scaled - result_summary$quantile_low_scaled) / 3.92

# add sample size
n_years_sample_size <- rbindlist(store_n_years)
colnames(n_years_sample_size) <- c("ID_block", "sample_size")
result_summary <- merge(result_summary, n_years_sample_size, all.x = T, by.x = "V2", by.y = "ID_block")

# create a ranked summary, to be used for forest plot
result_summary_sorted <- result_summary[order(result_summary$mean_scaled, decreasing = F),]
rownames(result_summary_sorted) <- NULL
result_summary_sorted$rank <- as.numeric(rownames(result_summary_sorted))

# forest plot 1
ggplot(result_summary_sorted, aes(x = mean, xmin = quantile_low, xmax = quantile_high, y = rank)) +
  geom_pointrange(size = 0.1, color = "black") +
  theme_bw()


# add species name to the summary of results
data_for_meta_regression <- merge(result_summary, unique(bird_data[,c(1,3)]), all.x = F, all.y = F, by.x = "V2", by.y = "ID_block")

# fix species names to match phylogeny

data_for_meta_regression$Binomial[data_for_meta_regression$Binomial=="Leiopicus_medius"] <- "Dendrocoptes_medius"
data_for_meta_regression$Binomial[data_for_meta_regression$Binomial=="Bubulcus_ibis"] <- "Ardea_ibis"
data_for_meta_regression$Binomial[data_for_meta_regression$Binomial=="Catharacta_skua"] <- "Stercorarius_skua"
data_for_meta_regression$Binomial[data_for_meta_regression$Binomial=="Hypotaenidia_sylvestris"] <- "Gallirallus_sylvestris"
data_for_meta_regression$Binomial[data_for_meta_regression$Binomial=="Charadrius_mongolus"] <- "Anarhynchus_mongolus"
data_for_meta_regression$Binomial[data_for_meta_regression$Binomial=="Charadrius_leschenaultii"] <- "Anarhynchus_leschenaultii"
data_for_meta_regression$Binomial[data_for_meta_regression$Binomial=="Larus_audouinii"] <- "Ichthyaetus_audouinii"
data_for_meta_regression$Binomial[data_for_meta_regression$Binomial=="Larus_ridibundus"] <- "Chroicocephalus_ridibundus"
data_for_meta_regression$Binomial[data_for_meta_regression$Binomial=="Cyanecula_svecica"] <- "Luscinia_svecica"
data_for_meta_regression$Binomial[data_for_meta_regression$Binomial=="Larus_genei"] <- "Chroicocephalus_genei"
data_for_meta_regression$Binomial[data_for_meta_regression$Binomial=="Charadrius_alexandrinus"] <- "Anarhynchus_alexandrinus"
data_for_meta_regression$Binomial[data_for_meta_regression$Binomial=="Accipiter_gentilis"] <- "Astur_gentilis"
data_for_meta_regression$Binomial[data_for_meta_regression$Binomial=="Bonasa_bonasia"] <- "Tetrastes_bonasia"
data_for_meta_regression$Binomial[data_for_meta_regression$Binomial=="Sylvia_melanocephala"] <- "Curruca_melanocephala"
data_for_meta_regression$Binomial[data_for_meta_regression$Binomial=="Sylvia_communis"] <- "Curruca_communis"
data_for_meta_regression$Binomial[data_for_meta_regression$Binomial=="Corvus_monedula"] <- "Coloeus_monedula"
data_for_meta_regression$Binomial[data_for_meta_regression$Binomial=="Catharacta_maccormicki"] <- "Stercorarius_maccormicki"
data_for_meta_regression$Binomial[data_for_meta_regression$Binomial=="Larus_melanocephalus"] <- "Ichthyaetus_melanocephalus"
data_for_meta_regression$Binomial[data_for_meta_regression$Binomial=="Larus_novaehollandiae"] <- "Chroicocephalus_novaehollandiae"
data_for_meta_regression$Binomial[data_for_meta_regression$Binomial=="Ixobrychus_sinensis"] <- "Botaurus_sinensis"
data_for_meta_regression$Binomial[data_for_meta_regression$Binomial=="Sylvia_cantillans"] <- "Curruca_cantillans"
data_for_meta_regression$Binomial[data_for_meta_regression$Binomial=="Sylvia_curruca"] <- "Curruca_curruca"
data_for_meta_regression$Binomial[data_for_meta_regression$Binomial=="Sylvia_hortensis"] <- "Curruca_hortensis"
data_for_meta_regression$Binomial[data_for_meta_regression$Binomial=="Sylvia_nisoria"] <- "Curruca_nisoria"
data_for_meta_regression$Binomial[data_for_meta_regression$Binomial=="Sylvia_undata"] <- "Curruca_undata"
data_for_meta_regression$Binomial[data_for_meta_regression$Binomial=="Sylvia_melanothorax"] <- "Curruca_melanothorax"
data_for_meta_regression$Binomial[data_for_meta_regression$Binomial=="Larus_ichthyaetus"] <- "Ichthyaetus_ichthyaetus"
data_for_meta_regression$Binomial[data_for_meta_regression$Binomial=="Nannopterum_auritus"] <- "Nannopterum_auritum"
data_for_meta_regression$Binomial[data_for_meta_regression$Binomial=="Manucerthia_mana"] <- "Loxops_mana"
data_for_meta_regression$Binomial[data_for_meta_regression$Binomial=="Haliaeetus_vocifer"] <- "Icthyophaga_vocifer"
data_for_meta_regression$Binomial[data_for_meta_regression$Binomial=="Ixobrychus_minutus"] <- "Botaurus_minutus"
data_for_meta_regression$Binomial[data_for_meta_regression$Binomial=="Charadrius_pecuarius"] <- "Anarhynchus_pecuarius"
data_for_meta_regression$Binomial[data_for_meta_regression$Binomial=="Taeniopygia_bichenovii"] <- "Stizoptera_bichenovii"
data_for_meta_regression$Binomial[data_for_meta_regression$Binomial=="Larus_maculipennis"] <- "Chroicocephalus_maculipennis"
data_for_meta_regression$Binomial[data_for_meta_regression$Binomial=="Alopecoenas_xanthonurus"] <- "Pampusana_xanthonura"
data_for_meta_regression$Binomial[data_for_meta_regression$Binomial=="Threskiornis_moluccus"] <- "Threskiornis_molucca"


###################
#### Save data ####
###################

getwd()
#write.table(data_for_meta_regression, "population_sensitivity_rainfall.txt")


###################
#### Figure S3 ####
###################

figure_s3 <- ggplot(result_summary_sorted, aes(y = rank)) +   
  
  geom_linerange(data = filter(result_summary_sorted, quantile_high_89 > 0) |>
                   transmute(rank, quantile_low_89 = pmax(0, quantile_low_89), quantile_high_89),
                 aes(xmin = quantile_low_89, xmax = quantile_high_89),
                 color = "#FC9272", linewidth = 0.15) +
  
  geom_linerange(data = filter(result_summary_sorted, quantile_low_89 < 0) |>
                   transmute(rank, quantile_low_89, quantile_high_89 = pmin(0, quantile_high_89)),
                 aes(xmin = quantile_low_89, xmax = quantile_high_89),
                 color = "#9ECAE1", linewidth = 0.15) +
  
  geom_linerange(data = filter(result_summary_sorted, quantile_high_50 > 0) |>
                   transmute(rank, quantile_low_50 = pmax(0, quantile_low_50), quantile_high_50),
                 aes(xmin = quantile_low_50, xmax = quantile_high_50),
                 color = "#CB181D", linewidth = 0.15) +
  
  geom_linerange(data = filter(result_summary_sorted, quantile_low_50 < 0) |>
                   transmute(rank, quantile_low_50, quantile_high_50 = pmin(0, quantile_high_50)),
                 aes(xmin = quantile_low_50, xmax = quantile_high_50),
                 color = "#2171B5", linewidth = 0.15) +
  
  geom_point(aes(x = mean), size = 0.04) +
  scale_x_continuous(breaks = seq(-1, 1, by = 0.5)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_vline(xintercept = 0, colour = "black", linewidth = 0.5) +
  xlab(expression(paste("Effect size of precipitation change on log ", lambda, " (per 100mm)")))+
  coord_cartesian(xlim = c(-1.2, 1.2), ylim = c(-5, 2383)) +
  theme_minimal() +
  theme(plot.title = element_text(size=10, hjust = 0.5),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 0.5, fill = NA),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=16),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank())

