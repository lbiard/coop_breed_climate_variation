############################################################################   
############################################################################
##                                                                        ##
##  R-code: Clean Living Planet data, selecting only time series meeting  ##
##          the requirements                                              ##
##                                                                        ##
##  This file use the whole Living Planet Database as an input, which is  ##
##  publicly available at https://www.livingplanetindex.org/data_portal   ##
##  keeping only suitable time series and estimating pop growth rate      ##
##                                                                        ##
############################################################################
############################################################################


#####################################################
#### clean working environment and load packages ####
#####################################################

rm(list = ls())

library(ggplot2)
library(brms)
library(data.table)
library(tidyverse)
library(patchwork)

######################
#### Load dataset ####
######################

setwd("")

lpdata <- read.csv("data/LPD_2024_public.csv")

###########################
#### Subset and format ####
###########################

# keep only birds
lpd_bird <- subset(lpdata, Class == "Aves")

# remove some unnecessary columns
lpd_bird <- subset(lpd_bird, select = -c(X, Included.in.LPR2024, Subspecies,
                                          Country, All_countries, Region, IPBES_region,
                                          IPBES_subregion, T_realm, T_biome, FW_realm,
                                          FW_biome, M_realm, M_ocean, M_biome, Native))

# exclude replicates
lpd_bird <- subset(lpd_bird, Replicate == 0)

# wide to long
lpd_bird <- melt(setDT(lpd_bird), measure.vars = colnames(lpd_bird)[17:87], variable.name = c("year"))

# remove data with NULL values (years with no observation)
lpd_bird <- subset(lpd_bird, value != "NULL")

# exclude total absence data (leads to issues when calculating growth rates)
lpd_bird <- subset(lpd_bird, value != 0)

# year as numerical
lpd_bird$year <- substr(lpd_bird$year, 2, 5)
lpd_bird$year <- as.numeric(lpd_bird$year)

# abundance as numerical
lpd_bird$value <- as.numeric(lpd_bird$value)

# check for record that are in log abundance, change to normal scale
unique(lpd_bird$Units)

lpd_bird$value[lpd_bird$ID=="23374"] <- exp(lpd_bird$value[lpd_bird$ID=="23374"])
lpd_bird$Units[lpd_bird$ID=="23374"] <- "mean colony counts"

lpd_bird$value[lpd_bird$ID=="6447"] <- exp(lpd_bird$value[lpd_bird$ID=="6447"])
lpd_bird$Units[lpd_bird$ID=="6447"] <- "numbers"

lpd_bird$value[lpd_bird$ID=="306"] <- exp(lpd_bird$value[lpd_bird$ID=="306"])
lpd_bird$Units[lpd_bird$ID=="306"] <- "number of nestlings"

lpd_bird$value[lpd_bird$ID=="308"] <- exp(lpd_bird$value[lpd_bird$ID=="308"])
lpd_bird$Units[lpd_bird$ID=="308"] <- "number of nestlings"


######################################################
#### Keep only data with enough consecutive years ####
######################################################

# Keeping only years from 1950 onwards

lpd_bird <- lpd_bird %>% 
  filter(year >= 1949) %>% 
  group_by(ID) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  filter(n >= 5) %>% 
  dplyr::select(-n)

# gaps in abundance
lpd_bird_gaps <- lpd_bird %>% 
  group_by(ID) %>% 
  group_modify(~{
    cyears = .$year
    diff_cyears = diff(cyears)
    cumsum_blocks = cumsum(c(1, diff_cyears != 1))
    
    summarise(., Binomial = Binomial[1],
              record_length = length(cyears),
              no_consecutive_blocks = n_distinct(cumsum_blocks),
              prop_1year_transitions = sum(diff_cyears == 1)/ length(diff_cyears),
              longest_block = max(table(cumsum_blocks)))
  }) %>% 
  ungroup()

# consecutive blocks
lpd_bird_blocks <- lpd_bird %>% 
  group_by(ID) %>%
  mutate(block = cumsum(c(1, diff(year) != 1)),
         max_block = max(block)) %>% 
  ungroup() %>% 
  dplyr::select(ID, Binomial, Order, value, 
                year, block, max_block) %>% 
  left_join(x = ., y = dplyr::select(lpd_bird_gaps, -c(Binomial, no_consecutive_blocks)),
            by = "ID") %>% 
  arrange(desc(longest_block)) %>% 
  mutate(ID = factor(ID, levels = unique(.$ID)))


#### Long blocks of 10 years or more ####

# IDs and blocks that we want to keep
ID_block_keep_long <- lpd_bird_blocks %>% 
  mutate(ID = as.numeric(as.character(ID))) %>% 
  group_by(ID, block) %>% 
  summarise(ID_block = paste0(ID[1],"_",block[1]),
            block_keep = if_else(n() >= 10, 1, 0)) %>% 
  ungroup() %>% 
  filter(block_keep == 1)

# Restricting the dataset
lpd_bird_IDblocks_10yr <- lpd_bird %>% 
  group_by(ID) %>%
  mutate(block = cumsum(c(1, diff(year) != 1)),
         ID_block = paste0(ID[1],"_",block)) %>% 
  ungroup() %>% 
  filter(ID_block %in% ID_block_keep_long$ID_block == T) 


# Summarize data 
lpd_bird_datasum10 <- data.frame(Dataset = c("Raw data", "Study data >= 10 years"),
                            Observations  = c(nrow(lpd_bird), 
                                              nrow(lpd_bird_IDblocks_10yr)),
                            Records = c(n_distinct(lpd_bird$ID), 
                                        n_distinct(lpd_bird_IDblocks_10yr$ID)),
                            Species = c(n_distinct(lpd_bird$Binomial),
                                        n_distinct(lpd_bird_IDblocks_10yr$Binomial)))


#### Calculate annual population growth rates ####

bird10 <- lpd_bird_IDblocks_10yr %>% 
  mutate(value2 = value) %>% 
  group_by(ID_block) %>% 
  group_modify(~{
    t0 <- .$value2[-(length(.$value2))] # get rid the last obs
    t1 <- .$value2[-1]                        # get rid of the first obs
    
    mutate(., pop_growth_rate = c(log(t1/t0),NA))
  }) %>% 
  ungroup() %>% 
  filter(is.na(pop_growth_rate) == F)

###################
#### Save data ####
###################

getwd()
#write.table(bird10, "data/growth_rate_populations.txt")


###################
#### Figure S1 ####
###################

# study length per record
n_records_10 <-  bird10 %>% 
  group_by(ID_block) %>% 
  summarise(study_length = n()) 

plot_record_length <- ggplot(n_records_10, aes(x = study_length)) + 
  geom_histogram(bins = 30,fill = "violetred4") +
  labs(x = "Number of years", y = "Number of records") +
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

# number of records species
n_records_per_sp_10 <- bird10 %>% 
  group_by(Binomial, ID_block) %>% 
  summarise(study_length = n()) %>% 
  group_by(Binomial) %>% 
  summarise(n_rec = n()) 

plot_record_species <- ggplot(n_records_per_sp_10, aes(x = log(n_rec))) + 
  geom_histogram(bins = 30,fill = "violetred4") +
  labs(x = "Log number of records", y = "Number of species") +
  theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

## combine figures
figure_s1 <- plot_record_length / plot_record_species
figure_s1

