############################################################################   
############################################################################
##                                                                        ##
##  R-code: Combine datasets to create final data for meta-regression     ##
##                                                                        ##
############################################################################
############################################################################

rm (list = ls ( ) ) 

setwd("")

########################
#### Import dataset ####
########################

population_sensitivity_temperature <- read.csv("population_sensitivity_temperature.txt", sep="") # output of script 3
population_sensitivity_rainfall <- read.csv("population_sensitivity_rainfall.txt", sep="")  # output of script 4

coop <- read.csv("sociality_data/cooperative_breeding_data.txt", sep="")
avonet <- read.csv("AVONET.csv")

latitude_data <- read.csv("growth_rate_climate.txt", sep="")  # output of script 2

# keep only columns that are needed
population_sensitivity_temperature <- subset(population_sensitivity_temperature, select = c(V2, 
                                                                                            mean,
                                                                                            sample_size,
                                                                                            Binomial))
colnames(population_sensitivity_temperature) <- c("id_block", "slope_temperature", "sample_size", "species")

population_sensitivity_rainfall <- subset(population_sensitivity_rainfall, select = c(V2, 
                                                                                      mean))
colnames(population_sensitivity_rainfall) <- c("id_block", "slope_rainfall")


# merge sensitivity to rainfall and temperature into a single dataset
combined_sensitivity <- merge(population_sensitivity_temperature, population_sensitivity_rainfall, all.x = T, by.x = "id_block", by.y = "id_block")

combined_sensitivity <- combined_sensitivity[,c(4,1,3,2,5)]

# add absolute values of sensitivity
combined_sensitivity$abs_slope_temperature <- abs(combined_sensitivity$slope_temperature)
combined_sensitivity$abs_slope_rainfall <- abs(combined_sensitivity$slope_rainfall)


############################
#### add sociality data ####
############################

coop <- subset(coop, select = c(tip_label,
                                fam_sys_known50,
                                fam_sys_inferred50))

combined_sensitivity <- merge(combined_sensitivity, coop, all.x = T, by.x = "species", by.y = "tip_label")
combined_sensitivity$coop_num <- NA
combined_sensitivity$coop_num[combined_sensitivity$fam_sys_inferred50=="no_fam"] <- 1
combined_sensitivity$coop_num[combined_sensitivity$fam_sys_inferred50=="family"] <- 2
combined_sensitivity$coop_num[combined_sensitivity$fam_sys_inferred50=="coop_families"] <- 3

combined_sensitivity$fam_sys_inferred50[combined_sensitivity$fam_sys_inferred50=="no_fam"] <- "1_no_fam"
combined_sensitivity$fam_sys_inferred50[combined_sensitivity$fam_sys_inferred50=="family"] <- "2_family"
combined_sensitivity$fam_sys_inferred50[combined_sensitivity$fam_sys_inferred50=="coop_families"] <- "3_coop_families"


#################################################
#### add life history data (e.g., mass, HWI) ####
#################################################

avonet <- subset(avonet, select = c(Species2,
                                    Hand.Wing.Index,
                                    Migration,
                                    Mass))

# first, let's make sure names match
avonet$Species2 <- sub(" ", "_", avonet$Species2)

setdiff(combined_sensitivity$species, avonet$Species2)

# correct names for unmatching taxonomy
avonet$Species2[avonet$Species2=="Charadrius_alexandrinus"] <- "Anarhynchus_alexandrinus"
avonet$Species2[avonet$Species2=="Charadrius_leschenaultii"] <- "Anarhynchus_leschenaultii"
avonet$Species2[avonet$Species2=="Charadrius_mongolus"] <- "Anarhynchus_mongolus"
avonet$Species2[avonet$Species2=="Charadrius_pecuarius"] <- "Anarhynchus_pecuarius"
avonet$Species2[avonet$Species2=="Bubulcus_ibis"] <- "Ardea_ibis"
avonet$Species2[avonet$Species2=="Ardea_intermedia"] <- "Ardea_plumifera"
avonet$Species2[avonet$Species2=="Accipiter_gentilis"] <- "Astur_gentilis"
avonet$Species2[avonet$Species2=="Ixobrychus_minutus"] <- "Botaurus_minutus"
avonet$Species2[avonet$Species2=="Ixobrychus_sinensis"] <- "Botaurus_sinensis"
avonet$Species2[avonet$Species2=="Lophochroa_leadbeateri"] <- "Cacatua_leadbeateri"
avonet$Species2[avonet$Species2=="Corvus_monedula"] <- "Coloeus_monedula"
avonet$Species2[avonet$Species2=="Charadrius_morinellus"] <- "Eudromias_morinellus"
avonet$Species2[avonet$Species2=="Haliaeetus_vocifer"] <- "Icthyophaga_vocifer"
avonet$Species2[avonet$Species2=="Dryobates_arizonae"] <- "Leuconotopicus_arizonae"
avonet$Species2[avonet$Species2=="Alopecoenas_xanthonurus"] <- "Pampusana_xanthonura"
avonet$Species2[avonet$Species2=="Streptopelia_chinensis"] <- "Spilopelia_chinensis"
avonet$Species2[avonet$Species2=="Apus_melba"] <- "Tachymarptis_melba"
avonet$Species2[avonet$Species2=="Thalassarche_chlororhynchos"] <- "Thalassarche_carteri"
avonet$Species2[avonet$Species2=="Calyptorhynchus_latirostris"] <- "Zanda_latirostris"

# add duplicate rows based on data from larus argentatus and diomedea exulans because they are considered subspecies in avonet
avonet <- rbind(avonet, c("Larus_smithsonianus", 54.9, 2, 1091.0), c("Diomedea_amsterdamensis", 56.4, 2, 6961.3))
avonet$Mass <- as.numeric(avonet$Mass)
avonet$Hand.Wing.Index <- as.numeric(avonet$Hand.Wing.Index)
avonet$Migration <- as.numeric(avonet$Migration)

setdiff(combined_sensitivity$species, avonet$Species2)

colnames(avonet) <- c("species", "hand_wing_index", "migration", "mass")

# merge
combined_sensitivity <- merge(combined_sensitivity, avonet, all.x = T, by.x = "species", by.y = "species")

combined_sensitivity$logmass <- log(combined_sensitivity$mass)


######################
#### Add latitude ####
######################

latitude_data <- subset(latitude_data, select = c("ID_block", "Latitude"))
latitude_data <- latitude_data[!duplicated(latitude_data$ID_block),]
colnames(latitude_data) <- c("ID_block", "latitude")

combined_sensitivity <- merge(combined_sensitivity, latitude_data, all.x = T, by.x = "id_block", by.y = "ID_block")
combined_sensitivity$abs_latitude <- abs(combined_sensitivity$latitude)



############################
#### Save final dataset ####
############################

getwd()
#write.table(combined_sensitivity, "final_dataset.txt")

