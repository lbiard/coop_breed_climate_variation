############################################################################   
############################################################################
##                                                                        ##
##  R-code: Extract the annual climatic variable for each year and each   ##
##          population                                                    ##
##                                                                        ##
##  This file extract the relevant annual temperature and precipitation   ##
##  for each population based on their latitude and longitude.            ##
##  The climatic data is publicly available on CDS                        ##
##  https://doi.org/10.24381/cds.f17050d7                                 ##
##                                                                        ##
############################################################################
############################################################################


#####################################################
#### clean working environment and load packages ####
#####################################################

rm(list = ls())

library(raster)
library(terra)
library(rasterVis)
library(tictoc)
library(exactextractr)

#######################
#### Load datasets ####
#######################

setwd("")

# Load growth rate data derived from LPI data
bird_data <- read.csv("growth_rate_populations.txt", sep="")  # output of script 1

# Load precipitation data
total_precipitation <- raster::brick("data_stream-moda_stepType-avgad.nc")
total_precipitation

# Load temperature data
mean_temperature <- raster::brick("data_stream-moda_stepType-avgua.nc")
mean_temperature

# Correct format of longitude (from -180 to 180, to 0 to 360) to match the format in the climate data rasters
bird_data$Longitude_corrected <- NA
bird_data$Longitude_corrected[bird_data$Longitude > 0] <- bird_data$Longitude[bird_data$Longitude > 0]
bird_data$Longitude_corrected[bird_data$Longitude < 0] <- 360 + bird_data$Longitude[bird_data$Longitude < 0]


bird_data$precipitation <- NA
bird_data$delta_precipitation <- NA
bird_data$temperature <- NA
bird_data$delta_temperature <- NA


#################################################################################
#### Extract climatic data from the rasters for each population in each year ####
#################################################################################

# Temperature is in Kelvin in the raster file, but I convert it to Celsius
# Rainfall is initially in a format of monthly average of daily total precipitation in meters
# which I convert into a total for the year

# extracting all the climatic data will take several hours

tic()
for (i in 1:dim(bird_data)[1]) {
  
  # total precipitation in year t
  bird_data$precipitation[i] <- 
    sum(raster::extract(total_precipitation, 
                        cbind(bird_data$Longitude_corrected[i],
                              bird_data$Latitude[i]))[1, (((bird_data$year[i] - 1950) * 12)+1) : (((bird_data$year[i] - 1950) * 12)+12)] *
        c(31,28,31,30,31,30,31,31,30,31,30,31) * 1000)
  
  
  # change in precipitation from year t to t+1
  bird_data$delta_precipitation[i] <- 
    sum(raster::extract(total_precipitation, 
                        cbind(bird_data$Longitude_corrected[i],
                              bird_data$Latitude[i]))[1, ((((bird_data$year[i]+1) - 1950) * 12)+1) : ((((bird_data$year[i]+1) - 1950) * 12)+12)] *
          c(31,28,31,30,31,30,31,31,30,31,30,31) * 1000) -
    bird_data$precipitation[i]
  
  # mean temperature in year t
  bird_data$temperature[i] <- mean(raster::extract(mean_temperature, 
                                              cbind(bird_data$Longitude_corrected[i],
                                                    bird_data$Latitude[i]))[1, (((bird_data$year[i] - 1950) * 12)+1) : (((bird_data$year[i] - 1950) * 12)+12)]) - 273.15
  
  # change in temperature from year t to t+1
  bird_data$delta_temperature[i] <- (mean(raster::extract(mean_temperature, 
                                                   cbind(bird_data$Longitude_corrected[i],
                                                         bird_data$Latitude[i]))[1, ((((bird_data$year[i]+1) - 1950) * 12)+1) : ((((bird_data$year[i]+1) - 1950) * 12)+12)]) - 273.15) -
    bird_data$temperature[i]
  
  if(i %% 100==0) {
    # track progress
    print(i)
  }
  
}
toc()


###################
#### Save data ####
###################

getwd()
#write.table(bird_data, "growth_rate_climate.txt")



