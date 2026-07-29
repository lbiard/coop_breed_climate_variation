# Populations of cooperatively breeding birds are buffered against annual variation in rainfall and temperature

This repository hosts data and code for Bliard L, Griesser M, Germain R. Populations of cooperatively breeding birds are buffered against annual variation in rainfall and temperature.

Preprint version: XXXX

## GENERAL INFORMATION

1. Title: Data and scripts from "Populations of cooperatively breeding birds are buffered against annual variation in rainfall and temperature".

2. Author Information:
	
        A.  Name: Louis Bliard
		Institution: Aarhus University
		Email: bliard.louis@gmail.com
		
        B.  Name: Michael Griesser
		Institution: University of Konstanz
		
        C.  Name: Ryan Germain
		Institution: Aarhus University
		
4. Date of data collection: NA

5. Geographic location of data collection: Global


## SHARING/ACCESS INFORMATION

1. Licenses/restrictions placed on the data: CC-BY 4.0

2. Links to publications that cite or use the data: XXXXX

3. Links to other publicly accessible locations of the data: https://livingplanetindex.org/data_portal ; https://doi.org/10.24381/cds.f17050d7 ; https://doi.org/10.1073/pnas.2409658122 ; https://doi.org/10.1111/ele.13898

4. Was data derived from another source?
- The demographic time series data was obtained from the Living Planet Database https://livingplanetindex.org/data_portal 
- The climatic data was obtained from the Copernicus Climate Data Store https://doi.org/10.24381/cds.f17050d7
Information on the specifics of the climatic data obtained for this project can be found in this repository in the file `receipt_climate_data.txt` 
- The maximum clade credibility tree was obtained from McTavish et al. 2025 v.1.6. https://doi.org/10.1073/pnas.2409658122 & https://github.com/McTavishLab/AvesData 
- The functional trait data was obtained from Avonet; Tobias et al. 2022 https://doi.org/10.1111/ele.13898
- Some of the code for extracting relevant data from Living Planet Database and for some analyses was based on Jackson et al 2022 article https://doi.org/10.7554/eLife.74161 and repository https://doi.org/10.5281/zenodo.6620489

5. Recommended citation for this dataset: Bliard L, Griesser M, Germain R. (XXXX) Populations of cooperatively breeding birds are buffered against annual variation in rainfall and temperature. [Data set].
If anyone wishes to use any of the data originating from the primary sources cited above, please cite the primary literature and not this repository.

## DATA & FILE OVERVIEW

1. File List:

a. data


b. code
b1. main

b2. supplements


2. Relationship between files, if important: 



## METHODOLOGICAL INFORMATION
 
1. Methods for processing the data: R

2. Instrument- or software-specific information needed to interpret the data: 
See section "SESSION INFORMATION" below

3. People involved with sample collection: Michael Griesser (for the data on cooperative breeding / family living)

4. People involved with data formatting and comparative analyses: Louis Bliard

## SESSION INFORMATION

```
R version 4.5.3 (2026-03-11)
Platform: aarch64-apple-darwin20
Running under: macOS Tahoe 26.5.2

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: Europe/Copenhagen
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] phangorn_2.12.1      ape_5.8-1            ggdist_3.3.3         mgcv_1.9-4           nlme_3.1-168        
 [6] DHARMa_0.4.7         exactextractr_0.10.1 tictoc_1.2.1         rasterVis_0.51.7     lattice_0.22-9      
[11] terra_1.9-11         raster_3.6-32        sp_2.2-1             patchwork_1.3.2      lubridate_1.9.5     
[16] forcats_1.0.1        stringr_1.6.0        dplyr_1.2.0          purrr_1.2.1          readr_2.2.0         
[21] tidyr_1.3.2          tibble_3.3.1         tidyverse_2.0.0      data.table_1.18.2.1  brms_2.23.0         
[26] Rcpp_1.1.1           ggplot2_4.0.2       

loaded via a namespace (and not attached):
 [1] tidyselect_1.2.1      viridisLite_0.4.3     farver_2.1.2          loo_2.9.0             S7_0.2.1             
 [6] tensorA_0.36.2.1      digest_0.6.39         timechange_0.4.0      lifecycle_1.0.5       sf_1.1-0             
[11] magrittr_2.0.4        posterior_1.6.1       compiler_4.5.3        rlang_1.1.7           tools_4.5.3          
[16] igraph_2.2.2          bridgesampling_1.2-1  interp_1.1-6          classInt_0.4-11       RColorBrewer_1.1-3   
[21] abind_1.4-8           KernSmooth_2.23-26    withr_3.0.2           grid_4.5.3            latticeExtra_0.6-31  
[26] e1071_1.7-17          scales_1.4.0          MASS_7.3-65           cli_3.6.5             mvtnorm_1.3-7        
[31] reformulas_0.4.4      generics_0.1.4        otel_0.2.0            RcppParallel_5.1.11-2 rstudioapi_0.18.0    
[36] tzdb_0.5.0            minqa_1.2.8           DBI_1.3.0             proxy_0.4-29          splines_4.5.3        
[41] bayesplot_1.15.0      parallel_4.5.3        matrixStats_1.5.0     vctrs_0.7.2           boot_1.3-32          
[46] Matrix_1.7-4          hms_1.1.4             jpeg_0.1-11           hexbin_1.28.5         units_1.0-1          
[51] glue_1.8.0            nloptr_2.2.1          codetools_0.2-20      distributional_0.7.0  stringi_1.8.7        
[56] gtable_0.3.6          deldir_2.0-4          quadprog_1.5-8        lme4_2.0-1            pillar_1.11.1        
[61] Brobdingnag_1.2-9     R6_2.6.1              Rdpack_2.6.6          rbibutils_2.4.1       png_0.1-9            
[66] backports_1.5.0       rstantools_2.6.0      class_7.3-23          fastmatch_1.1-8       coda_0.19-4.1        
[71] checkmate_2.3.4       zoo_1.8-15            pkgconfig_2.0.3      
```

## DATA-SPECIFIC INFORMATION FOR: `XXXX.txt`

1. Number of variables: 

2. Number of cases/rows: 
   Each row correspond to one population

3. Variable List: 

4. Missing data codes: NA
