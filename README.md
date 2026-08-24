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
─ Session info ──────────────────────────────────────────────────────────────────────────────────────────────────
 setting  value
 version  R version 4.5.3 (2026-03-11)
 os       macOS Tahoe 26.5.2
 system   aarch64, darwin20
 ui       RStudio
 language (EN)
 collate  en_US.UTF-8
 ctype    en_US.UTF-8
 tz       Europe/Copenhagen
 date     2026-07-29
 rstudio  2026.01.1+403 Apple Blossom (desktop)
 pandoc   NA
 quarto   1.8.25 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto

─ Packages ──────────────────────────────────────────────────────────────────────────────────────────────────────
 package        * version  date (UTC) lib source
 abind            1.4-8    2024-09-12 [1] CRAN (R 4.5.0)
 ape            * 5.8-1    2024-12-16 [1] CRAN (R 4.5.0)
 backports        1.5.0    2024-05-23 [1] CRAN (R 4.5.0)
 bayesplot        1.15.0   2025-12-12 [1] CRAN (R 4.5.2)
 boot             1.3-32   2025-08-29 [1] CRAN (R 4.5.3)
 bridgesampling   1.2-1    2025-11-19 [1] CRAN (R 4.5.2)
 brms           * 2.23.0   2025-09-09 [1] CRAN (R 4.5.0)
 Brobdingnag      1.2-9    2022-10-19 [1] CRAN (R 4.5.0)
 checkmate        2.3.4    2026-02-03 [1] CRAN (R 4.5.2)
 class            7.3-23   2025-01-01 [1] CRAN (R 4.5.3)
 classInt         0.4-11   2025-01-08 [1] CRAN (R 4.5.0)
 cli              3.6.5    2025-04-23 [1] CRAN (R 4.5.0)
 coda             0.19-4.1 2024-01-31 [1] CRAN (R 4.5.0)
 codetools        0.2-20   2024-03-31 [1] CRAN (R 4.5.3)
 data.table     * 1.18.2.1 2026-01-27 [1] CRAN (R 4.5.2)
 DBI              1.3.0    2026-02-25 [1] CRAN (R 4.5.2)
 deldir           2.0-4    2024-02-28 [1] CRAN (R 4.5.0)
 DHARMa         * 0.4.7    2024-10-18 [1] CRAN (R 4.5.0)
 digest           0.6.39   2025-11-19 [1] CRAN (R 4.5.2)
 distributional   0.7.0    2026-03-17 [1] CRAN (R 4.5.2)
 dplyr          * 1.2.0    2026-02-03 [1] CRAN (R 4.5.2)
 e1071            1.7-17   2025-12-18 [1] CRAN (R 4.5.2)
 exactextractr  * 0.10.1   2025-12-01 [1] CRAN (R 4.5.2)
 farver           2.1.2    2024-05-13 [1] CRAN (R 4.5.0)
 fastmatch        1.1-8    2026-01-17 [1] CRAN (R 4.5.2)
 forcats        * 1.0.1    2025-09-25 [1] CRAN (R 4.5.0)
 generics         0.1.4    2025-05-09 [1] CRAN (R 4.5.0)
 ggdist         * 3.3.3    2025-04-23 [1] CRAN (R 4.5.0)
 ggplot2        * 4.0.2    2026-02-03 [1] CRAN (R 4.5.2)
 glue             1.8.0    2024-09-30 [1] CRAN (R 4.5.0)
 gtable           0.3.6    2024-10-25 [1] CRAN (R 4.5.0)
 hexbin           1.28.5   2024-11-13 [1] CRAN (R 4.5.0)
 hms              1.1.4    2025-10-17 [1] CRAN (R 4.5.0)
 igraph           2.2.2    2026-02-12 [1] CRAN (R 4.5.2)
 interp           1.1-6    2024-01-26 [1] CRAN (R 4.5.0)
 jpeg             0.1-11   2025-03-21 [1] CRAN (R 4.5.0)
 KernSmooth       2.23-26  2025-01-01 [1] CRAN (R 4.5.3)
 lattice        * 0.22-9   2026-02-09 [1] CRAN (R 4.5.3)
 latticeExtra     0.6-31   2025-09-10 [1] CRAN (R 4.5.0)
 lifecycle        1.0.5    2026-01-08 [1] CRAN (R 4.5.2)
 lme4             2.0-1    2026-03-05 [1] CRAN (R 4.5.2)
 loo              2.9.0    2025-12-23 [1] CRAN (R 4.5.2)
 lubridate      * 1.9.5    2026-02-04 [1] CRAN (R 4.5.2)
 magrittr         2.0.4    2025-09-12 [1] CRAN (R 4.5.0)
 MASS             7.3-65   2025-02-28 [1] CRAN (R 4.5.3)
 Matrix           1.7-4    2025-08-28 [1] CRAN (R 4.5.3)
 matrixStats      1.5.0    2025-01-07 [1] CRAN (R 4.5.0)
 mgcv           * 1.9-4    2025-11-07 [1] CRAN (R 4.5.3)
 minqa            1.2.8    2024-08-17 [1] CRAN (R 4.5.0)
 mvtnorm          1.3-7    2026-04-15 [1] CRAN (R 4.5.2)
 nlme           * 3.1-168  2025-03-31 [1] CRAN (R 4.5.3)
 nloptr           2.2.1    2025-03-17 [1] CRAN (R 4.5.0)
 otel             0.2.0    2025-08-29 [1] CRAN (R 4.5.0)
 patchwork      * 1.3.2    2025-08-25 [1] CRAN (R 4.5.0)
 phangorn       * 2.12.1   2024-09-17 [1] CRAN (R 4.5.0)
 pillar           1.11.1   2025-09-17 [1] CRAN (R 4.5.0)
 pkgconfig        2.0.3    2019-09-22 [1] CRAN (R 4.5.0)
 png              0.1-9    2026-03-15 [1] CRAN (R 4.5.2)
 posterior        1.6.1    2025-02-27 [1] CRAN (R 4.5.0)
 proxy            0.4-29   2025-12-29 [1] CRAN (R 4.5.2)
 purrr          * 1.2.1    2026-01-09 [1] CRAN (R 4.5.2)
 quadprog         1.5-8    2019-11-20 [1] CRAN (R 4.5.0)
 R6               2.6.1    2025-02-15 [1] CRAN (R 4.5.0)
 raster         * 3.6-32   2025-03-28 [1] CRAN (R 4.5.0)
 rasterVis      * 0.51.7   2025-09-01 [1] CRAN (R 4.5.0)
 rbibutils        2.4.1    2026-01-21 [1] CRAN (R 4.5.2)
 RColorBrewer     1.1-3    2022-04-03 [1] CRAN (R 4.5.0)
 Rcpp           * 1.1.1    2026-01-10 [1] CRAN (R 4.5.2)
 RcppParallel     5.1.11-2 2026-03-05 [1] CRAN (R 4.5.2)
 Rdpack           2.6.6    2026-02-08 [1] CRAN (R 4.5.2)
 readr          * 2.2.0    2026-02-19 [1] CRAN (R 4.5.2)
 reformulas       0.4.4    2026-02-02 [1] CRAN (R 4.5.2)
 rlang            1.1.7    2026-01-09 [1] CRAN (R 4.5.2)
 rstantools       2.6.0    2026-01-10 [1] CRAN (R 4.5.2)
 rstudioapi       0.18.0   2026-01-16 [1] CRAN (R 4.5.2)
 S7               0.2.1    2025-11-14 [1] CRAN (R 4.5.2)
 scales           1.4.0    2025-04-24 [1] CRAN (R 4.5.0)
 sessioninfo      1.2.4    2026-06-04 [1] CRAN (R 4.5.2)
 sf               1.1-0    2026-02-24 [1] CRAN (R 4.5.2)
 sp             * 2.2-1    2026-02-13 [1] CRAN (R 4.5.2)
 stringi          1.8.7    2025-03-27 [1] CRAN (R 4.5.0)
 stringr        * 1.6.0    2025-11-04 [1] CRAN (R 4.5.0)
 tensorA          0.36.2.1 2023-12-13 [1] CRAN (R 4.5.0)
 terra          * 1.9-11   2026-03-26 [1] CRAN (R 4.5.2)
 tibble         * 3.3.1    2026-01-11 [1] CRAN (R 4.5.2)
 tictoc         * 1.2.1    2024-03-18 [1] CRAN (R 4.5.0)
 tidyr          * 1.3.2    2025-12-19 [1] CRAN (R 4.5.2)
 tidyselect       1.2.1    2024-03-11 [1] CRAN (R 4.5.0)
 tidyverse      * 2.0.0    2023-02-22 [1] CRAN (R 4.5.0)
 timechange       0.4.0    2026-01-29 [1] CRAN (R 4.5.2)
 tzdb             0.5.0    2025-03-15 [1] CRAN (R 4.5.0)
 units            1.0-1    2026-03-11 [1] CRAN (R 4.5.2)
 vctrs            0.7.2    2026-03-21 [1] CRAN (R 4.5.2)
 viridisLite      0.4.3    2026-02-04 [1] CRAN (R 4.5.2)
 withr            3.0.2    2024-10-28 [1] CRAN (R 4.5.0)
 zoo              1.8-15   2025-12-15 [1] CRAN (R 4.5.2)

 [1] /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library
 * ── Packages attached to the search path.
```

## DATA-SPECIFIC INFORMATION FOR: `XXXX.txt`

1. Number of variables: 

2. Number of cases/rows: 
   Each row correspond to one population

3. Variable List: 

4. Missing data codes: NA
