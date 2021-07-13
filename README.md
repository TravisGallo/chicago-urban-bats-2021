
# chicago-urban-bats-2021


A repository that contains the data and code for:

Lehrer, E.W, T. Gallo, M. Fidino, R.J. Kilgour, P.J. Wolff, S. Magle. **2021**. Urban bat occupancy is highly influenced by noise and the location of water: Considerations for nature-based urban planning. Landscape and Urban Planning. 210:104063.

This `README` file includes information on the various scripts and datasets used for this analysis. Not every data source is saved in this repository (e.g., LiDAR and GIS processing). However, the manuscript includes methodology and citation on where spatial datasets came from or how they were created.

**Note:** All of these files must be within your working directory for the analysis to work. Our analysis was done in parallel and used all but two of the cores in a computer. Therefore, if you have two or less cores on your computer you will need to adjust the settings.

---

<div align="center"> <h3>data</h3> </div>

---

**2020-06-26_Lehrer_etal_yarray.rds:** the number of days each species was detected. This data is organized in a three-dimensionl array with the first array being sites, the second dimension being species, and the third dimension being the sampling season.

**2020-06-26_Lehrer_etal_jmatrix.rds:** the number of days sampled at each site. These data are organized as a matrix with the sampling sites as rows and the sampling season as columns. The values are the number of days a bat detector was recording at each site and season. If a site was not sampled a zero is reported.

**2020-06-26_Lehrer_etal_HabitatCoV_2017.csv:** The proportion of impervious cover `imperv_1000`, canopy cover `canopy_1000`, road density `roads_1000`, edge density `edge_1000` within a 1-km fixed radius buffer around each site. Distance to nearest water source from each site `dist2water` and whether we considered a site urban or rural `trt`. The `trt` data was not used in this analysis.

**2020-06-26_Lehrer_etal_DetCovs2013_2017.csv:** Average monthly rainfall `AvgPrec...` for each sampling period and the model number of the recording device `Det...` that was deployed at each site during each sampling season. These data were used to estimate the detection probability of each bat species.

**2020-06-26_Lehrer_etal_building_xyz_allsites.rds:** The height and X,Y,Z coordinate of each building located within 1-km fixed radius buffer around each site. This data was generated from a building shapefile maintained by the City of Chicago and also generated from LiDAR data provided by Cook County (neither located in this repository).

**2020-06-26_Lehrer_etal_median_soundlevelpressure.rds:** Median sound pressure level at each site. These data were calculated using .wav files that are not located in this repository.

---

<div align="center"> <h3>scripts</h3> </div>

---

**2020-06-26_Lehrer_etal_bats2018_utility_funcs.R:** script that loads utility functions. Sourced in `2020-06-26_Lehrer_etal_urbanbatanalysis.R`.

**2020-06-26_Lehrer_etal_urbanbatanalysis.R:** the only file that needs to be run. Sources all files listed to load data sets, format data for JAGS model, run JAGS model, summarize posterior distributions, calculate derived parameters.

**2020-06-26_Lehrer_etal_JAGS_occu_LASSO.R:** JAGS model used to estimate the probability of occupancy for each bat species. Read in `2020-06-26_Lehrer_etal_urbanbatanalysis.R`

