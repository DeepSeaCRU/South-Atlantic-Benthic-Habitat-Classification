# Regional benthic habitat classification of the South Atlantic

The code and files contained in this repository support replication of a regional benthic habitat classification of the South Atlantic produced by McQuaid et al. (in prep). We used statistical clustering algorithms to classify broad-scale (10km2) environmental data into distinct habitat classes, which reflect variation in physical conditions and we assume support distinct biological communities. 

******* LINE ON CONDITIONS OF USE? *******

Process for running the script:
<br /> Step 1: Run [insert script name]
<br /> Step 2: Run [insert script name]
<br /> Step 3: Run [insert script name], using classification outputs of Step 2 as inputs
<br /> Step 3: Run [insert script name], using classification outputs of Step 3 as inputs

Folder/file descriptions:
| File/folder name     | Description |    
| ------------- | ------------- | 
| scripts     |  Contains scripts to run the analyses |
| layers     |  Contains mask for South Atlantic used in this study, biogeography layer produced for this study, and output habitat classifications for the hierarchical and non-hierarchical approaches |
| script name 1    |  script function |
| script name 2    |  script function |
| script name 3    |  script function |
| script name 4   |  script function |
... etc

Input raster datasets are available at the following sites:
| Variable     | Description      | Manipulation     |  Source |
| ------------- | ------------- | -------- |  -------- |
| Depth  | Continuous bathymetric model | Resampled to 5 arc min resolution |  https://www.gebco.net/data_and_products/gridded_bathymetry_data/gebco_2020/ |
| BBPI  | Measure of where a referenced location is relative to the locations surrounding it |Created in ArcGIS from depth data using Benthic Terrain Modeler extension. Inner radius 1, outer radius 10, scale factor is ∼100km |  https://www.gebco.net/data_and_products/gridded_bathymetry_data/gebco_2020/ |
| FBPI  | Measure of where a referenced location is relative to the locations surrounding it | Created in ArcGIS from depth data using Benthic Terrain Modeler extension. Inner radius 1, outer radius 2, scale factor is ∼10km | https://www.gebco.net/data_and_products/gridded_bathymetry_data/gebco_2020/  |
| Slope  | Gradient, or rate of maximum change in z-value | Created in ArcGIS from depth data using Benthic Terrain Modeler extension | https://www.gebco.net/data_and_products/gridded_bathymetry_data/gebco_2020/   |
| Salinity | Mean benthic salinity at mean depth for the period 2000-2014 | NA | https://www.bio-oracle.org/ |
| Temperature | Mean benthic temperature at mean depth for the period 2000-2014 | NA | https://www.bio-oracle.org/  |
| Productivity| Mean particulate organic carbon ﬂux to seafloor for the period 2006-2015 | Output from the MEDUSA model (Yool et al. 2013) regridded from ORCA0083 to NEMO 5 arc min | TBC |
| Dissolved oxygen | Mean benthic dissolved oxygen concentration at mean depth for the period 2000-2014 | NA | https://www.bio-oracle.org/ |
| Nitrate | Mean benthic nitrate concentration at mean depth for the period 2000-2014 | NA | https://www.bio-oracle.org/  |
| Phosphate | Mean benthic phosphate concentration at mean depth for the period 2000-2014| NA | https://www.bio-oracle.org/ |
| Silicate | Mean benthic silicate concentration at mean depth for the period 2000-2014 | NA | https://www.bio-oracle.org/  |
| Current velocity |  Mean benthic current velocity at mean depth for the period 2000-2014 | NA | https://www.bio-oracle.org/  |
| Biogeography | Biogeographic provinces | Developed using the outputs of water mass structure analysis (see non-hierarchical approach in paper) and published biogeographic classifications where possible (e.g. Vinogradova, 1997; Watling et al., 2013; Zezina, 1997) |  https://www.bio-oracle.org/  |


