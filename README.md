# Broad-scale benthic habitat classification of the South Atlantic

The code and files contained in this repository support replication of a broad-scale benthic habitat classification of the South Atlantic produced by McQuaid et al. (2023). We used statistical clustering algorithms to classify broad-scale (10km2) environmental data into distinct habitat classes, which reflect variation in physical conditions and we assume support distinct biological communities. 

We request that any use of the input data is referenced as per the table below, and that classification outputs are referenced as:
McQuaid K. A. Bridges A. E. H., Howell K. L., Gandra T. B. R., de Souza V., Currie J. C., Hogg O. T., Pearman T. R. R., Bell J. B. B., Atkinson L. J., Baum D., Bonetti J., Carranza A., Defeo O., Furey T., Gasalla M. A., Golding N, Hampton S. L., Horta S., Jones D. O. B., Lombard A. T., Manca E., Marin Y., Martin S., Mortensen P., Passdore C., Piechaud N., Sink K. J. & Yool A. 2023. Broad-scale benthic habitat classification of the South Atlantic. Progress in Oceanography. DOI: https://doi.org/10.1016/j.pocean.2023.103016.

<br />__Process for running the scripts:__
* Step 1: Run the non-hierarchical approach using the script "NHApproach NP Clustering & Final Classification.R" and input data layers detailed in the table below
* Step 2: Run Level 1 of the hierarchical approach using the script "HApproach Level1 - PCAs.R" and input data layers detailed in the table below
* Step 3: Produce a confidence layer for Level 1 of the hierarchical approach using the script "HApproach Level1 - confusion Index.R"
* Step 4: Run Level 2 of the hierarchical approach using the script "HApproach_Level2,3_RScript.R" and the classification outputs of Step 3 (Level 1) as inputs
* Step 5: Adapt the script "HApproach_Level2,3_RScript.R" to run Level 3 of the hierarchical approach using the classification outputs of Step 4 (Level 2) as inputs


<br /> __Folder/file descriptions:__
| Folder/file name     | Description |    
| ------------- | ------------- | 
| scripts     |  Folder containing scripts to run the analyses |
| inputs     |  Folder containing mask for South Atlantic used in this study and biogeography layer produced and used in this study |
| outputs     |  Folder containing a link to download the final habitat classifications produced in this study |
| NHApproach NP Clustering & Final Classification.R    |  R script for data preparation, clustering, final classification and production of confidence layers for the non-hierarchical approach |
| HApproach Level1 - PCAs NP.R   |  R script for data preparation, Principal Component Analysis and clustering for Level 1 of the hierarchical approach |
| HApproach Level1 - confusion Index.R   |  R script for production of confidence layers using a confusion index for Level 1 outputs of the hierarchical approach |
| HApproach_Level2,3_RScript NP.R  |  R script for data preparation, Principal Component Analysis and clustering for Level 2 of the hierarchical approach |
| Mask_SouthAtlantic.shp | Mask of the South Atlantic used in this study. The boundaries come from the International Hydrographic Organization’s definition of the South Atlantic (Flanders Marine Institute, 2020), but are extended north, south, and east to include the entire EEZs of Brazil, South Georgia and the South Sandwich Islands, and South Africa, respectively. |
| BiogeographicRegions.tif | Developed using the outputs of water mass structure analysis (see non-hierarchical approach in paper) and published biogeographic classifications where possible (e.g. Vinogradova, 1997; Watling et al., 2013; Zezina, 1997) | 


<br /> __Input raster datasets are available at the following sites:__
| Variable     | Description      | Manipulation     |  Source |
| ------------- | ------------- | -------- |  -------- |
| Depth  | Continuous bathymetric model | Resampled to 5 arc min resolution |  https://www.gebco.net/data_and_products/gridded_bathymetry_data/gebco_2020/ |
| BBPI  | Measure of where a referenced location is relative to the locations surrounding it |Created in ArcGIS from depth data using Benthic Terrain Modeler extension. Inner radius 1, outer radius 10, scale factor is ∼100km |  https://www.gebco.net/data_and_products/gridded_bathymetry_data/gebco_2020/ |
| FBPI  | Measure of where a referenced location is relative to the locations surrounding it | Created in ArcGIS from depth data using Benthic Terrain Modeler extension. Inner radius 1, outer radius 2, scale factor is ∼10km | https://www.gebco.net/data_and_products/gridded_bathymetry_data/gebco_2020/  |
| Slope  | Gradient, or rate of maximum change in z-value | Created in ArcGIS from depth data using Benthic Terrain Modeler extension | https://www.gebco.net/data_and_products/gridded_bathymetry_data/gebco_2020/   |
| Salinity | Mean benthic salinity at mean depth for the period 2000-2014 | NA | https://www.bio-oracle.org/ |
| Temperature | Mean benthic temperature at mean depth for the period 2000-2014 | NA | https://www.bio-oracle.org/  |
| Productivity| Mean particulate organic carbon ﬂux to seafloor for the period 2006-2015 | Output from the MEDUSA model (Yool et al. 2013) regridded from ORCA0083 to NEMO 5 arc min | https://zenodo.org/record/6513616#.Yn6TGx1BzIU |
| Dissolved oxygen | Mean benthic dissolved oxygen concentration at mean depth for the period 2000-2014 | NA | https://www.bio-oracle.org/ |
| Nitrate | Mean benthic nitrate concentration at mean depth for the period 2000-2014 | NA | https://www.bio-oracle.org/  |
| Phosphate | Mean benthic phosphate concentration at mean depth for the period 2000-2014| NA | https://www.bio-oracle.org/ |
| Silicate | Mean benthic silicate concentration at mean depth for the period 2000-2014 | NA | https://www.bio-oracle.org/  |
| Current velocity |  Mean benthic current velocity at mean depth for the period 2000-2014 | NA | https://www.bio-oracle.org/  |
| Biogeography | Biogeographic regions | Developed using the outputs of water mass structure analysis (see non-hierarchical approach in paper) and published biogeographic classifications where possible (e.g. Vinogradova, 1997; Watling et al., 2013; Zezina, 1997) |  See folder of input data layers  |
