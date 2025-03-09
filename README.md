# Analysing & comparing disturbances in set aside forest areas for conservation with mananged fores.
This repository holds code for running the analysis of comparing disturbance regimes in set aside forests versus comparable managed forests across Germany. For a detailed description of the purpose, methodology and results, see the following paper: 

Krüger, K., Senf, C., Hagge, J., Seidl, R. (2025). Setting aside areas for conservation does not increase disturbances in temperate forests. Journal of Applied Ecology. 

The data to reproduce the analysis can be obtained under the following link: (link)

All data to reproduce the results are available, all other layers can be generated with the code provided in this repository. The vector files delineating set aside forest areas can just partially be made available due to data usage agreements with the respective regional authorities - for reproducibility, please request the respective layers at the regional authorities, while we provide links to those you can openly access. 

## platforms

software used for development and implementation : R version 4.3.3

## Scrips:
Scripts are named in order (01_, etc.). 

- `01_preprocess_environmental_data`: Process all environmental data layers used in the analysis (e.g. align resolutions, extract specific sub-information or change from vector to raster format).
- `02_agregate_tree_species`: Script to process the national tree species map from Blickensdörfer et al. (2024). We aggregate to more coarse forest types, which serve as input variable in the matching of comparable landscapes.
- `03_prepare_matching_rasters`: Aggregating the environmental information layers with focal windows to classify the landscape characteristics of this focal window in the centre cell. Those aggregated layers are converted into a data frame and serve as input for the matching of comparable landscapes. Additionally, in this script we also randomly draw set aside landscapes, for which we search comparable managed landscapes.
- `04_matching`: Script to match managed landscapes to pre-selected set aside sublandscapes.
- `05_get_disturbance_patchID_for_extraction`: Extract disturbance patches and assign unique patch IDs, which serve as input for the forest disturbance information extraction to calculate the maximum patch size.
- `06_extract_disturbance_metrics`: Extracting the geospecific forest disturbance information for each sublandscape. 
- `07_define_pulse_dist_years`: Script to identify years with very high occurence of forest distrubances, which we label as "pulse disturbance years".
- `08_modelling`: Model fitting for all forest distrubance metrics to identify the effect on management and setting aside forests on the distrubance rate, maximum patch size, patch density, disturance frequency and high severity rate of disturbances.
- `9_detection_accuracy`: Validation of probabilies to detetct forest disturbances initiated by natural disturbance agents in managed and set aside forest areas.
- `10_npv_comparison`: We benchmark the composition or manifested forest types in set aside forest areas against the Natural Pontential Vegetation (after Bohn et al. 2004).
