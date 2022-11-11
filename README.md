# Buonanduci_etal_2022_LandscEcol

Data accompanying the manuscript 'Fine-scale spatial heterogeneity shapes compensatory responses of a subalpine forest to severe bark beetle outbreak' by Buonanduci, Morris, Agne, Battaglia, and Harvey published in Landscape Ecology. See the main text of the manuscript for complete descriptions of data collection and processing.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7314248.svg)](https://doi.org/10.5281/zenodo.7314248)

[![CC BY 4.0][cc-by-shield]][cc-by]

This information is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by]. Any user of these data ("User" hereafter) is required to cite it appropriately in any publication that results from its use. These data may be actively used by others for ongoing research, so coordination may be necessary to prevent duplicate publication. The User is urged to contact the authors of these data for questions about methodology or results.  The User is encouraged to consider collaboration or co-authorship with authors where appropriate. Misinterpretation of data may occur if used out of context of the original study. Substantial efforts are made to ensure accuracy of the data and documentation, however complete accuracy of data sets cannot be guaranteed. All data are made available as is. Data may be updated periodically and it is the responsibility of the User to check for new versions of the data. The authors and the repository where these data were obtained shall not be liable for damages resulting from any use or misinterpretation of the data.

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg


There are two data files made available for reproducibility purposes:
- input_individual_trees.csv
- input_aggregate_raster.csv


### input_individual_trees.csv
This file contains individual tree-level data necessary to run the tree-scale analyses. The following columns are included:
- **Plot**: unique plot identifier.
- **TreeID**: unique tree identifier.
- **Tag**: tree tag number.
- **Species**: four letter code indicating the species of each measured tree. 
  - ABLA = *Abies lasiocarpa*, subalpine fir
  - PICO = *Pinus contorta*, lodgepole pine
  - PIEN = *Picea engelmannii*, Engelmann spruce
- **X**: easting of each measured tree in NAD_1983_UTM_Zone_13N.
- **Y**: northing of each measured tree in NAD_1983_UTM_Zone_13N.
- **DBH_89**: tree diameter at breast height (DBH; cm), measured in 1989. NA's indicate tree DBH was less than minimum size threshold in 1989.
- **DBH_04**: tree diameter at breast height (DBH; cm), measured in 2004. NA's indicate tree DBH was less than minimum size threshold in 2004.
- **DBH_18**: tree diameter at breast height (DBH; cm), measured in 2018. NA's indicate either (a) tree DBH was less than minimum size threshold in 2018 or (b) tree was recorded dead in 2018.
- **Dia_89_04**: annualized diameter increment (cm/year) between 1989 and 2004, calculated as (DBH_04 - DBH_89) / (2004 - 1989).
- **Dia_04_18**: annualized diameter increment (cm/year) between 2004 and 2018, calculated as (DBH_18 - DBH_04) / (2018 - 2004).
- **Dia_dif**: difference in annualized diameter increment between the two growth periods, calculated as Dia_04_18 - Dia_89_04.
- **NDI_89_live**: neighborhood density index (NDI) for all live trees in 1989. Calculated within a 10 m radius of focal tree using DBH recorded in 1989.
- **NDI_89_consp**: neighborhood density index (NDI) for all conspecific live trees in 1989. Calculated within a 10 m radius of focal tree using DBH recorded in 1989.
- **NDI_89_mort**: neighborhood density index (NDI) for all trees recorded as live in 1989 and dead in 2004. Calculated within a 10 m radius of focal tree using DBH recorded in 1989.
- **NDI_89_propconsp**: proportion of NDI contributed by the same species as the focal tree, calculated as NDI_89_consp / NDI_89_live.
- **NDI_89_propmort**: proportion of NDI that was killed during the growth period, calculated as NDI_89_mort / NDI_89_live.
- **NDI_04_live**: neighborhood density index (NDI) for all live trees in 2004. Calculated within a 10 m radius of focal tree using DBH recorded in 2004.
- **NDI_04_consp**: neighborhood density index (NDI) for all conspecific live trees in 2004. Calculated within a 10 m radius of focal tree using DBH recorded in 2004.
- **NDI_04_mort**: neighborhood density index (NDI) for all trees recorded as live in 2004 and dead in 2018. Calculated within a 10 m radius of focal tree using DBH recorded in 2004.
- **NDI_04_propconsp**: proportion of NDI contributed by the same species as the focal tree, calculated as NDI_04_consp / NDI_04_live.
- **NDI_04_propmort**: proportion of NDI that was killed during the growth period, calculated as NDI_04_mort / NDI_04_live.
- **Slope**: topographic slope (degrees).
- **TPI**: topographic position index (TPI; m) calculated within a 10 m radius of the focal tree.


### input_aggregate_raster.csv
This file contains 10 m raster cell data necessary to run the within-stand aggregate growth analyses. The following columns are included:
- **Plot**: unique plot identifier.
- **X**: easting of raster cell in NAD_1983_UTM_Zone_13N.
- **Y**: northing of each raster cell in NAD_1983_UTM_Zone_13N.
- **grow_BA_89_04**: annualized basal area (BA) growth of live trees between 1989 and 2004 (m<sup>2</sup>ha<sup>-1</sup>yr<sup>-1</sup>), calculated using only those trees that were alive and surveyed in both 1989 and 2004.
- **grow_BA_04_18**: annualized basal area (BA) growth of live trees between 2004 and 2018 (m<sup>2</sup>ha<sup>-1</sup>yr<sup>-1</sup>), calculated using only those trees that were alive and surveyed in both 2004 and 2018.
- **grow_BA_dif**: difference in BA growth between the two growth periods, calculated as grow_BA_04_18 - grow_BA_89_04.
- **ingrow_stems_89_04**: overstory recruitment between 1989 and 2004 (number of stems), calculated as the number of trees that grew to exceed the minimum survey diameter between 1989 and 2004
- **ingrow_stems_04_18**: overstory recruitment between 2004 and 2018 (number of stems), calculated as the number of trees that grew to exceed the minimum survey diameter between 2004 and 2018
- **ingrow_stems_dif**: : difference in overstory recruitment between the two growth periods, calculated as ingrow_stems_04_18 - ingrow_stems_89_04.
- **live_BA_89**: live basal area (BA) in 1989 (m<sup>2</sup>/ha).
- **live_BA_04**: live basal area (BA) in 2004 (m<sup>2</sup>/ha).
- **live_stems_89**: live stems in 1989 (number of stems).
- **live_stems_04**: live stems in 2004 (number of stems).
- **live_BA_LS_89**: live basal area (BA) of late-seral species in 1989 (m<sup>2</sup>/ha).
- **live_BA_LS_04**: live basal area (BA) of late-seral species in 2004 (m<sup>2</sup>/ha).
- **surv_BA_89**: basal area (BA) of trees recorded as live in both 1989 and 2004 (m<sup>2</sup>/ha). Calculated using DBH recorded in 1989.
- **surv_BA_04**: basal area (BA) of trees recorded as live in both 2004 and 2018 (m<sup>2</sup>/ha). Calculated using DBH recorded in 2004.
- **surv_stems_89**: number of trees recorded as live in both 1989 and 2004.
- **surv_stems_04**: number of trees recorded as live in both 2004 and 2018.
- **surv_BA_LS_89**: basal area (BA) of late-seral species trees recorded as live in both 1989 and 2004 (m<sup>2</sup>/ha). Calculated using DBH recorded in 1989.
- **surv_BA_LS_04**: basal area (BA) of late-seral species trees recorded as live in both 2004 and 2018 (m<sup>2</sup>/ha). Calculated using DBH recorded in 2004.
- **mort_BA_89_04**: basal area (BA) of trees recorded as live in 1989 and dead in 2004 (m<sup>2</sup>/ha). Calculated using DBH recorded in 1989.
- **mort_BA_04_18**: basal area (BA) of trees recorded as live in 2004 and dead in 2018 (m<sup>2</sup>/ha). Calculated using DBH recorded in 2004.
- **live_BA_LS_89_prop**: proportion of live basal area (BA) composed of late-seral species, calculated as live_BA_LS_89 / live_BA_89.
- **live_BA_LS_04_prop**: proportion of live basal area (BA) composed of late-seral species, calculated as live_BA_LS_04 / live_BA_04.
- **surv_BA_LS_89_prop**: proportion of surviving basal area (BA) composed of late-seral species, calculated as surv_BA_LS_89 / surv_BA_89.
- **surv_BA_LS_04_prop**: proportion of surviving basal area (BA) composed of late-seral species, calculated as surv_BA_LS_04 / surv_BA_04.
- **mort_BA_89_04_prop**: proportion of live basal area (BA) that died during the growth period, calculated as mort_BA_89_04 / live_BA_89.
- **mort_BA_04_18_prop**: proportion of live basal area (BA) that died during the growth period, calculated as mort_BA_04_18 / live_BA_04.
- **Slope**: topographic slope (degrees).
- **TPI**: topographic position index (TPI; m) calculated within a 10 m radius of the raster cell center coordinates.

