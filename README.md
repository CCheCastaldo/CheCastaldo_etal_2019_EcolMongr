
### Disentangling herbivore impacts in primary succession by refocusing the plant stress and vigor hypotheses on phenology

Che-Castaldo, C., C. M. Crisafulli J. G. Bishop, E. F. Zipkin, and W. F. Fagan. 2019. *Ecological Monographs*.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3265488.svg)](https://doi.org/10.5281/zenodo.3265488)

Please contact the first author for questions about the code or data: Christian Che-Castaldo (<chrischecastaldo@icloud.com>)

------------------------------------------------------------------------

### Abstract

The plant stress and plant vigor hypotheses are widely used to explain the distribution and abundance of insect herbivores across their host plants. These hypotheses are the subject of contentious debate within the plant herbivore research community, with several studies finding simultaneous support for both hypotheses for the same plant-herbivore interaction. We address the question of how such support is possible using dynamic site-occupancy models to quantify the attack dynamics of *Cryptorynchus lapathi* (poplar-willow weevil) on *Salix sitchensis* (Sitka willow), a dioecious shrub colonizing Mount St. Helens after the 1980 eruption, in relation to host plant stress, vigor, and sex. We also introduce several scaling criteria as a rigorous test of the plant vigor hypothesis and demonstrate why modeling insect detection is important in plant-insect studies. Weevils responded positively to water stress associated with seasonal dry-downs, and this response was phenologically compartmentalized by larval feeding mode. Weevils preferentially attacked large and/or flowering stems, imposing an ecological cost on willow reproduction via increased stem mortality and susceptibility to future attack. We propose that the dual response to host plant stress and vigor is due to the synchronization between young weevil larval feeding and willow nutrient pulses that are mediated by environmental stress. In turn, this process drives successional dynamics, causing the juvenilization of upland willow plants and possibly delaying establishment of a willow-dominated upland sere. These results highlight the common, but often overlooked, phenological basis of the plant stress and plant vigor hypotheses, which both focus on how stress changes the quality of plant resources available to immature insects.

------------------------------------------------------------------------

### Repository Directory

1.  `Library` contains the 12 datasets used in this paper. The metadata for each dataset are listed below.

2.  `ModelBuild` contains the R and bash scripts used to fit the willow-weevil occupancy model in JAGS on AWS, as well as all model output.

3.  `FiguresTablesEtc` contains the R scripts used to make all figures and tables.

4.  `StemBiomass` contains the R script for the model of willow stem biomass as a function of willow basal stem diameter used to estimate the scaling coefficient *n*<sub>1</sub>.

5.  `StemMeasurementError` contains the R script for the model that estimates the routine error associated with measuring willow basal stem diameters using calipers.

------------------------------------------------------------------------

### Library Metadata

#### SalixBorderWide

Metadata for the willow-weevil dataset described in the Methods and Appendices S4-S7.

| field     | description                                                                                        |
|:----------|:---------------------------------------------------------------------------------------------------|
| transect  | name of transect, see methods and Appendix S4                                                      |
| point     | meter mark along transect, see Methods and Appendix S4                                             |
| plantid   | willow plant identifier used during field surveys, see Methods and Appendix S4                     |
| stemid    | willow 1<sup>st</sup> order stem identifier used during field surveys, see Methods and Appendix S5 |
| habitat   | habitat of willow plant: 0 = riparian, 1 = upland                                                  |
| sex       | sex of willow plant: 0 = male, 1 = female                                                          |
| plant     | willow plant identifier used in occupancy model                                                    |
| stem      | willow stem identifier used in occupancy model                                                     |
| site09    | status of willow stem in 2009: 0 = dead or did not exist, 1 = alive, see Methods and Appendix S5   |
| site10    | status of willow stem in 2010: 0 = dead or did not exist, 1 = alive, see Methods and Appendix S5   |
| site11    | status of willow stem in 2011: 0 = dead or did not exist, 1 = alive, see Methods and Appendix S5   |
| br09      | number of 2<sup>nd</sup> order branches on willow stem in 2009, see Methods and Appendix S5        |
| br10      | number of 2<sup>nd</sup> order branches on willow stem in 2010, see Methods and Appendix S5        |
| br11      | number of 2<sup>nd</sup> order branches on willow stem in 2011, see Methods and Appendix S5        |
| tbd09     | estimated total basal stem diameter in mm (tbd) in 2009, see Methods and Appendices S6-S7          |
| tbd10     | estimated total basal stem diameter in mm (tbd) in 2010, see Methods and Appendices S6-S7          |
| tbd11     | estimated total basal stem diameter in mm (tbd) in 2011, see Methods and Appendices S6-S7          |
| gr09      | relative growth rate (rgr) of tbd from 2009 to 2010, see Methods and Appendix S7                   |
| gr10      | relative growth rate (rgr) of tbd from 2010 to 2011, see Methods and Appendix S7                   |
| tbd09.sd  | standardized 2009 tbd, see Methods                                                                 |
| tbd10.sd  | standardized 2010 tbd, see Methods                                                                 |
| tbd11.sd  | standardized 2011 tbd, see Methods                                                                 |
| gr09.sd   | standardized 2009 to 2010 rgr, see Methods                                                         |
| gr10.sd   | standardized 2010 to 2011 rgr, see Methods                                                         |
| repro09   | stem reproductive status in 2009: 0 = non-reproductive, 1 = reproductive, see Methods              |
| repro10   | stem reproductive status in 2010: 0 = non-reproductive, 1 = reproductive, see Methods              |
| repro11   | stem reproductive status in 2011: 0 = non-reproductive, 1 = reproductive, see Methods              |
| borer1.09 | weevil larvae visual survey by observer 1 in 2009: 0 = no, 1 = yes, see Methods and Fig. S1-S2B    |
| borer2.09 | weevil larvae visual survey by observer 2 in 2009: 0 = no, 1 = yes, see Methods and Fig. S1-S2B    |
| borer1.10 | weevil larvae visual survey by observer 1 in 2010: 0 = no, 1 = yes, see Methods and Fig. S1-S2B    |
| borer2.10 | weevil larvae visual survey by observer 2 in 2010: 0 = no, 1 = yes, see Methods and Fig. S1-S2B    |
| borer1.11 | weevil larvae visual survey by observer 1 in 2011: 0 = no, 1 = yes, see Methods and Fig. S1-S2B    |
| borer2.11 | weevil larvae visual survey by observer 2 in 2011: 0 = no, 1 = yes, see Methods and Fig. S1-S2B    |
| tbd1.mean | mean tbd of upland (habitat = 1) 2009 stems, see Methods for explanation                           |
| tbd1.sd   | standard deviation of tbd of upland (habitat = 1) 2009 stems, see Methods for explanation          |
| tbd2.mean | mean tbd of 2009 riparian (habitat = 0) and all stems in 2010 and 2011                             |
| tbd2.sd   | standard deviation of 2009 riparian (habitat = 0) and all stems in 2010 and 2011                   |
| rgr.mean  | mean rgr of all stems                                                                              |
| rgr.sd    | standard deviation of rgr of all stems                                                             |

<br>

#### stemBiomass

Metadata for the harvested upland willow stem dataset described in the Methods and Appendix S8.

| field | description                  |
|:------|:-----------------------------|
| stem  | stem identifier              |
| soma  | aboveground stem biomass (g) |
| bd    | basal stem diameter (mm)     |

<br>

#### pairedStemMeasurements

Metadata for the caliper measurement error dataset described in the Methods and Appendix S6.

| field         | description                                     |
|:--------------|:------------------------------------------------|
| bd\_observer1 | basal stem diameter (mm) measured by observer 1 |
| bd\_observer2 | basal stem diameter (mm) measured by observer 2 |

<br>

#### transectPoints

Metadata for the `sf` dataset of upland and riparian transect point locations referenced in the willow-weevil dataset.

| field    | description                                            |
|:---------|:-------------------------------------------------------|
| transect | name of transect, see methods and Appendix S4          |
| point    | meter mark along transect, see Methods and Appendix S4 |
| type     | type of transect point, riparian vs. upland            |
| geometry | transect point coordinates in WGS84 (epsg:4326)        |

<br>

#### weevilPhenology

Metadata for the willow-phenology dataset described in the Methods and Appendix S9.

| field    | description                                                            |
|:---------|:-----------------------------------------------------------------------|
| stem\_id | willow stem identifier                                                 |
| date     | date willow stem harvested                                             |
| month    | month willow stem harvested                                            |
| year     | year willow stem harvested                                             |
| habitat  | habitat willow stem harvested from: riparian vs. upland                |
| stage    | life stage of weevils found inside willow stem: eggs vs. early instars |
| count    | number of weevil eggs or early instars per willow stem\_id             |

<br>

#### soilReleaseCurve

Metadata for the soil release curve dataset described in Appendices S2 and S3.

| field | description                                                      |
|:------|:-----------------------------------------------------------------|
| ec5   | raw sensor voltage (mV)                                          |
| vwc   | variable water content (m<sup>3</sup> water m<sup>-3</sup> soil) |
| kpa   | soil water potential (-kPa)                                      |

<br>

#### soilProbes

Metadata for the soil moisture dataset described and/or plotted in Appendices S2 and S3.

| field     | description                                                      |
|:----------|:-----------------------------------------------------------------|
| unit      | Decagon logger id                                                |
| port      | Decagon EC-5 soil moisture sensor port id                        |
| treatment | sensor emplaced in dry vs. irrigated soil                        |
| year      | year of sensor reading                                           |
| date      | date and time of sensor reading                                  |
| raw       | raw sensor voltage (mV)                                          |
| vwc       | variable water content (m<sup>3</sup> water m<sup>-3</sup> soil) |
| pfc       | percent soil field capacity                                      |

<br>

#### dailyPrecipitation

Metadata for the [Spirit Lake SNOTEL 777](https://wcc.sc.egov.usda.gov/nwcc/site?sitenum=777) daily precipitation dataset plotted in Appendix S2.

| field  | description                                          |
|:-------|:-----------------------------------------------------|
| date   | date and time of measurement                         |
| precip | daily precipitation (cm) based on the PREC.I-1 field |

<br>

#### monthlySnotel

Metadata for the [Spirit Lake SNOTEL 777](https://wcc.sc.egov.usda.gov/nwcc/site?sitenum=777) monthly precipitation dataset plotted in Appendix S2.

| field  | description                                                              |
|:-------|:-------------------------------------------------------------------------|
| year   | year of measurement                                                      |
| month  | month of measurement                                                     |
| precip | monthly precipitation (cm) based on the precipitation accumulation field |

<br>

#### hoboPumicePlain

Metadata for the hobo temperature and relative humidity dataset plotted in Appendix S2.

| field   | description                                                     |
|:--------|:----------------------------------------------------------------|
| month   | month                                                           |
| habitat | habitat Hobo loggers were deployed in: 0 = riparian, 1 = upland |
| temp.m  | mean daily temperature (°C)                                     |
| temp.sd | standard deviation of daily temperature (°C)                    |
| RH.m    | mean minimum daily relative humidity (percent)                  |
| RH.sd   | standard deviation minimum daily relative humidity (percent)    |

<br>

#### leafGasExchange

Metadata for the leaf gas exchange dataset described in Appendix S3.

| field          | description                                                                                                   |
|:---------------|:--------------------------------------------------------------------------------------------------------------|
| plant\_id      | plant tag identifier                                                                                          |
| rill           | artificial stream identifier, W vs. E, see Appendix S3                                                        |
| block          | identifier for 3 m subunits of stream containing 3 male and 3 female plants, see Appendix S3                  |
| section        | identifier for section of stream with similar soil water content in response to experimental irrigation       |
| pfc            | percent soil field capacity of each section on dates indicated below, see Appendix S3                         |
| leaf\_id       | leaf identifier unique to each plant                                                                          |
| year           | year of gas exchange measurements and leaf harvesting                                                         |
| date           | date of gas exchange measurements and leaf harvesting                                                         |
| time           | time of gas exchange measurements and leaf harvesting                                                         |
| chamber\_area  | area of leaf area inside Licor 6400 portable photosynthesis system leaf chamber (mm<sup>2</sup>)              |
| leaf\_area     | area of entire leaf measured in imageJ (mm<sup>2</sup>)                                                       |
| leaf\_mass     | mass of dried leaf (g)                                                                                        |
| lma            | leaf mass to area ratio (g mm<sup>-2</sup>)                                                                   |
| transpiration  | instantaneous transpiration rate (mol H<sub>2</sub>O m<sup>-2</sup> s<sup>-1</sup>) measured with Licor 6400  |
| photosynthesis | instantaneous photosynthetic rate (mol CO<sub>2</sub> m<sup>-2</sup> s<sup>-1</sup>) measured with Licor 6400 |
| conductance    | stomatal conductance (mol H<sub>2</sub>O m<sup>-2</sup> s<sup>-1</sup>) measured with Licor 6400              |
| ci             | intercellular CO<sub>2</sub> concentration measured with Licor 6400                                           |
| wue            | instantaneous water use efficiency calculated as photosynthesis / conductance                                 |

<br>

#### leafP

Metadata for the leaf gas exchange dataset described in Appendix S3.

| field     | description                                                                                             |
|:----------|:--------------------------------------------------------------------------------------------------------|
| plant\_id | plant tag identifier                                                                                    |
| rill      | artificial stream identifier, W vs. E, see Appendix S3                                                  |
| block     | identifier for 3 m subunits of stream containing 3 male and 3 female plants, see Appendix S3            |
| section   | identifier for section of stream with similar soil water content in response to experimental irrigation |
| pfc       | percent soil field capacity of each section on dates indicated below, see Appendix S3                   |
| year      | year of gas exchange measurements and leaf harvesting                                                   |
| P         | leaf percent phosphorus                                                                                 |
