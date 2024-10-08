# Software for Heliospheric Science: IBEX Ribbon Separation
Code associated with IBEX ribbon separation method in Beesley et al. (2023) as part of a suite of "Software for Heliospheric Science" code. See https://github.com/lanl/Theseus for software related to 2-degree Theseus map-making using IBEX data. 

---
## Introduction

NASA’s Interstellar Boundary Explorer (IBEX) satellite collects data on energetic neutral atoms (ENAs) that can provide insight into the heliosphere boundary between our solar system and interstellar space. From these data, researchers construct maps of the ENA intensities (often, expressed in terms of flux) observed by latitude and longitude. The ENA flux observed in these maps is believed to come from at least two distinct sources: one source which manifests as a ribbon of concentrated ENA flux and one source (or possibly several) that manifest as smoothly-varying globally-distributed flux (GDF). Each ENA source type and its corresponding ENA intensity map is of separate scientific interest.

This code demonstrates the ribbon separation method proposed in Beesley et al. (2023). This method takes a pre-make ENA intensity or flux map and corresponding uncertainties as the input and outputs a ribbon/GDF partitioned version of this same map. Beesley et al. (2023) defines a two-stage separation algorithm, where the first stage defines a ribbon "mask" region (i.e., a region possibly containing ribbon flux) and the second stage predicts the GDF flux underneath the ribbon in the mask region. This separation procedure is repeated multiple times for a given input map, and the final separation is obtained as a weighted combination of the best 25% of separations according to a goodness-of-separation heuristic. This code also implements the ribbon center estimation algorithm proposed in Beesley et al. (2023). Uncertainty estimation code is not provided, but the results of applying the Monte Carlo-based uncertainty estimation to simulated and real IBEX data are provided. 

## System requirements

The code is supported on all operating systems for which the requisite downloads (see below) are possible. The example code was tested on a MacBook Pro running macOS Big Sur 11.7.8, using R version 4.2.0.

## Installation

To downloading and install software and packages:
 - R (>= 2.14.0) follow instructions at https://www.r-project.org/

Installation should take less than 15 minutes on a normal desktop computer.

## Demonstration

The ribbon separation method is demonstrated using two real data products: IBEX Science Operations Center (ISOC) public-release 6-degree maps and higher-resolution 2-degree maps using the recently-published Theseus method (https://arxiv.org/abs/2210.12005). These input data products are provided in **input_isoc.csv** and **input_theseus.csv**, respectively. The resulting ribbon separations and corresponding ribbon center estimates (Theseus only) are stored in ecliptic coordinates in **output_isoc_ecliptic.csv** and **output_theseus_ecliptic_part1/2.csv**. Results in ribbon-centric coordinates are provided in **output_isoc_rc.csv** and **output_theseus_rc.csv**. Results using method in Reisenfeld et al. (2021) are stored in **output_isoc_xx_reisenfeld.csv** and **output_theseus_xx_reisenfeld.csv** files. 

Code used to generate these real data ribbon separations is provided in **.R** files **run_separations_isoc.R** and **run_separations_theseus.R**. Users wishing to reproduce these separations locally should expect the code to take ~ 12 hours to run. Although code for performing the simulated data separations is not provided, results of these map separations are provided in **output_simulations_ecliptic.csv** and **output_simulations_rc.csv**. These files also contain the separations using the method in Reisenfeld et al. (2021). Results using method in Reisenfeld et al. (2021) are stored in **output_simulations_xx_reisenfeld.csv** files. 

The **R** scripts **plot_simulated.R** and **plot_realdata.R** will reproduce many of the figures in the manuscript and supplement. 

NOTE: Anyone wishing to use these maps **for space science** (not statistical methodological development) should contact the LANL IBEX team first (email dreisenfeld@lanl.gov). Additional data corrections (e.g., for the Compton-Getting effect) are advised.

## Instructions for use

After R is installed, run **plot_simulated.R** to reproduce many simulated data results, or run **plot_realdata.R** to reproduce all of the real data results in the manuscript and more. Users wishing to reproduce the ribbon separations stored in **output_isoc_xx.csv** and/or **output_theseus_xx.csv** may run **run_separations_isoc.R** and/or **run_separations_theseus.R**. 

## Input Data Details

**input_isoc.csv** is an 6 column data frame containing the ISOC data. These maps were obtained from the ISOC public release 16 (https://ibex.princeton.edu/DataRelease) are represent ram-only yearly maps that have been survival probability corrected but not Compton-Getting effect corrected. Fluxes have been converted to ENA rates for this demonstration:
- lon: is the ecliptic longitude (between 0 and 360)
- lat: is the ecliptic latitude (between -90 and 90)
- ena_rate: estimated ENA rate provided by ISOC
- sd_ena_rate: standard error for the ENA rate
- ESA: energy step (2-6) for map
- Time_Group: year of map

**input_theseus.csv** is an 6 column data frame containing Theseus method (https://arxiv.org/abs/2210.12005) map separations. These maps correspond to combined ram/antiram maps for the first 6 months of each year from 2009-2019. No additional data corrections (e.g., flux, survival probability or Compton-Getting effect) were applied. Code for Theseus map generation is provided at https://github.com/lanl/Theseus. 
- lon: is the ecliptic longitude (between 0 and 360)
- lat: is the ecliptic latitude (between -90 and 90)
- ena_rate: estimated ENA rate generated by Theseus method
- sd_ena_rate: standard error for the ENA rate
- ESA: energy step (2-6) for map
- Time_Group: year of map. 'A' denotes that maps correspond to the first 6 months of each year.

## Output Data Details

**output_isoc_rc.csv** and **output_isoc_ecliptic.csv** are 10 column data frames containing the ribbon-separated ISOC data in ribbon-centered and ecliptic coordinates, respectively. A ribbon center of longitude 221 and latitude 39 was assumed. **output_theseus_rc.csv** and **output_theseus_ecliptic.csv** are also 8 column data frames containing the ribbon-separated Theseus data in ribbon-centered and ecliptic coordinates, respectively. For Theseus, ribbon-centric rotational frames were determined using the ESA 4 total ENA rate map for each year.
- lon: is the ecliptic or ribbon-centered longitude (between 0 and 360)
- lat: is the ecliptic or ribbon-centered (between -90 and 90)
- est_gdf: estimated ENA rate contribution from the globally-distributed flux (GDF)
- se_gdf: standard error of GDF ENA rate
- est_ribbon: estimated ENA rate contribution from the ribbon
- se_ribbon: standard error of ribbon ENA rate
- input_data: ENA rate in the input data, by definition equals est_gdf + est_ribbon (ISOC only)
- input_data_se: standard error of ENA rate in the input data, provided alongside input data (ISOC only)
- ESA: energy step (2-5) for map
- Time_Group: time period of map, where "A" denotes the first 6 months of the year for Theseus maps

**output_simulations_rc.csv** and **output_simulations_ecliptic.csv** are 13 column data frames containing the ribbon-separated simulated data in ribbon-centered and ecliptic coordinates, respectively. The ribbon center was fixed to the simulation truth of of longitude 221.5 and latitude 39.
- lon: is the ecliptic or ribbon-centered longitude (between 0 and 360)
- lat: is the ecliptic or ribbon-centered (between -90 and 90)
- est_gdf: estimated ENA rate contribution from the globally-distributed flux (GDF)
- se_gdf: standard error of GDF ENA rate
- est_ribbon: estimated ENA rate contribution from the ribbon
- se_ribbon: standard error of ribbon ENA rate
- input_data: ENA rate in the input data, by definition equals est_gdf + est_ribbon
- input_data_se: standard error of ENA rate in the input data, provided alongside input data
- truth_gdf: simulation truth GDF-only ENA rate
- truth_ribbon: simulation truth ribbon-only ENA rate
- ESA: ESA (2-6) of real data map used to generate the simulated binned data structure used by Thesesus. Equals "noesa" for simulation truth maps (i.e., not estimated by Theseus).
- MAP_CLASS: simulated ribbon profile associated with map (weakscattering, spatialretention, or varyingprofile)
- MAP_STAT: indicates were map corresponds to simulation truth, a Theseus map estimate using simulated data with standard follow-up times, or a Theseus map estimate using simulated data with triple the observed data follow-up times (truth, amountstandard, or amount3x)

**ibex_palette** is a 4 column data frame containing the pixel color information to plot sky maps:
- red: red numeric value
- green: green numeric value
- blue: blue numeric value
- hex: hex value
  
## Attribution

If you use these ribbon separations or code in your research work, please cite the following paper:

Lauren J. Beesley, Dave Osthus, Kelly R. Moran, Madeline A. Ausdemore, Grant David Meadors,  Thomas K. Kim, Sung Jun Noh, Nehpreet K. Walia, Paul H. Janzen, Eric J. Zirnstein, Brian P. Weaver, Daniel B. Reisenfeld. Statistical methods for partitioning ribbon and globally-distributed flux using data from the Interstellar Boundary Explorer. [Under consideration at JASA Applications and Case Studies](https://arxiv.org/abs/2302.03089).

---
Copyright 2023 for **O4627**

This program is Open-Source under the BSD-3 License.
 
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and
the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions
and the following disclaimer in the documentation and/or other materials provided with the
distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse
or promote products derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
