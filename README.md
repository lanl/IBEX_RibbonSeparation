# IBEX_RibbonSeparation
Minimal working example for IBEX ribbon separation method in Beesley et al. (2023)

---
## Introduction

NASA’s Interstellar Boundary Explorer (IBEX) satellite collects data on energetic neutral atoms (ENAs) that can provide insight into the heliosphere boundary between our solar system and interstellar space. From these data, researchers construct maps of the ENA intensities (often, expressed in terms of flux) observed by latitude and longitude. The ENA flux observed in these maps is believed to come from at least two distinct sources: one source which manifests as a ribbon of concentrated ENA flux and one source (or possibly several) that manifest as smoothly-varying globally-distributed flux (GDF). Each ENA source type and its corresponding ENA intensity map is of separate scientific interest.

This code demonstrates the ribbon separation method proposed in Beesley et al. (2023). This method takes a pre-make ENA intensity or flux map and corresponding uncertainties as the input and outputs a ribbon/GDF partitioned version of this same map with corresponding standard errors, which propagate uncertainty in both the input map and uncertainty in the map partitioning. Beesley et al. (2023) defines a two-stage separation algorithm, where the first stage defines a ribbon "mask" region (i.e., a region possibly containing ribbon flux) and the second stage predicts the GDF flux underneath the ribbon in the mask region. This separation procedure is repeated multiple times for a given input map, and the final separation is obtained as a weighted combination of the best 25% of separations according to a goodness-of-separation heuristic. This code also implements the ribbon center estimation algorithm proposed in Beesley et al. (2023).

## System requirements

The code is supported on all operating systems for which the requisite downloads (see below) are possible. The example code was tested on a MacBook Pro running macOS Monterey 12.6.3, using R version 4.2.2.

## Installation

To downloading and install software and packages:
 - R (>= 2.14.0) follow instructions at https://www.r-project.org/

Installation should take less than 15 minutes on a normal desktop computer.

## Demonstration

The ribbon separation method is demonstrated using three separate data products: IBEX Science Operations Center (ISOC) public-release 6-degree maps, higher-resolution 2-degree maps using the recently-published Theseus method (https://arxiv.org/abs/2210.12005), and simulated data maps. These input data products are provided in **inputs_isoc.RDS**, **inputs_theseus.RDS**, and **inputs_simulated.RDS**, respectively. 

The **.R** files **run_separations_isoc.R**, **run_separations_theseus.R**, and **run_separations_simulated.R** will run the ribbon separation algorithm on each of the three data products. These ribbon separations may take several hours to run, but the results are pre-populated. 



The resulting ribbon separations are stored in **separations_isoc.RDS**, **separations_theseus.RDS**, and **separations_simulated.RDS**. The script **run_centers_simulated.R** will apply the proposed ribbon center estimation method to simulated data, and the outputs are stored in **centers_simulated.RDS**.

The example code **generate_figures_simulated.R** will reproduce Figures XX-XX in the manuscript. Script **generate_figures_realdata.R** will reproduce Figures XX-XX. 

NOTE: Anyone wishing to use these maps **for space science** (not statistical methodological development) should contact the LANL IBEX team first (email dreisenfeld@lanl.gov). Additional data corrections (e.g., for the Compton-Getting effect) are advised.

## Reproducing results in Beesley et al. (2023)



## Instructions for use

After R is installed, run **generate_figures_simulated.R** to reproduce Figures XX-XX, or run **generate_figures_realdata.R** to reproduce Figures XX-XX. 

Users may need to setwd('Path/to/Theseus/Directory/') in line 34 of both **.R** files.


## Data Details

### data_illustration.RDS

There are **5** data products in **data_illustration.RDS**.

**Xobs** is an 11 column data frame containing the simulated binned direct event data:
- obs_id: is a unique label for the binned direct event data, from 1 to the number of rows of Xobs
- ecliptic_lon: is the ecliptic longitude (between 0 and 360)
- ecliptic_lat: is the ecliptic latitude (between -90 and 90)
- ecliptic_lon_center: is the ecliiptic longitude in "nose centered" frame and is used for plotting purposes
- x, y, and z: are the sperical coordinates for the ecliptic longitude and latitude and are used for PPR and GAM regression
- orbit_number: is the orbit number from the data collection and is used in the bootstrap sampling
- counts: are the number of direct events
- time: is the exposure time (in seconds)
- background: is the background rate (background particles per second)

**Xpix** is an 8 column data frame containing the sky map pixel locations and other relevant pixel information:
- pix_id: is a unique label for each pixel, from 1 to the number of rows of Xpix (for a 2 degree map, that's 16,200)
- ecliptic_lon: is the ecliptic longitude (between 0 and 360)
- ecliptic_lat: is the ecliptic latitude (between -90 and 90)
- ecliptic_lon_center: is the ecliiptic longitude in "nose centered" frame and is used for plotting purposes
- x, y, and z: are the sperical coordinates for the ecliptic longitude and latitude and are used for PPR and GAM regression
- wt_pix: is proportional to the size of the pixel (area on a unit sphere) where larger pixels are near ecliptic_latitude 0 and smaller pixels are near ecliptic latitudes -90 and 90 and is used in `glmnet()`

**Xpsf** is a nrow(Xobs) by nrow(Xpix) sparse matrix where each entry is non-neagative and each row sums to 1. It is used to relate the unblurred sky map to the blurred binned direct events.

**ibex_palette** is a 4 column data frame containing the pixel color information to plot sky maps:
- red: red numeric value
- green: green numeric value
- blue: blue numeric value
- hex: hex value

### data_results.RDS 

There are **2** data products in **data_results.RDS**.

**ibex_data** is an 11 column data frame containing the binned direct event data for all ESA 4 "A" maps between 2010 and 2021:
- esa: is the energy step for the IBEX data. Only ESA 4 is provided
- map: is the 6-month map id (e.g., "2013A"). Only "A" maps are provided
- ecliptic_lon: is the ecliptic longitude (between 0 and 360)
- ecliptic_lat: is the ecliptic latitude (between -90 and 90)
- ecliptic_lon_center: is the ecliiptic longitude in "nose centered" frame and is used for plotting purposes
- x, y, and z: are the sperical coordinates for the ecliptic longitude and latitude
- counts: are the number of direct events
- time: is the exposure time (in seconds)
- background: is the background rate (background particles per second)

**real_sky_maps** is an 8 column data frame containing the ISOC and Theseus sky map estimates with UQ for all ESA 4 "A" maps between 2010 and 2021:
- esa: is the energy step for the IBEX data. Only ESA 4 is provided
- map: is the 6-month map id (e.g., "2013A"). Only "A" maps are provided
- ecliptic_lon: is the ecliptic longitude (between 0 and 360)
- ecliptic_lat: is the ecliptic latitude (between -90 and 90)
- ecliptic_lon_center: is the ecliiptic longitude in "nose centered" frame and is used for plotting purposes
- method: the sky map estimation method. Either "theseus" of "isoc"
- quantity: the type of response, either "mean" or "ci" for 95% confidence interval width
- value: either the ENA rate (ENAs/sec when q

  
## Attribution

If you use these ribbon separations or code in your research work, please cite the following paper:

Lauren J. Beesley, Dave Osthus, Kelly R. Moran, Madeline A. Ausdemore, Grant David Meadors, Paul H. Janzen, Eric J. Zirnstein, Brian P. Weaver, Daniel B. Reisenfeld. Statistical methods for partitioning ribbon and globally-distributed flux using data from the Interstellar Boundary Explorer. [Under consideration at JASA Applications and Case Studies](https://arxiv.org/abs/2302.03089).

---
Copyright 2023 for **CO4627**

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
