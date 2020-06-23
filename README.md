# mostela-dieren

Repository for data and script for the manuscript: Jeroen Mos & Tim R. Hofmeester (accepted) The Mostela: an adjusted camera-trapping device as a promising non-invasive tool to study and monitor small mustelids. Mammal Research.

<b>Abstract</b>

In spite of their potential important role in shaping small mammal population dynamics, weasel (Mustela nivalis) and stoat (Mustela erminea) are understudied due to the difficulty of detecting these species. Furthermore, their conservation status in many countries is unknown due to lack of monitoring techniques. There is thus an important need for a method to detect these small mustelids. In this study, we tested the efficiency of a recently developed camera-trapping device, the Mostela, as a new technique to detect mustelids in a study area near Dieren, the Netherlands. We placed Mostelas in linear landscape features, and other microhabitats thought to be frequently visited by weasels, from March – October 2017 and February – October 2018. We tested for yearly and monthly differences in site use and detectability, as well as the effect of entrance tube size, using an occupancy-modelling framework. We found large seasonal differences in site use and detectability of weasels with the highest site use in June to October and highest detection probability in August and September. Detection probability was approximately two times higher for Mostelas with a 10 cm entrance tube compared to 8 cm. Furthermore, we were able to estimate activity patterns based on the time of detection, identify the sex in most detections (69.5%) and distinguish several individuals. Concluding, the Mostela seems promising as a non-invasive monitoring tool to study the occurrence and ecology of small mustelids. Further development of individual recognition from images would enable using the Mostela for density estimates applying capture-recapture models.

<b>Data</b>

The occ- files contain the detection history and covariate data as used in the manuscript.
The weasel-activity-data contains the time of detection data needed for the activity pattern analysis.

The analysis is split in two R-scripts:
Occ-Mostela-Dieren-MS.R contains the description of the model and code to run the occupancy model using the different occ- files.
Activity-analysis-Mostela-Dieren-MS.R contains the code needed to run the activity analysis using the weasel-activity-data file.

<b>Repository</b>

The published version of the data and code can be found on Zenodo

<a href="https://zenodo.org/badge/latestdoi/242999580"><img src="https://zenodo.org/badge/242999580.svg" alt="DOI"></a>
