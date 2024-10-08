# Code for 'Wearable reproductive trackers: quantifying a key life history event remotely'
[Link to paper ](https://animalbiotelemetry.biomedcentral.com/articles/10.1186/s40317-022-00298-8) 

[![DOI](https://zenodo.org/badge/548860592.svg)](https://zenodo.org/doi/10.5281/zenodo.10993580)

This repository holds the code for a paper published in *Animal Biotelemetry* entitled; 'Wearable reproductive trackers: quantifying a key life history event remotely'. Briefly it uses biologging data, GPS and accelerometer, to classify incubation events in birds using a simple rule based classification system. We also re-sample our biologging data, altering the temporal resolution, to find the optimum settings for biologging tags to balance battery life and accurate classification of incubation.

## _Authors_

- Luke Ozsanlav-Harris <a itemprop="sameAs" content="https://orcid.org/0000-0003-3889-6722" href="https://orcid.org/0000-0003-3889-6722" target="orcid.widget" rel="noopener" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" alt="ORCID iD icon" target="_blank" style="width:1em;margin-right:.5em;"/></a>
- Larry Griffin
- Mitch Weegman <a itemprop="sameAs" content="https://orcid.org/0000-0003-1633-0920" href="https://orcid.org/0000-0003-1633-0920" target="orcid.widget" rel="me noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" alt="ORCID iD icon" style="width:1em;margin-right:.5em;"/></a>
- Lei Cao
- Geoff Hilton <a itemprop="sameAs" content="https://orcid.org/0000-0001-9062-3030" href="https://orcid.org/0000-0001-9062-3030" target="orcid.widget" rel="me noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" alt="ORCID iD icon" style="width:1em;margin-right:.5em;"/></a>
- Stuart Bearhop <a itemprop="sameAs" content="https://orcid.org/0000-0002-5864-0129" href="https://orcid.org/0000-0002-5864-0129" target="orcid.widget" rel="me noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" alt="ORCID iD icon" style="width:1em;margin-right:.5em;"/></a>


## Applying this method
If users want to apply this method to their tracking dataset then please use the script entitled `Exemplar_Classify_Incubations_GPS+ACC.R`. Unfortunately the data to run this script was not able to be archived but this script provides a cleaned, well commented version of the method to follow along to. If users want help applying this method then please get in contact (see paper for my e-mail). I should also be able to provide a subset of the tracking data upon request to make this Exemplar script run in full and interpretation of the method easier.


## Manuscript Status

*Major revisions* received: 31/05/2022

*Minor revisions* received: 18/07/2022

*Accepted*: 29/07/2022

Manuscript available in *Animal Biotelemetry* [https://animalbiotelemetry.biomedcentral.com/articles/10.1186/s40317-022-00298-8](https://animalbiotelemetry.biomedcentral.com/articles/10.1186/s40317-022-00298-8)

## Code description
- `Exemplar_Classify_Incubations_GPS+ACC.R`: PLEASE USE, Cleaned and commented version of the classification method
- `1) Joint Incubation Classifer.R`: Classify incubations using GPS and accelerometer data. Requires each individual to have both data streams and then a set of known breeders from which training values can be calculated.
- `2) GPS-only Incubation Classifer`: Classify incubations using GPS data only
- `3) Acc-only Incubation Classifer`: Classify incubations using accelerometer data only
- `4) Produce Chapter Plots`: Produce of all of the figures in the manuscript


## Data description
- `TrainingValues/Training values.csv`: The training values created for each of the different sampling schedules. 

## Abstract
Advancements in biologging technology allow terabytes of data to be collected that record the location of individuals but also their direction, speed and acceleration. These multi-stream data sets allow researchers to infer movement and behavioural patterns at high spatiotemporal resolutions and in turn quantify fine-scale changes in state along with likely ecological causes and consequences. The scope offered by such data sets is increasing and there is potential to gain unique insights into a suite of ecological and life history phenomena. We use multi-stream data from global positioning system (GPS) and accelerometer (ACC) devices to quantify breeding events remotely in an Arctic breeding goose. From a training set of known breeders we determine the movement and overall dynamic body acceleration patterns indicative of incubation and use these to classify breeding events in individuals with unknown reproductive status. Given that researchers are often constrained by the amount of biologging data they can collect due to device weights, we carry out a sensitivity analysis. Here we explore the relative merits of GPS vs ACC data and how varying the temporal resolution of the data affects the accuracy of classifying incubation for birds. Classifier accuracy deteriorates as the temporal resolution of GPS and ACC are reduced but the reduction in precision (false positive rate) is larger in comparison to recall (false negative rate). Precision fell to 95.5% whereas recall rarely fell below 98% over all sampling schedules tested. Our dataset could have been reduced by c.95% while maintaining precision and recall >98%. The GPS-only classifier generally outperformed the ACC-only classifier across all accuracy metrics but both performed worse than the combined GPS and ACC classifier. GPS and ACC data can be used to reconstruct breeding events remotely, allowing unbiased, 24-hour monitoring of individuals. Our resampling-based sensitivity analysis of classifier accuracy has important implications with regards to both device design and sampling schedules for study systems where device size is constrained. It will allow researchers with similar aims to optimize device battery, memory usage and lifespan to maximise the ability to correctly quantify life history events. 

