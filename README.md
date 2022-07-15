# Code for 'Wearable reproductive trackers: quantifying a key life history event remotely'
This repository holds the code for a publication that is in review at *Animal Biotelemetry* entitled; 'Wearable reproductive trackers: quantifying a key life history event remotely'

_Authors_:

- Luke Ozsanlav-Harris </a>
- Larry Griffin
- Mitch Weegman
- Lei Cao
- Geoff Hilton
- Stuart Bearhop


## Manuscrip Status
A preprint of our MS can be found on research square here: (https://www.researchsquare.com/article/rs-1489736/v1)

*Major revison* received 31/05/2022

## Code description
- `1_`: Classify incubations using GPS and accelerometer data. Requires each individual to have both data streams and then a set of known breeders from which training values can be calcualted.

## Abstract
Advancements in biologging technology allow terabytes of data to be collected that record the location of individuals but also their direction, speed and acceleration. These multi-stream data sets allow researchers to infer movement and behavioural patterns at high spatiotemporal resolutions and in turn quantify fine-scale changes in state along with likely ecological causes and consequences. The scope offered by such data sets is increasing and there is potential to gain unique insights into a suite of ecological and life history phenomena. We use multi-stream data from global positioning system (GPS) and accelerometer (ACC) devices to quantify breeding events remotely in an Arctic breeding goose. From a training set of known breeders we determine the movement and overall dynamic body acceleration patterns indicative of incubation and use these to classify breeding events in individuals with unknown reproductive status. Given that researchers are often constrained by the amount of biologging data they can collect due to device weights, we carry out a sensitivity analysis. Here we explore the relative merits of GPS vs ACC data and how varying the temporal resolution of the data affects the accuracy of classifying incubation for birds. Classifier accuracy deteriorates as the temporal resolution of GPS and ACC are reduced but the reduction in precision (false positive rate) is larger in comparison to recall (false negative rate). Precision fell to 95.5% whereas recall rarely fell below 98% over all sampling schedules tested. Our dataset could have been reduced by c.95% while maintaining precision and recall >98%. The GPS-only classifier generally outperformed the ACC-only classifier across all accuracy metrics but both performed worse than the combined GPS and ACC classifier. GPS and ACC data can be used to reconstruct breeding events remotely, allowing unbiased, 24-hour monitoring of individuals. Our resampling-based sensitivity analysis of classifier accuracy has important implications with regards to both device design and sampling schedules for study systems where device size is constrained. It will allow researchers with similar aims to optimize device battery, memory usage and lifespan to maximise the ability to correctly quantify life history events. 

