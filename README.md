# dmrs-neonates

Code coming with the bioRxiv: https://www.biorxiv.org/content/10.1101/2023.10.16.562599v1

## A_dMRS_processing/

This folder contains the files related to the processing of the spectra: 
- coil combination based on water signal
- outlier removal based on kurtosis tensor estimation of the upfield signal mean Mup for high b values data
- phase & frequency correction prior to summing
- eddy current suppression & water residual removal whenever necessary
