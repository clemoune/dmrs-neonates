# dmrs-neonates

Code coming with the bioRxiv: https://www.biorxiv.org/content/10.1101/2023.10.16.562599v1

## A_dMRS_processing/

This folder contains the files related to the processing of the spectra: 
- coil combination based on water signal
- outlier removal based on kurtosis tensor estimation of the upfield signal mean Mup for high b values data
- phase & frequency correction prior to summing
- eddy current suppression & water residual removal whenever necessary


## B_Astrosticks_Spheres_model/

This folder contains the files related to the analytical modelling of the data ("astrosticks + spheres" model):
- .csv files containing the signal attenuations and ADCs at long diffusion times coming from the LCModel quantification
- Jupyter notebook containing the code (DIVE is needed as indicated in the header of the notebook): dMRS_neonate_astrosticks-and-spheres.ipynb
- Outputs/ contains the .csv files produced by the analysis (hence the **"Fitting Data"** step of the notebook can be ignored)


