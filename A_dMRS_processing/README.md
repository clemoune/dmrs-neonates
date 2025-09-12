## Processing dMRS raw data

This folder contains the scripts used to process the Bruker raw data (Paravision 360.1).

The main script is *Protocol_reading_pups.m*, it calls the others scripts *BSB_...*.  

The script *BSB_HighB_Pups.m* processes the data acquired at high b values, it includes the DKI-based outliers suppression.
The script *BSB_LongTM_Pups.m* processes the data acquired at long diffusion times. The call to the script *BSB_PhaseEstimated.m* is needed the high b-values & the long diffusion times data are not treated simulataneously.

The folder **support_functions** contains the functions to read the acquisition parameters of each Bruker folder (*read_parameters.m*), to handle the conversion from Bruker format to matlab matrixes (*load_array_FID2*), to extract signals from different coils (*coil_decomposition.m*)to correct for phase and frequency individual repetitions (all the *phaser...* functions, and associated *cout...* functions), to remove the water residual (*svdfid.m*), and to run the DKI fit (*diffmodelfit.m*).

The folder **MatDataLight** contains the raw data acquired at high b values (simply coil combined and ordered by direction & b-value) as well as the preprocessed data at long diffusion times (generated with *BSB_LongTM_Pups.mat*). If you load any of the *ReorderedLight.mat*, you can run steps 3 to 6C in the *BSB_HighB_Pups.m* script for one animal at a given age.

For storage purposes, the folder is hosted on Zenodo (10.5281/zenodo.17107499).

