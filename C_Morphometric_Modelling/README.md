## Morphometric modelling and analysis of real cells 

This folder contains multiple sub-folders. 

### Morphometric modelling 

The main scripts are: 
1) ```MachineLearningPipeline/A_Fit_Morphometric_Model.m```: runs the morphometric modelling 
2) ```MachineLearningPipeline/B_Visualisation_morphometry_results.m```: displays the results

These scripts rely on the ```MISST_1``` and ```Data``` folders too. 

If you use any part of the code written in ```MachineLearningPipeline```, please cite the to following papers: 

Palombo, Marco, et al. "New paradigm to assess brain cell morphology by diffusion-weighted MR spectroscopy in vivo." Proceedings of the National Academy of Sciences 113.24 (2016): 6671-6676.

and 

Palombo, Marco, et al. "A generative model of realistic brain cells with application to numerical simulation of the diffusion-weighted MR signal." NeuroImage 188 (2019): 391-402.

The folder ```MISST_1``` is an older/specific version of MISST apckage from UCL (http://mig.cs.ucl.ac.uk/index.php?n=Tutorial.MISST). It is needed to correctly the main scripts. 

The folder ```Data``` contains the processed and quantified data for metabolites at each time point. 

### Analyse real cells
Folder ```analyse_real_cells```

This folder contains the code and data used to generate the analysis of Figure 6.

