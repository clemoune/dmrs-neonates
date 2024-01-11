function strings = GetParameterStrings(modelname)
%
% function strings = GetParameterStrings(modelname)
%
% Given an input modelname, this function returns the names of the model
% parameters
%
% author: Andrada Ianus (a.ianus.11@ucl.ac.uk), Gary Hui Zhang (gary.zhang@ucl.ac.uk)
%
% Volume fraction inside axons (ficvf)
% Intrinsic diffusion inside restricted compartments (di)
% Radius of the axons (cylinders) (rad)
% Radius of the cells (spheres) (rads)
% Volume fraction of CSF (Ball) (fiso)
% Diffusion constant in CSF (Ball) (diso)
% Angles determining fibre direction (theta, phi)
% Angles for tensor (theta, phi, psi)

if (strcmp(modelname, 'AstroCylinders'))
    strings = {'di', 'rad'};

elseif (strcmp(modelname, 'AstroSticks'))
    strings = {'di'};  
    
elseif (strcmp(modelname, 'Ball'))
    strings = {'diso'}; 
    
elseif (strcmp(modelname, 'Cylinder'))
    strings = {'di', 'rad','theta', 'phi'};
    
elseif (strcmp(modelname, 'FiniteAstroCylinders'))
    strings = {'di', 'rad','lcyl'};    
    
elseif (strcmp(modelname, 'FiniteCylinder'))
    strings = {'di', 'rad','lcyl','theta', 'phi'};     
    
elseif (strcmp(modelname, 'Sphere'))
    strings = {'di','rads'};    
    
elseif (strcmp(modelname, 'SphereAstroCylinders'))
    strings = {'di','rads','rad'};    
    
elseif (strcmp(modelname, 'Stick'))
    strings = {'di','theta', 'phi'};  
    
elseif (strcmp(modelname, 'Tensor'))    
    strings = {'d1','d2','d3','theta', 'phi','psi'};     
    
elseif (strcmp(modelname, 'TortZeppelinCylinder'))    
    strings = {'ficvf','di', 'rad','theta', 'phi'};    
    
elseif (strcmp(modelname, 'TortZeppelinCylinderBall'))    
    strings = {'ficvf','di','rad','fiso','diso','theta', 'phi'};      
    
elseif (strcmp(modelname, 'ZeppelinCylinder'))    
    strings = {'ficvf','di','dh', 'rad','theta', 'phi'};  
    
elseif (strcmp(modelname, 'Zeppelin'))    
    strings = {'di','dh','theta', 'phi'};  

elseif (strcmp(modelname, 'BallSphereAstroSticksAstroZeppelins'))    
    strings = {'fic','fee','fstr', 'dic', 'dvasc', 'Rs', 'dstr_par', 'dstr_ort'};  

elseif (strcmp(modelname, 'BallSphereAstroSticks'))    
    strings = {'fic','fee', 'dic', 'dvasc', 'Rs'};  

elseif (strcmp(modelname, 'BallNormalSpheresAstroSticks'))    
    strings = {'fic','fee', 'dic', 'dvasc', 'Rmean', 'Rvar'};      
elseif (strcmp(modelname, 'BallNormalSpheresAstroSticks_T2weighted'))    
    strings = {'fic','fee', 'dic', 'dvasc', 'Rmean', 'Rvar'};      
else
    error(['Parameter strings yet to be defined for this model:', modelname]);
end

