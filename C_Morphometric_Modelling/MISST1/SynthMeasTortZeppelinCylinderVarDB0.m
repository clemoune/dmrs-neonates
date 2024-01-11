function [E, J] = SynthMeasTortZeppelinCylinderVarDB0(params,protocol,params_deriv)
%
% camino.m--------------------------------------------------------------
% Normalized diffusion signal simulating a TortZeppelinCylinder model.
% 
% [E,J]=SynthMeasTortZeppelinCylinderVarD(params, protocol,params_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a two compartment diffusion model:
% Intracellular compartment - parallel cylinders with one radius
% Extracellular compartment - hindered diffusion with tortuosity constraint
% diffusivity depends on frequency
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% params - size 5 vector of model parameters in SI units for Cylinder:
%       params(1) - f: intracellular volume fraction 
%       params(2) - di: intrinsic diffusivity = parallel diffusivity 
%       params(3) - rad: cylinder radius
%       params(4) - alpha: coefficient for linear dependence of diffusivity
%       params(5) - theta: polar angle theta in spherical coordinates desbribing the fibre
% direction
%       params(6) - phi: azimuthal angle phi in spherical coordinates describing the
% fibre direction
% protocol - structure which includes all the information related to the 
%        diffusion protocol and required to generate the signal.    
% params_deriv - a vector of 0s and 1s with the same size as params, indicating
%       which parameters are considered in the Jacobian;
%
%------------------------------------------------------------------------
% This file is part of the camino.m toolbox.
% Copyright (c) 2015, UCL Microstructure Imaging Group (MIG), All rights reserved.
% Distributed under the Modified BSD Licence (see: LICENSE.pdf).
%
% Authors:
%   Andrada Ianus (a.ianus.11@ucl.ac.uk)
%   Daniel C. Alexander (d.alexander@ucl.ac.uk)
% 

f = params(1); % volume fraction of cylinder compartment

if strcmp(protocol.pulseseq,'TWOGSE')
    D0_par = params(2);
    alpha = params(4);
    NT = floor(protocol.smalldel.*protocol.omega./pi+0.00000000001);
    Dpar = D0_par.*ones(length(protocol.smalldel),1);
    Dpar(NT>1) = D0_par + alpha.*protocol.omega(NT>1)';
else error('VarD model implemented for TWOGSE only');
end

params_zep = {Dpar, D0_par*(1-f), params(5),params(6)}; % turtuosity constraint
params_cyl = {Dpar, params(3), params(5),params(6)}; % parameters for cylinder
b0=params(7);
if nargout == 1           
        E = b0*((1-f)*Zeppelin(params_zep,protocol)+f*Cylinder(params_cyl,protocol));
else 
    error ('Derivative not implemented') 

end
end