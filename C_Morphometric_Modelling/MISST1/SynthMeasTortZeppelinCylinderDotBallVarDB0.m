function [E, J] = SynthMeasTortZeppelinCylinderDotBallVarDB0(params,protocol,params_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a TortZeppelinCylinderDotBallB0 model.
% 
% [E,J]=SynthMeasTortZeppelinCylinderDotBallB0(params, protocol,params_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a four compartment diffusion model:
% Intracellular compartment - parallel cylinders with one radius
% Extracellular compartment - hindered diffusion with tortuosity constraint
% Static water - fully restricted diffusion  
% CSF - free diffusion
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% params - size 9 vector of model parameters in SI units for Cylinder:
%       params(1) - intracellular volume fraction as fraction of anisotropic part
%       params(2) - intrinsic diffusivity = parallel diffusivity 
%       params(3) - cylinder radius
%       params(4) - volume fraction of restricted water
%       params(5) - volume fraction of CSF
%       params(6) - free diffusivity in CSF 
%       params(7) - alpha: coefficient for linear dependence of diffusivity
%       params(8) - polar angle theta in spherical coordinates desbribing the fibre
% direction
%       params(9) - azimuthal angle phi in spherical coordinates describing the
% fibre direction
%       params(10) - B0 signal intensity without diffusion weighting
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
ficvf = params(1); % volume fraction of cylinder compartment as fraction of anisotropic part
di0 = params(2);
rad = params(3);
irfrac = params(4);
fiso = params(5);
diso = params(6);
alpha = params(7);
theta = params(8);
phi = params(9);
b0 = params(10);

if strcmp(protocol.pulseseq,'TWOGSE')
    NT = floor(protocol.smalldel.*protocol.omega./pi+0.00000000001);
    di = di0.*ones(length(protocol.smalldel),1);
    di(NT>1) = di0 + alpha.*protocol.omega(NT>1)';
else error('VarD model implemented for TWOGSE only');
end

params_zep = {di, (1-ficvf)*di0, theta, phi}; % turtuosity constraint
params_cyl = {di, rad, theta, phi}; % parameters for cylinder
params_ball = diso; % parameters for ball


if nargout == 1           
        E = b0*((1-ficvf)*(1-irfrac)*(1-fiso)*Zeppelin(params_zep,protocol)+ficvf*(1-irfrac)*(1-fiso)*Cylinder(params_cyl,protocol)+ ...
            +irfrac*(1-fiso)*Dot(protocol)+fiso*Ball(params_ball,protocol));
else 
    error('derivatives not implemented')
end
 