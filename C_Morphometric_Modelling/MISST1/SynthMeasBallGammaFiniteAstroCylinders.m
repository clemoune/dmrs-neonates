function [E, J] = SynthMeasBallGammaFiniteAstroCylinders(params,protocol,params_deriv)
%
% camino.m--------------------------------------------------------------
% Normalized diffusion signal simulating a BallGammaFiniteAstroCylinders model.
% 
% [E,J]=SynthMeasBallGammaFiniteAstroCylinders(params, protocol,params_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a two compartment diffusion model: 
% Intracellular compartment - GammaFiniteAstroCylinders
% Extracellular compartment - isotropic diffusion
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% params - size 6 vector of model parameters in SI units for Cylinder:
%       params(1) - intracellular volume fraction
%       params(2) - intrinsic diffusivity for cylinders 
%       params(3) - hindered diffusivity of the isotropic extracellular
%       space < intrinsic diffusivity
%       params(4) - cylinder mean radius
%       params(5) - shape of the gamma distribution
%       params(6) - cylinder eccentricity (ratio between length and diameter)
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

index1 = 3; % parameters for ball
index2 = [2 4 5 6]; % parameters for GammaFiniteAstroCylinders
f = params(1); % volume fraction of GammaFiniteAstroCylinders

if nargout == 1           
        E = (1-f)*Ball(params(index1),protocol)+f*GammaFiniteAstroCylinders(params(index2),protocol);
else 
    if nargin <3
        params_deriv = ones(size(params));
    end
    [E1 J1] = Ball(params(index1),protocol,params_deriv(index1));
    [E2 J2] = GammaFiniteAstroCylinders(params(index2),protocol,params_deriv(index2));
    E = (1-f)*E1+f*E2;
    J = zeros(length(E),nnz(params_deriv));
    it = 0;
    if params_deriv(1) it = it+1; J(:,it) = E2 - E1; end % volume fraction
    if params_deriv(2) it = it+1; J(:,it) = f*J2(:,1); end % Di
    if params_deriv(3) it = it+1; J(:,it) = (1-f)*J1; end % Dh
    if params_deriv(4) it = it+1; J(:,it) = f*J2(:,2); end %  mean radius
    if params_deriv(5) it = it+1; J(:,it) = f*J2(:,3); end %  shape of the gamma distribution
    if params_deriv(6) it = it+1; J(:,it) = f*J2(:,4); end %  eccentricity 
    
   
end
end