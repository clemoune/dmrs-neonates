function [E, J] = SynthMeasZeppelinCylinder(params,protocol,params_deriv)
%
% camino.m--------------------------------------------------------------
% Normalized diffusion signal simulating a ZeppelinCylinder model.
% 
% [E,J]=SynthMeasZeppelinCylinder(params, protocol,params_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a two compartment diffusion model: 
% Intracellular compartment - parallel cylinders with one radius
% Extracellular compartment - hindered diffusion
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% params - size 6 vector of model parameters in SI units for Cylinder:
%       params(1) - intracellular volume fraction
%       params(2) - intrinsic diffusivity for cylinders = parallel
%       diffusivity for Zeppelin
%       params(3) - perpendicular (hindered) diffusivity
%       params(4) - cylinder radius
%       params(5) - polar angle theta in spherical coordinates desbribing the fibre
% direction
%       params(6) - azimuthal angle phi in spherical coordinates describing the
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

index1 = [2 3 5 6]; % parameters for zeppelin
index2 = [2 4 5 6]; % parameters for cylinder
f = params(1); % volume fraction of cylinder compartment

if nargout == 1           
        E = (1-f)*Zeppelin(params(index1),protocol)+f*Cylinder(params(index2),protocol);
else 
    if nargin <3
        params_deriv = ones(size(params));
    end
    [E1 J1] = Zeppelin(params(index1),protocol,params_deriv(index1));
    [E2 J2] = Cylinder(params(index2),protocol,params_deriv(index2));
    E = (1-f)*E1+f*E2;
    J = zeros(length(E),nnz(params_deriv));
    it = 0;
    if params_deriv(1) it = it+1; J(:,it) = E2 - E1; end % volume fraction
    if params_deriv(2) it = it+1; J(:,it) = (1-f)*J1(:,1) + f*J2(:,1); end % Dpar
    if params_deriv(3) it = it+1; J(:,it) = (1-f)*J1(:,2); end % Dperp
    if params_deriv(4) it = it+1; J(:,it) = f*J2(:,2); end % radius
    if params_deriv(5) it = it+1; J(:,it) = (1-f)*J1(:,3) + f*J2(:,3); end % theta
    if params_deriv(6) it = it+1; J(:,it) = (1-f)*J1(:,4) + f*J2(:,4); end % phi 
    
   
end
end