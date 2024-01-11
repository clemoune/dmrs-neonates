function [E, J] = SynthMeasZeppelinCylinderStick(params,protocol,params_deriv)
%
% camino.m--------------------------------------------------------------
% Normalized diffusion signal simulating a ZeppelinCylinderStick model.
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
%       params(5) - stick volume fraction
%       params(6) - polar angle theta in spherical coordinates desbribing the fibre
% direction
%       params(7) - azimuthal angle phi in spherical coordinates describing the
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

index1 = [2 3 6 7]; % parameters for zeppelin
index2 = [2 4 6 7]; % parameters for cylinder
index3 = [2 6 7]; % parameters for stick
f1 = params(1); % volume fraction of cylinder compartment
f2 = params(5);

if nargout == 1           
        E = (1-f1-f2)*Zeppelin(params(index1),protocol)+f1*Cylinder(params(index2),protocol) +...
            f2*Stick(params(index2),protocol);
else 
    if nargin <3
        params_deriv = ones(size(params));
    end
    [E1 J1] = Zeppelin(params(index1),protocol,params_deriv(index1));
    [E2 J2] = Cylinder(params(index2),protocol,params_deriv(index2));
    [E3 J3] = Stick(params(index3),protocol,params_deriv(index3));
    E = (1-f1-f2)*E1+f1*E2 + f2*E3;
    J = zeros(length(E),nnz(params_deriv));
    it = 0;
    if params_deriv(1) it = it+1; J(:,it) = E2 - E1; end % volume fraction 1 
    if params_deriv(2) it = it+1; J(:,it) = (1-f1-f2)*J1(:,1) + f1*J2(:,1) + f2*J3(:,1); end % Dpar
    if params_deriv(3) it = it+1; J(:,it) = (1-f1)*J1(:,2); end % Dperp
    if params_deriv(4) it = it+1; J(:,it) = f1*J2(:,2); end % radius
    if params_deriv(5) it = it+1; J(:,it) = E3 - E1; end % volume fraction of sticks
    if params_deriv(6) it = it+1; J(:,it) = (1-f1-f2)*J1(:,3) + f1*J2(:,3) +f2*J2(:,2); end % theta
    if params_deriv(7) it = it+1; J(:,it) = (1-f1-f2)*J1(:,4) + f1*J2(:,4)+f2*J2(:,3); end % phi 
    
   
end
end