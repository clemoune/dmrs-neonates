function [E, J] = SynthMeasTensorDotB0(params,protocol,params_deriv)
%
% camino.m--------------------------------------------------------------
% Normalized diffusion signal simulating a TensorDotB0 model.
% 
% [E, J] = SynthMeasTensorDotB0(params,protocol,params_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a two compartment diffusion model: 
% Stationary water compartment - dot
% Anisotropic diffusion - tensor
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% params - size 6 vector of model parameters in SI units for Cylinder:
%       params(1) - volume fraction of stationary water
%       params(2) - d1: diffusivity of the material along the first direction
%       params(3) - d2: diffusivity of the material along the second direction
%       params(4) - d3: diffusivity of the material along the third direction
%       params(5) - theta: angle from the z direction of the main direction of the
% tensor
%       params(6) - phi: azymuthal angle of the main direction of the
% tensor
%       params(7) - psi: angle of rotation of the second direction 
%       params(8) - b0: signal intensity without diffusion weighting
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

index1 = 2:7; % parameters for tensor
f = params(1); % volume fraction of stationary water
B0 = params(end);

if nargout == 1           
        E = ((1-f)*Tensor(params(index1),protocol)+f).*B0;
else 
    if nargin <3
        params_deriv = ones(size(params));
    end
    [E1 J1] = Tensor(params(index1),protocol,params_deriv(index1));
    E = ((1-f)*E1+f).*B0;
    J = zeros(length(E),nnz(params_deriv));
    it = 0;
    if params_deriv(1) it = it+1; J(:,it) = (1-E1).*B0; end % volume fraction
    if params_deriv(2) it = it+1; J(:,it) = (1-f)*J1(:,1).*B0 ; end % Dpar
    if params_deriv(3) it = it+1; J(:,it) = (1-f)*J1(:,2).*B0 ; end % Dperp1
    if params_deriv(4) it = it+1; J(:,it) = (1-f)*J1(:,3).*B0 ; end % Dperp2   
    if params_deriv(5) it = it+1; J(:,it) = (1-f)*J1(:,4).*B0 ; end % Dperp2 
    if params_deriv(6) it = it+1; J(:,it) = (1-f)*J1(:,5).*B0 ; end % Dperp2 
    if params_deriv(7) it = it+1; J(:,it) = (1-f)*J1(:,6).*B0 ; end % Dperp2 
    if params_deriv(8) it = it+1; J(:,it) = ((1-f)*E1+f); end % B0

    
   
end
end