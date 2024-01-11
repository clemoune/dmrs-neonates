function [E, J] = SynthMeasBinghamCylinders(params,protocol,params_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a BinghamCylinders model.
%
% [E,J]=SymthMeasBinghamCylinders(params, protocol,params_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate of cylinders with a Bingham 
% distribution of orientations and a diffusion protocol specified in the input 
%
% Substrate: Cylinders with a Bingham distribution of orientations
% Diffusion pulse sequence: wrapper for various sequences
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% params - size 7 vector of model parameters in SI units for AstroCylinders:
%       params(1) - free diffusivity of the material inside the cylinders.
%       params(2) - radius of the cylinders.
%       params(3) - kappa: concentration parameter of Bingham distribution
%       params(4) - beta: concentration parameter of Bingham distribution
%       params(5) - theta: polar angle in spherical coordinates desbribing the fibre
% direction
%       params(6) - phi: azimuthal angle in spherical coordinates describing the
% fibre direction
%       params(7) - psi: 3rd euler angle of the Bingham's distribution 
% protocol - structure which includes all the information related to the 
%       diffusion protocol and required to generate the signal. 
% params_deriv - a vector of 0s and 1s with the same size as params, indicating
%       which parameters are considered in the Jacobian;
%
%------------------------------------------------------------------------
% This file is part of the camino.m toolbox.
% Copyright (c) 2015, UCL Microstructure Imaging Group (MIG), All rights reserved.
% Distributed under the Modified BSD Licence (see: LICENSE.pdf).
%
% Authors:   
%   Gary Hui Zhang (gary.zhang@ucl.ac.uk)
%   Maira Tariq (maira.tariq.11@ucl.ac.uk)
% 


if nargout == 1           
    E = BinghamCylinders(params,protocol);
else
    if nargin == 3
        [E Jinit] = BinghamCylinders(params,protocol,params_deriv);
        J = zeros(length(E),nnz(params_deriv));
        it = 0;        
        for i = 1:length(params_deriv)
            if params_deriv(i)
                it = it+1;
                J(:,it) = Jinit(:,i);
            end
        end            
    
    else     
    [E J] = BinghamCylinders(params,protocol);    
    end
end
end