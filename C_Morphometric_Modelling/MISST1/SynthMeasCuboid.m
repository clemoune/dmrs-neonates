function [E, J] = SynthMeasCuboid(params,protocol,params_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a cuboid/ensemble of cuboids model.
% 
% [E,J]=SymthMeasCuboid(params, protocol,params_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate consisting of a cuboid/ 
% an ensemble of cuboids and a diffusion protocol specified in the input
% Substrate: cuboid / ensemble of cuboids
% Diffusion pulse sequence: wrapper for various sequences
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% params - size 4 vector of model parameters in SI units for Cuboid:
%       params(1) - free diffusivity of the material inside the cuboids.
%       params(2) - length of the cuboid (x dir)
%       params(3) - width of the cuboid (y dir)
%       params(4) - height of the cuboid (z dir) 
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


if nargout == 1           
    E = Cuboid(params,protocol);
else
    if nargin == 3
    [E Jinit] = Cuboid(params,protocol,params_deriv);
        J = zeros(length(E),nnz(params_deriv));
        it = 0;        
        for i = 1:length(params_deriv)
            if params_deriv(i)
                it = it+1;
                J(:,it) = Jinit(:,i);
            end
        end            
    else     
    [E J] = Cuboid(params,protocol);    
    end
end
end