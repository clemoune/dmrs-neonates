function [E, J] = SynthMeasKurtosis(params,protocol,params_deriv)
%
% camino.m--------------------------------------------------------------
% Free isotropic normalized diffusion signal.
% 
% [E,J]=Ball(params, protocol,params_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for free isotropic diffusion
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% params - size 1 vector of model parameters in SI units for Ball:
%       params(1) - diso: free diffusivity of the material 
%       params(2) - K: kurtosis from exp(-bD + b^2 D^2 K/6)    
% protocol - structure which includes all the information related to the 
%        diffusion protocol and required to generate the signal.
% params_deriv - a vector of 0s and 1s with the same size as x, indicating
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
    E = Kurtosis(params,protocol);
else
    if nargin == 3
    [E Jinit] = Kurtosis(params,protocol,params_deriv);
        J = zeros(length(E),nnz(params_deriv));
        it = 0;        
        for i = 1:length(params_deriv)
            if params_deriv(i)
                it = it+1;
                J(:,it) = Jinit(:,i);
            end
        end            
    else     
    [E J] = Kurtosis(params,protocol);    
    end
end
end