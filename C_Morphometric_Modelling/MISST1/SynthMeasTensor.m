function [E, J] = SynthMeasTensor(params,protocol,params_deriv)
%
% camino.m--------------------------------------------------------------
% Normalized diffusion signal simulating the Tensor model.
% 
% [E,J]=Tensor(params, protocol,params_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a diffusion tensor model and a 
% diffusion protocol specified in the input
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% params - size 6 vector of model parameters in SI units for Stick:
%       params(1) - d1: diffusivity of the material along the first direction
%       params(2) - d2: diffusivity of the material along the second direction
%       params(3) - d3: diffusivity of the material along the third direction
%       params(4) - theta: angle from the z direction of the main direction of the
% tensor
%       params(5) - phi: azymuthal angle of the main direction of the
% tensor
%       params(6) - psi: angle of rotation of the second direction 
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
    E = Tensor(params,protocol);
else
    if nargin == 3
    [E Jinit] = Tensor(params,protocol,params_deriv);
        J = zeros(length(E),nnz(params_deriv));
        it = 0;        
        for i = 1:length(params_deriv)
            if params_deriv(i)
                it = it+1;
                J(:,it) = Jinit(:,i);
            end
        end            
    else     
    [E J] = Tensor(params,protocol);    
    end
end
end