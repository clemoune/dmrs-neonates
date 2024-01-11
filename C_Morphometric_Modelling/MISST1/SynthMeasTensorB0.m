function [E, J] = SynthMeasTensorB0(params,protocol,params_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating the Tensor model.
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
% params - size 7 vector of model parameters in SI units for Stick:
%       params(1) - d1: diffusivity of the material along the first direction
%       params(2) - d2: diffusivity of the material along the second direction
%       params(3) - d3: diffusivity of the material along the third direction
%       params(4) - theta: angle from the z direction of the main direction of the
% tensor
%       params(5) - phi: azymuthal angle of the main direction of the
% tensor
%       params(6) - psi: angle of rotation of the second direction 
%       params(7) - b0: signal intensity without diffusion weighting
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
    E = params(end).*Tensor(params(1:end-1),protocol);
else
    S0 = params(end);
    if nargin == 3
    [Enorm Jnorm] = Tensor(params(1:end-1),protocol,params_deriv(1:end-1));
        J = zeros(length(Enorm),nnz(params_deriv));
        it = 0;        
        for i = 1:length(params_deriv)-1
            if params_deriv(i)
                it = it+1;
                J(:,it) = S0*Jnorm(:,i);
            end
        end
        if params_deriv(end)
            J(:,end) = Enorm;
        end
    else     
    [Enorm Jnorm] = Tensor(params(1:end-1),protocol); 
    J = [S0*Jnorm Enorm];
    end
    E = S0*Enorm;
end
end