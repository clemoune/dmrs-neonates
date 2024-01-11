function [E, J] = SynthMeasFiniteAstroCylinders(params,protocol,params_deriv)
%
% camino.m--------------------------------------------------------------
% Normalized diffusion signal simulating the FiniteAstroCylinder model.
% 
% [E,J]=SynthMeasFiniteAstroCylinder(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for the a model of isotropically oriented
% finite cylinders and protocol specified in the input
% 
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% params - size 3 vector of model parameters in SI units for AstroCylinders:
%       params(1) - di: free diffusivity of the material inside the cylinders.
%       params(2) - rad: radius of the cylinders.       
%       params(3) - ecc: eccentricity (ratio between cylinder length and diameter)   
% protocol - structure which includes the information related to the 
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
    E = FiniteAstroCylinders(params,protocol);
else
    if nargin == 3
    [E Jinit] = FiniteAstroCylinders(params,protocol,params_deriv);
        J = zeros(length(E),nnz(params_deriv));
        it = 0;        
        for i = 1:length(params_deriv)
            if params_deriv(i)
                it = it+1;
                J(:,it) = Jinit(:,i);
            end
        end            
    else     
    [E J] = FiniteAstroCylinders(params,protocol);    
    end
end
end