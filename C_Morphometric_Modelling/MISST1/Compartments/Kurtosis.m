function [E,J]=Kurtosis(x, protocol, x_deriv)
%
% camino.m--------------------------------------------------------------
% Free isotropic diffusion signal.
% 
% [E,J]=Kurtosis(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for an isotropic kurtosis model
% Substrate: isotropic diffusion and kurtosis
% Diffusion pulse sequence: wrapper for various sequences
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% x - size 1 vector of model parameters in SI units for Ball:
%       x(1) - free diffusivity of the material 
%       x(2) -isotropic kurtosis of the material 
% protocol - structure which includes all the information related to the 
%        diffusion protocol and required to generate the signal.
%       The basic information for a PGSE sequence is the following:
%       protocol.grad_dirs - is the gradient direction for each measurement.
%           It has size [N 3] where N is the number of measurements.
%       protocol.G - gradient strength, size [1 N]
%       protocol.delta - pulse separation, size [1 N]
%       protocol.smalldel - pulse duration, size [1 N]
%       Other diffusion sequences might have additional fields.
% x_deriv - a vector of 0s and 1s with the same size as x, indicating
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


D=x(1);
K = x(2);
bval = GetBvalues(protocol)';


% Signal
logE=(-bval.*D + bval.^2.*D.^2.*K./6);
E = exp(logE);                   
  

 % Compute the Jacobian matrix; computed numerically
if(nargout>1)
     J = zeros(length(E), 1);      
    if nargin < 3               
        dEdD = (-bval+2*bvals.^2.*D.*K/6).*exp(-bval.*D + bval.^2.*D.^2.*K./6); 
        dEdK = (bvals.^2.*D.^2/6).*exp(-bval.*D + bval.^2.*D.^2.*K./6); 
        J(:,1) = dEdD;
        J(:,2) = dEdK;
    else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian         
        
            if x_deriv(1) ~= 0                             
                dEdD = (-bval+2*bvals.^2.*D.*K/6).*exp(-bval.*D + bval.^2.*D.^2.*K./6); 
                
                J(:,1) = dEdD;
            end  
             if x_deriv(2) ~= 0  
                dEdK = (bvals.^2.*D.^2/6).*exp(-bval.*D + bval.^2.*D.^2.*K./6); 
            J(:,2) = dEdK;
             end

    end   
 end
