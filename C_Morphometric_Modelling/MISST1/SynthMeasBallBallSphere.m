function [E, J] = SynthMeasBallBallSphere(params,protocol,params_deriv)
%
% camino.m--------------------------------------------------------------
% Normalized diffusion signal simulating the BallBallSphere model.
% 
% [E,J]=BallBallSphere(params, protocol,params_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a three compartment diffusion model: 
% Intracellular compartment - sphere
% Extracellular compartment - ball
% Vascular compartment - ball
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% params - size 6 vector of model parameters in SI units for BallBallSphere:
%       params(1) - f1: intracellular volume fraction
%       params(2) - f2: extracellular volume fraction
%       params(3) - di: diffusivity inside the sphere
%       params(4) - rads: sphere radius
%       params(5) - dv: vascular (pseudo) diffusivity
%       params(6) - diso: extracellular diffusivity
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
%   Colleen Bailey (colleen.bailey@ucl.ac.uk)
%   Andrada Ianus (a.ianus.11@ucl.ac.uk)
%   Daniel C. Alexander (d.alexander@ucl.ac.uk)
%  



index1 = [3 4]; % intracellular diffusivity, radius
index2 = [6]; % extracellular diffusivity,
index3 = [5]; % pseudo diffusion

fic = params(1); % volume fraction of sphere compartment / intracellular
fee = params(2); % volume fraction of ball compartment / extracellular exravascular
fv = 1-fic-fee; % volume fraction of stick compartment / vascular

if nargout == 1           
        E = fic*Sphere(params(index1),protocol)+fee*Ball(params(index2),protocol) + fv*Ball(params(index3),protocol);
              
else 
    if nargin <3
        params_deriv = ones(size(params));
    end
    [E1 J1] = Sphere(params(index1),protocol,params_deriv(index1));
    [E2 J2] = Ball(params(index2),protocol,params_deriv(index2));
    [E3 J3] = Ball(params(index3),protocol,params_deriv(index3));
    E = fic*E1+fee*E2+fv*E3;
    J = zeros(length(E),nnz(params_deriv));
    it = 0;

    if params_deriv(1) it = it+1; J(:,it) = E1 - E3; end % volume fraction wrt fic
    if params_deriv(2) it = it+1; J(:,it) = E2 - E3; end % volume fraction wrt fec
    if params_deriv(3) it = it+1; J(:,it) = fic*J1(:,1) + fee*J2(:,1); end % di
    if params_deriv(4) it = it+1; J(:,it) = fv*J3(:,1); end % Dp 
    if params_deriv(5) it = it+1; J(:,it) = fic*J1(:,2); end % R
  
end
end