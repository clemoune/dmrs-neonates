function [E, J] = SynthMeasBallSphereIndependent(params,protocol,params_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating the BallSphere model with independent di and diso.
% 
% [E,J]=BallSphereAstroStick(params, protocol,params_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a three compartment diffusion model: 
% Intracellular compartment - sphere
% Extracellular compartment - ball
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% params - size 7 vector of model parameters in SI units for BallSphereStick:
%       params(1) - f1: intracellular volume fraction
%       params(2) - di: diffusivity inside the sphere
%       params(3) - Rs: sphere radius
%       params(4) - diso diffusivity outside the sphere
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
% model parameters are: f1, di, D, R, diso
%Adapted to BallSphere by James Breen-Norris (james.breen-norris.13@ucl.ac.uk) 


index1 = [2 3]; % free diffusivity inside, radius
index2 = [4]; % free diffusivity,


fic = params(1); % volume fraction of sphere compartment / intracellular
fee = 1 - fic; % volume fraction of ball compartment / extracellular exravascular



if nargout == 1           
        E = fic*Sphere(params(index1),protocol)+fee*Ball(params(index2),protocol);
              
else 
    if nargin <3
        params_deriv = ones(size(params));
    end
    [E1 J1] = Sphere(params(index1),protocol,params_deriv(index1));
    [E2 J2] = Ball(params(index2),protocol,params_deriv(index2));
    E = fic*E1+fee*E2;
    J = zeros(length(E),nnz(params_deriv));
    it = 0;

    if params_deriv(1) it = it+1; J(:,it) = E1 - E2; end % volume fraction wrt fic
    if params_deriv(2) it = it+1; J(:,it) = fic*J1(:,1); end % Di
    if params_deriv(3) it = it+1; J(:,it) = fic*J1(:,2); end % R
    if params_deriv(4) it = it+1; J(:,it) = fee*J2(:,1); end % Diso

  
end
end
