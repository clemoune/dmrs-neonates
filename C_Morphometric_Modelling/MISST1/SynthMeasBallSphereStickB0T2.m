function [E, J] = SynthMeasBallSphereStickB0T2(params,protocol,params_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating the BallSphereStick model.
% 
% [E,J]=BallSphereAstroStickB0T2(params, protocol,params_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a three compartment diffusion model: 
% Intracellular compartment - sphere
% Extracellular compartment - ball
% Vascular compartment - Stick
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% params - size 7 vector of model parameters in SI units for BallSphereStick:
%       params(1) - f1: intracellular volume fraction
%       params(2) - f2: extracellular volume fraction
%       params(3) - di: diffusivity inside the sphere
%       params(4) - dv: vascular (pseudo) diffusivity
%       params(5) - Rs: sphere radius
%       params(6) - theta: polar angle theta in spherical coordinates desbribing 
% the direction of the vessels
%       params(7) - phi: azimuthal angle phi in spherical coordinates describing 
% the direction of the vessels
%       params(8) - b0: S0 signal without diffusion weighting at T2 = 0
%       params(9) - t2: T2 value
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
% model parameters are: f1, f2, D, Dp, R, theta, phi


index1 = [3 5]; % free diffusivity inside, radius
index2 = [3]; % free diffusivity,
index3 = [4 6 7]; % pseudo diffusion, theta, phi 

fic = params(1); % volume fraction of sphere compartment / intracellular
fee = params(2); % volume fraction of ball compartment / extracellular exravascular
fv = 1-fic-fee; % volume fraction of stick compartment / vascular

S0 = params(8);
T2 = params(9);

if nargout == 1           
        E = fic*Sphere(params(index1),protocol)+fee*Ball(params(index2),protocol) + fv*Stick(params(index3),protocol);
        E = E.*S0.*exp(-protocol.TE'./T2);      
else 
    if nargin <3
        params_deriv = ones(size(params));
    end
    [E1 J1] = Sphere(params(index1),protocol,params_deriv(index1));
    [E2 J2] = Ball(params(index2),protocol,params_deriv(index2));
    [E3 J3] = Stick(params(index3),protocol,params_deriv(index3));
    E = S0.*exp(-protocol.TE'./T2).*(fic*E1+fee*E2+fv*E3);
    J = zeros(length(E),nnz(params_deriv));
    it = 0;

    if params_deriv(1) it = it+1; J(:,it) = S0.*exp(-protocol.TE'./T2).*(E1 - E3); end % volume fraction wrt fic
    if params_deriv(2) it = it+1; J(:,it) = S0.*exp(-protocol.TE'./T2).*(E2 - E3); end % volume fraction wrt fec
    if params_deriv(3) it = it+1; J(:,it) = S0.*exp(-protocol.TE'./T2).*(fic*J1(:,1) + fee*J2(:,1)); end % D
    if params_deriv(4) it = it+1; J(:,it) = S0.*exp(-protocol.TE'./T2).*(fv*J3(:,1)); end % Dp 
    if params_deriv(5) it = it+1; J(:,it) = S0.*exp(-protocol.TE'./T2).*(fic*J1(:,2)); end % R
    if params_deriv(6) it = it+1; J(:,it) = S0.*exp(-protocol.TE'./T2).*(fv*J3(:,2)); end % theta
    if params_deriv(7) it = it+1; J(:,it) = S0.*exp(-protocol.TE'./T2).*(fv*J3(:,3)); end % phi 
    if params_deriv(8) it = it+1; J(:,it) = exp(-protocol.TE'./T2).*(fic*E1+fee*E2+fv*E3); end % S0 
    if params_deriv(9) it = it+1; J(:,it) = S0.*protocol.TE'./T2.^2.*exp(-protocol.TE'./T2).*(fic*E1+fee*E2+fv*E3); end % T2 
end
end