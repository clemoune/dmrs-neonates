function [E,J]=WatsonZeppelin(x, protocol,x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a WatsonZeppelin compartment.
%
% [E,J]=WatsonZeppelin(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a compartment of anisotropic hindered 
% diffusion with diffusivities dependent on Watson distribution and 
% protocol specified in the input 
%
% Substrate: Anisotropic hindered diffusion with defusivities depending on
% the Watson distribution. 
% Diffusion pulse sequence: wrapper for various sequences
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% x - size 4 vector of model parameters in SI units for AstroCylinders:
%       x(1) - initial parallel diffusivity
%       x(2) - initial perpendicular (hindered) diffusivity
%       x(3) is the angle the polar angle theta in spherical coordinates desbribing the fibre
% direction
%       x(4) is the azimuthal angle phi in spherical coordinates describing the
% fibre direction
% protocol - structure which includes all the information related to the 
%       diffusion protocol and required to generate the signal. 
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
%   Gary Hui Zhang (gary.zhang@ucl.ac.uk)
%   Maira Tariq (maira.tariq.11@ucl.ac.uk)
%   Andrada Ianus (a.ianus.11@ucl.ac.uk)
% 

if length(x) ~= 5
    error('the model should have exactly five parameters');
end

dPar = x(1);
dPerp = x(2);
kappa = x(3);
theta = x(4);
phi = x(5);

% get the equivalent diffusivities
if (nargout == 1)
    dw = WatsonHinderedDiffusionCoeff(dPar, dPerp, kappa);
else
    [dw, Jdw] = WatsonHinderedDiffusionCoeff(dPar, dPerp, kappa);
end

xZep = [dw(1) dw(2) theta phi];


% compute the signal from the zeppelin compartment, given the diffusion
% coefficients from Watson distribution

% fun = str2func(['Zeppelin_' protocol.pulseseq]);

% E = fun(xZep, protocol);
E = Zeppelin(xZep,protocol);

% else
%     xZep_deriv = [x_deriv(1) x_deriv(2) x_deriv(4) x_deriv(5)];
%     [E, JZep] = fun(xZep, protocol, xZep_deriv);

% end

% % Compute the Jacobian matrix
% if(nargout>1)
%     % Construct the jacobian matrix.
%     J = zeros(size(E, 1), 3);
%     if x_deriv(1) ~= 0 J(:,1) = JZep(:,1)*Jdw(1,1) + JZep(:,2)*Jdw(2,1); end;
%     if x_deriv(2) ~= 0 J(:,2) = JZep(:,1)*Jdw(1,2) + JZep(:,2)*Jdw(2,2); end;
%     if x_deriv(3) ~= 0 J(:,3) = JZep(:,1)*Jdw(1,3) + JZep(:,2)*Jdw(2,3); end;
% end

% Compute the Jacobian matrix; computed numerically
if(nargout>1)
    dx = 0.00001;
    J = zeros(length(E), 5);
    if nargin < 3
        for i = 1:5; % compute the derivatives for all model parameters
            xpert = x;
            xpert(i) = xpert(i)*(1+dx);
            Epert = WatsonZeppelin(xpert, protocol);
            dEdx = (Epert - E)/(xpert(i)*dx);
            J(:,i) = dEdx;
        end
    else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
        
        for i = 1:length(x_deriv);
            if x_deriv(i) ~= 0
                xpert = x;
                xpert(i) = xpert(i)*(1+dx);
                Epert = WatsonZeppelin(xpert, protocol);
                dEdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEdx;
            end
        end
    end
end


