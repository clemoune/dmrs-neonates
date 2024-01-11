function [E, J] = SynthMeasTortWatsonZeppelinCylinderBallB0(params,protocol,params_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a TortWatsonZeppelinCylinderDotBall model.
% 
% [E,J]=SynthMeasTortZeppelinCylinderDotBall(params, protocol,params_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a three compartment diffusion model:
% Intracellular compartment - cylinders with one radius and a Watson
% distribution of orientations
% Extracellular compartment - hindered diffusion with tortuosity
% constraint depending on the concentration parameter of the Watson
% distribution
% CSF - free diffusion
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% params - size 9 vector of model parameters in SI units for Cylinder:
%       params(1) - ficvf: intracellular volume fraction as fraction of anisotropic part
%       params(2) - di: free diffusivity of the material inside and outside the
% cylinders
%       params(3) - rad: cylinder radius
%       params(4) - kappa: concentration parameter of the Watson's distribution.
%       params(5) - fiso: volume fraction of isotropic compartment
%       params(6) - diso: diffusivity of isotropic compartment 
%       params(7) - theta: polar angle theta in spherical coordinates desbribing the fibre
% direction
%       params(8) - phi: azimuthal angle phi in spherical coordinates describing the
% fibre direction
%       params(9) - b0: signal intensity without diffusion weighting
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
%   Maira Tariq (maira.tariq.11@ucl.ac.uk)
%   Gary Hui Zhang (gary.zhang@ucl.ac.uk)
% 


if nargin <3
        params_deriv = ones(size(params));
end
    
% Volume fractions
f = params(1); % volume fraction of intra-neurite compartment
fiso = params(5); % volume fraction of the isotropic compartment
S0 = params(9); % b=0 measurement

% tortuosity to work out dperp:
dpar=params(2);
dperp=dpar*(1-f);

% Parameters for WatsonZeppelin - including tortuosity model
% parameters for WatsonZeppelin: [dpar dperp kappa theta phi]
params_WatZep = [dpar dperp params(4) params(7:8)];

% parameters for WatsonCylinders: [di rad kappa theta phi]
params_WatCyl = [params(2) params(3) params(4) params(7:8)] ;

% parameters for Ball: [diso]
params_Ball = [params(6)] ;

if nargout == 1
    % anisotropic signal:
    Eaniso = f*WatsonCylinders(params_WatCyl,protocol) + (1-f)*WatsonZeppelin(params_WatZep,protocol);
    % isotropic signal
    Eiso = Ball(params_Ball,protocol);
    % total signal
    Enorm = fiso*Eiso + (1-fiso)*Eaniso;
else
    
    if nargin <3
        params_deriv = ones(size(params));
    end
    params_deriv_WatZep = [params_deriv(2) (params_deriv(1) || params_deriv(2)) params_deriv(4) params_deriv(7:8)];
    params_deriv_WatCyl = [params_deriv(2) params_deriv(3)  params_deriv(4) params_deriv(7:8)];
    params_deriv_Ball = [params_deriv(6)];

    [E1, J1] = WatsonZeppelin(params_WatZep,protocol,params_deriv_WatZep);
    [E2, J2] = WatsonCylinders(params_WatCyl,protocol,params_deriv_WatCyl);
    [E3, J3] = Ball(params_Ball,protocol,params_deriv_Ball);
    
    % Anisotropic signal
    Eaniso = f*E2 + (1-f)*E1;
    
    % total signal
    Enorm = fiso*E3 + (1-fiso)*Eaniso;
    
    
    % compute the jacobian
    Jnorm = zeros(length(Enorm),nnz(params_deriv));
    it = 0;
    if params_deriv(1) it = it+1; Jnorm(:,it) = (1-fiso)*(E2 - E1 + (1-f)*J1(:,2)*(-params_WatZep(1))); end % volume fraction, including turtuosity term
    if params_deriv(2) it = it+1; Jnorm(:,it) = (1-fiso)* ((1-f)*(J1(:,1)+J1(:,2)) + f*J2(:,1) ); end % Dpar, including turtuosity term
    if params_deriv(3) it = it+1; Jnorm(:,it) = (1-fiso)*f*J2(:,2);  end % radius
    if params_deriv(4) it = it+1; Jnorm(:,it) = (1-fiso)*((1-f)*J1(:,3) + f*J2(:,3)); end % kappa
    if params_deriv(5) it = it+1; Jnorm(:,it) = -(1-f)*E1 - f* E2 + E3; end % fiso
    if params_deriv(6) it = it+1; Jnorm(:,it) = fiso*J3(:,1); end %diso
    if params_deriv(7) it = it+1; Jnorm(:,it) = (1-fiso)*((1-f)*J1(:,4) + f*J2(:,4)); end % theta
    if params_deriv(8) it = it+1; Jnorm(:,it) = (1-fiso)*((1-f)*J1(:,5) + f*J2(:,5)); end % phi
    
    J = Jnorm*S0;
    if params_deriv(9) it = it+1; J(:,it) = Enorm; end %b0
    
end

% multiply the normalised signal by S0
E = Enorm*S0;

end