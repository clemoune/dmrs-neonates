function [E, J] = SynthMeasTortZeppelinStick(params,protocol,params_deriv)
%
% camino.m--------------------------------------------------------------
% Normalized diffusion signal simulating a TortZeppelinStick model.
% 
% [E,J]=SynthMeasTortZeppelinStick(params, protocol,params_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a two compartment diffusion model:
% Intracellular compartment - parallel sticks 
% Extracellular compartment - hindered diffusion with tortuosity constraint
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% params - size 4 vector of model parameters in SI units:
%       params(1) - f: intracellular volume fraction 
%       params(2) - di: intrinsic diffusivity = parallel diffusivity 
%       params(3) - theta: polar angle theta in spherical coordinates desbribing the fibre
% direction
%       params(4) - phi: azimuthal angle phi in spherical coordinates describing the
% fibre direction
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

f = params(1); % volume fraction of cylinder compartment

params_zep = [params(2) params(2)*(1-f) params(3:4)]; % turtuosity constraint
params_stick = params(2:end); % parameters for stick

if nargout == 1           
        E = (1-f)*Zeppelin(params_zep,protocol)+f*Stick(params_stick,protocol);
else 
    if nargin <3
        params_deriv = ones(size(params));
    end
    params_deriv_zep = [params_deriv(2) (params_deriv(1) || params_deriv(2)) params_deriv(4:5)];
    params_deriv_stick = params_deriv(2:end);
    [E1 J1] = Zeppelin(params_zep,protocol,params_deriv_zep);
    [E2 J2] = Stick(params_stick,protocol,params_deriv_stick);
    E = (1-f)*E1+f*E2;
    J = zeros(length(E),nnz(params_deriv));
    it = 0;
    if params_deriv(1) it = it+1; J(:,it) = E2 - E1 + (1-f)*J1(:,2)*(-params_zep(1)); end % volume fraction, including turtuosity term
    if params_deriv(2) it = it+1; J(:,it) = (1-f)*J1(:,1) + f*J2(:,1) + (1-f).^2*J1(:,2) ; end % Dpar, including turtuosity term
    if params_deriv(4) it = it+1; J(:,it) = (1-f)*J1(:,3) + f*J2(:,2); end % theta
    if params_deriv(5) it = it+1; J(:,it) = (1-f)*J1(:,4) + f*J2(:,3); end % phi 
    
   
end
end