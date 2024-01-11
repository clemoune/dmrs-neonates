function [E, J] = SynthMeasTortZeppelinCylinderDotBallB0(params,protocol,params_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a TortZeppelinCylinderDotBallB0 model.
% 
% [E,J]=SynthMeasTortZeppelinCylinderDotBallB0(params, protocol,params_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a four compartment diffusion model:
% Intracellular compartment - parallel cylinders with one radius
% Extracellular compartment - hindered diffusion with tortuosity constraint
% Static water - fully restricted diffusion  
% CSF - free diffusion
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% params - size 9 vector of model parameters in SI units for Cylinder:
%       params(1) - intracellular volume fraction as fraction of anisotropic part
%       params(2) - intrinsic diffusivity = parallel diffusivity 
%       params(3) - cylinder radius
%       params(4) - volume fraction of restricted water
%       params(5) - volume fraction of CSF
%       params(6) - free diffusivity in CSF 
%       params(7) - polar angle theta in spherical coordinates desbribing the fibre
% direction
%       params(8) - azimuthal angle phi in spherical coordinates describing the
% fibre direction
%       params(9) - B0 signal intensity without diffusion weighting
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
ficvf = params(1); % volume fraction of cylinder compartment as fraction of anisotropic part
di = params(2);
rad = params(3);
irfrac = params(4);
fiso = params(5);
diso = params(6);
theta = params(7);
phi = params(8);
b0 = params(9);

params_zep = [di (1-ficvf)*di theta phi]; % turtuosity constraint
params_cyl = [di rad theta phi]; % parameters for cylinder
params_ball = diso; % parameters for ball


if nargout == 1           
        E = b0*((1-ficvf)*(1-irfrac)*(1-fiso)*Zeppelin(params_zep,protocol)+ficvf*(1-irfrac)*(1-fiso)*Cylinder(params_cyl,protocol)+ ...
            +irfrac*(1-fiso)*Dot(protocol)+fiso*Ball(params_ball,protocol));
else 
    if nargin <3
        params_deriv = ones(size(params));
    end
    params_deriv_zep = [params_deriv(2) (params_deriv(1) || params_deriv(2)) params_deriv(7:8)];
    params_deriv_cyl = [params_deriv(2) params_deriv(3) params_deriv(7:8)];
    params_deriv_ball = [params_deriv(6)];
    
    [E1 J1] = Zeppelin(params_zep,protocol,params_deriv_zep);
    [E2 J2] = Cylinder(params_cyl,protocol,params_deriv_cyl);
    [E3 J3] = Dot(protocol);
    [E4 J4] = Ball(params_ball,protocol,params_deriv_ball);
    
    E = b0*((1-ficvf)*(1-irfrac)*(1-fiso)*E1+ficvf*(1-irfrac)*(1-fiso)*E2+ ...
        irfrac*(1-fiso)*E3 + fiso*E4);
    J = zeros(length(E),nnz(params_deriv));
    it = 0;
    if params_deriv(1) % ficvf volume fraction, including turtuosity term
        it = it+1; 
        J(:,it) = b0*((1-irfrac)*(1-fiso)*E2 - (1-irfrac)*(1-fiso)*E1 + (1-irfrac)*(1-fiso)*(1-ficvf)*J1(:,2)*(-params_zep(1))); 
    end 
    if params_deriv(2) % di, including turtuosity term
        it = it+1; 
        J(:,it) = b0*((1-irfrac)*(1-ficvf)*(1-fiso)*J1(:,1) + ficvf*(1-irfrac)*(1-fiso)*J2(:,1) + (1-irfrac)*(1-ficvf).^2*(1-fiso)*J1(:,2)) ; 
    end 
    if params_deriv(3) % radius
        it = it+1;
        J(:,it) = b0*(ficvf*(1-irfrac)*(1-fiso)*J2(:,2)); 
    end 
    if params_deriv(4) % irfrac
        it = it+1;
        J(:,it) =b0*(-(1-ficvf)*(1-fiso)*E1 - ficvf*(1-fiso)*E2 + (1-fiso)*E3);
    end
    if params_deriv(5) % fiso
        it = it+1; 
        J(:,it) = b0*(-(1-ficvf)*(1-ficvf)*E1 - ficvf*(1-irfrac)*E2 -irfrac*E3 +E4); 
    end 
    if params_deriv(6) % diso
        it = it+1; 
        J(:,it) = b0*(fiso*J4(:,1)); 
    end 
    if params_deriv(7) % theta
        it = it+1; 
        J(:,it) = b0*((1-irfrac)*(1-ficvf)*(1-fiso)*J1(:,3) + (1-irfrac)*ficvf*(1-fiso)*J2(:,3)); 
    end 
    if params_deriv(8) % phi
        it = it+1; 
        J(:,it) = b0*((1-irfrac)*(1-ficvf)*(1-fiso)*J1(:,4) + (1-irfrac)*ficvf*(1-fiso)*J2(:,4)); 
    end  
    if params_deriv(9) % b0
        it = it+1; 
        J(:,it) = (1-ficvf)*(1-irfrac)*(1-fiso)*E1+ficvf*(1-irfrac)*(1-fiso)*E2+ ...
        irfrac*(1-fiso)*E3 + fiso*E4;
    end  
    
   
end
end