function [E, J] = SynthMeasTortZeppelinCylinderBall(params,protocol,params_deriv)
%
% camino.m--------------------------------------------------------------
% Normalized diffusion signal simulating a TortZeppelinCylinderBall model.
% 
% [E,J]=SynthMeasTortZeppelinCylinderBall(params, protocol,params_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a three compartment diffusion model:
% Intracellular compartment - parallel cylinders with one radius
% Extracellular compartment - hindered diffusion with tortuosity constraint
% CSF - free diffusion
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% params - size 7 vector of model parameters in SI units for Cylinder:
%       params(1) - ficvf: intracellular volume fraction as fraction of anisotropic part
%       params(2) - di: intrinsic diffusivity = parallel diffusivity 
%       params(3) - rad: cylinder radius
%       params(4) - fiso: volume fraction of CSF
%       params(5) - diso: free diffusivity in CSF 
%       params(6) - theta: polar angle theta in spherical coordinates desbribing the fibre
% direction
%       params(7) - phi: azimuthal angle phi in spherical coordinates describing the
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

ficvf = params(1); % volume fraction of cylinder compartment as fraction of anisotropic part
di = params(2);
rad = params(3);
fiso = params(4);
diso = params(5);
theta = params(6);
phi = params(7);

params_zep = [di (1-ficvf)*di theta phi]; % turtuosity constraint
params_cyl = [di rad theta phi]; % parameters for cylinder
params_ball = diso; % parameters for ball


if nargout == 1           
        E = (1-ficvf)*(1-fiso)*Zeppelin(params_zep,protocol)+ficvf*(1-fiso)*Cylinder(params_cyl,protocol)+ fiso*Ball(params_ball,protocol);
else 
    if nargin <3
        params_deriv = ones(size(params));
    end
    params_deriv_zep = [params_deriv(2) (params_deriv(1) || params_deriv(2)) params_deriv(6:7)];
    params_deriv_cyl = [params_deriv(2) params_deriv(3) params_deriv(6:7)];
    params_deriv_ball = [params_deriv(5)];
    
    [E1 J1] = Zeppelin(params_zep,protocol,params_deriv_zep);
    [E2 J2] = Cylinder(params_cyl,protocol,params_deriv_cyl);
    [E3 J3] = Ball(params_ball,protocol,params_deriv_ball);
     E = (1-ficvf)*(1-fiso)*E1+ficvf*(1-fiso)*E2+ fiso*E3;
    J = zeros(length(E),nnz(params_deriv));
    it = 0;
    if params_deriv(1) % ficvf volume fraction, including turtuosity term
        it = it+1; 
        J(:,it) = (1-fiso)*E2 - (1-fiso)*E1 + (1-fiso)*(1-ficvf)*J1(:,2)*(-params_zep(1)); 
    end 
    if params_deriv(2) % di, including turtuosity term
        it = it+1; 
        J(:,it) = (1-ficvf)*(1-fiso)*J1(:,1) + ficvf*(1-fiso)*J2(:,1) + (1-ficvf).^2*(1-fiso)*J1(:,2) ; 
    end 
    if params_deriv(3) % radius
        it = it+1;
        J(:,it) = ficvf*(1-fiso)*J2(:,2); 
    end 
    if params_deriv(4) % fiso
        it = it+1; 
        J(:,it) = -(1-ficvf)*E1 - ficvf* E2 + E3; 
    end 
    if params_deriv(5) % diso
        it = it+1; 
        J(:,it) = fiso*J3(:,1); 
    end 
    if params_deriv(6) % theta
        it = it+1; 
        J(:,it) = (1-ficvf)*(1-fiso)*J1(:,3) + ficvf*(1-fiso)*J2(:,3); 
    end 
    if params_deriv(7) % phi
        it = it+1; 
        J(:,it) = (1-ficvf)*(1-fiso)*J1(:,4) + ficvf*(1-fiso)*J2(:,4); 
    end  
    
   
end
end