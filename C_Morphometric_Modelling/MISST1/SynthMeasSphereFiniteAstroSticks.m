function [E, J] = SynthMeasSphereFiniteAstroSticks(params,protocol,params_deriv)
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
%       params(1) - ficvf: volume fraction of cylinders 
%       params(2) - di: free diffusivity of the material inside the cylinders.    
%       params(3) - L: Length of the sticks  
%       params(4) - rads: radius of the sphere.  
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
ficvf = params(1); % volume fraction of cylinders
di = params(2);
L = params(3);
rads = params(4);

index1 = [2 4];
index2 = [2 3];
params_sphere = params(index1);
params_cyl = params(index2);

if nargout == 1           
    E = ficvf* FiniteAstroSticks(params_cyl,protocol) + (1-ficvf)*Sphere(params_sphere,protocol);
else
    if nargin <3
        params_deriv = ones(size(params));
    end
    [E1, J1] = Sphere(params(index1),protocol,params_deriv(index1));
    [E2, J2] = FiniteAstroSticks(params(index2),protocol,params_deriv(index2));
    E = ficvf*E2+(1-ficvf)*E1;
    J = zeros(length(E),nnz(params_deriv));
    it = 0;

    if params_deriv(1) it = it+1; J(:,it) = E2 - E1; end % volume fraction wrt fic
    if params_deriv(2) it = it+1; J(:,it) = (1-ficvf)*J1(:,1) + ficvf*J2(:,1) ; end % di
    if params_deriv(3) it = it+1; J(:,it) = ficvf*J2(:,2); end % L
    if params_deriv(4) it = it+1; J(:,it) = (1-ficvf)*J1(:,2); end % R
  
end
end