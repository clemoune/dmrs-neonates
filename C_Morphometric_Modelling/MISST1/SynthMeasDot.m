function [E, J] = SynthMeasDot(~,protocol,~)
%
% camino.m--------------------------------------------------------------
% Fully restriced diffusion signal.
% 
% [E,J]=SynthMeasDot(~, protocol,~)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a one compartment model cosisting of
% fully restricted diffusion signal
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% protocol - structure which includes all the information related to the 
%        diffusion protocol and required to generate the signal.    
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


if nargout == 1           
    E = Dot(protocol);
else
    E = Dot(protocol);
    J  = [];  
end
end