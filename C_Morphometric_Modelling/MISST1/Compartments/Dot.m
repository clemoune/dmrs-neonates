function [E,J]=Dot(protocol)
%
% camino.m--------------------------------------------------------------
% Diffusion signal for fully restricted diffusion
% 
% [E,J]=Dot(protocol)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for fully restricted diffusion
% Substrate: fully restricted diffusion
% Diffusion pulse sequence: -
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
E = ones(length(protocol.G),1);
if nargout >1
    J = [];
end