function [E J]=CuboidSym(x,protocol,x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a cuboid/ensemble of cuboids with 2 equal
% sides
% 
% [E,J]=CuboidSym(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate consisting of a cuboid/ 
% an ensemble of cuboids with 2 equal sides and a diffusion protocol 
% specified in the input
% Substrate: cuboid / ensemble of cuboids with 2 equal sides
% Diffusion pulse sequence: wrapper for various sequences
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% x - size 3 vector of model parameters in SI units for CuboidSym:
%       x(1) - free diffusivity of the material inside the cuboids.
%       x(2) - length of the cuboid (lx = ly)
%       x(3) - eccentricity (ratio between height of the cuboid lz and 
%       length lx)
% protocol - structure which includes all the information related to the 
%        diffusion protocol and required to generate the signal.    
% x_deriv - a vector of 0s and 1s with the same size as x, indicating
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
if isfield(protocol,'approx')
    approx = protocol.approx;
else
    approx = 'GPD';
end

fun_name = ['CuboidSym_' approx '_' protocol.pulseseq];
fun = str2func(fun_name);

if nargin <3
    if nargout == 1
        E = fun(x,protocol);
    else
        [E J] = fun(x,protocol);

    end
else
    if nargout == 1
        E = fun(x,protocol, x_deriv);
    else
        [E J] = fun(x,protocol, x_deriv);
    end
end