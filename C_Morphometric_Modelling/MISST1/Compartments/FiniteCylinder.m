function [E J]=FiniteCylinder(x,protocol,x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating the FiniteCylinder compartment.
% 
% [E,J]=AstroCylinders(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate consisting of isotropically 
% oriented finite cylinders and a diffusion protocol specified in the input
% Substrate: Parallel impermeable finite cylinders
% Diffusion pulse sequence: wrapper for various sequences
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% x - size 5 vector of model parameters in SI units for FiniteCylinder:
%       x(1) - free diffusivity of the material inside the cylinders.
%       x(2) - radius of the cylinders.
%       x(3) - eccentricity (ratio between length and diameter)
%       x(4) - polar angle theta in spherical coordinates desbribing the fibre
% direction
%       x(5) - azimuthal angle phi in spherical coordinates describing the
% fibre direction
% protocol - structure which includes all the information related to the 
%       diffusion protocol and required to generate the signal. 
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

fun_name = ['FiniteCylinder_' approx '_' protocol.pulseseq];
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