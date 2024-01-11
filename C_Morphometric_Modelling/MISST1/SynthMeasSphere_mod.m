function E = SynthMeasSphere_mod(params, protocol)
% function [E J] = SynthMeasSphere(model, protocol, params_deriv)
% returns the diffusion signal and signal jacobian for the given protocol
% and Sphere tissue model
%
% author: Andrada Ianus (a.ianus.11@ucl.ac.uk)
%
% $Id$

fun_name = ['Sphere_' protocol.pulseseq '_mod'];
fun = str2func(fun_name);
E = fun(protocol);
end