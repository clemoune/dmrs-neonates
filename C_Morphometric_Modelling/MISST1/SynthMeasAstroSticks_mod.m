function E = SynthMeasAstroSticks_mod(params,protocol)
% function [E J] = SynthMeasAstroSticks(model, protocol, params_deriv)
% returns the diffusion signal and signal jacobian for the given protocol
% and AstroSticks tissue model
%
% author: Andrada Ianus (a.ianus.11@ucl.ac.uk)
%
% $Id$

fun_name = ['AstroSticks_' protocol.pulseseq '_mod'];

fun = str2func(fun_name);

E = fun(params,protocol);
end