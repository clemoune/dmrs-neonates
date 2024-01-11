function E = SynthMeas_mod(model,protocol)
% function [E J] = SynthMeas(model, protocol)
% returns the diffusion signal and signal jacobian for the given model and protocol
%
% author: Andrada Ianus (a.ianus.11@ucl.ac.uk)
%
% $Id$

% add MM constants if necessary  
    if ~isfield(protocol,'gstep') % in case the user did not call before MM constants
         protocol = MMConstants(model,protocol);          
    end
fun = str2func(['SynthMeas' model.name '_mod']);
 
    E = fun(model.params,protocol);

end