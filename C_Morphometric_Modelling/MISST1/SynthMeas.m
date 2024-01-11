function [E, J] = SynthMeas(model,protocol)
% function [E J] = SynthMeas(model, protocol)
% returns the diffusion signal and signal jacobian for the given model and protocol
%
% author: Andrada Ianus (a.ianus.11@ucl.ac.uk)
%
% $Id$

% add MM constants if necessary 
if strcmp(protocol.pulseseq,'GEN')    
    if ~isfield(protocol,'gstep') % in case the user did not call before MM constants
         protocol = MMConstants(model,protocol);          
    end
end
fun = str2func(['SynthMeas' model.name]);
if nargout == 1    
    E = fun(model.params,protocol);
else
    if ~isfield(model,'params_deriv')
         [E, J] = fun(model.params,protocol);
    else
    [E, J] = fun(model.params,protocol,model.params_deriv);
    end
end
end