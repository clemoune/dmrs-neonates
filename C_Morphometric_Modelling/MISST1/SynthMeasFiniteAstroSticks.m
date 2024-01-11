function [E, J] = SynthMeasFiniteAstroSticks(params,protocol,params_deriv)
% function [E J] = SynthMeasAstroCylinders(model, protocol, params_deriv)
% returns the diffusion signal and signal jacobian for the given protocol
% and FiniteAstroCylinders tissue model
%
% author: Andrada Ianus (a.ianus.11@ucl.ac.uk)
%
% $Id$

%fun_name = ['FiniteAstroCylinders_' protocol.pulseseq];
fun_name = 'FiniteAstroCylinders';

params(2) = 1e-12;

fun = str2func(fun_name);
if nargout == 1           
    E = fun(params,protocol);
else
    if nargin == 3
        [E Jinit] = fun(params,protocol,params_deriv);
        J = zeros(length(E),nnz(params_deriv));
        it = 0;        
        for i = 1:length(params_deriv)
            if params_deriv(i)
                it = it+1;
                J(:,it) = Jinit(:,i);
            end
        end            
    
    else     
    [E J] = fun(params,protocol);    
    end
end
end