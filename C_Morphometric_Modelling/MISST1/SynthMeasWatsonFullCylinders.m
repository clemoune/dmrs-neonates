function [E, J] = SynthMeasWatsonFullCylinders(params,protocol,params_deriv)
% The function returns the diffusion signal and signal jacobian for the
% given substrate and protocol
%
%
% $Id$


if nargout == 1           
    E = WatsonCylinders(params,protocol);
else
    if nargin == 3
        [E Jinit] = WatsonFullCylinders(params,protocol,params_deriv);
        J = zeros(length(E),nnz(params_deriv));
        it = 0;        
        for i = 1:length(params_deriv)
            if params_deriv(i)
                it = it+1;
                J(:,it) = Jinit(:,i);
            end
        end            
    
    else     
    [E J] = WatsonFullCylinders(params,protocol);    
    end
end
end