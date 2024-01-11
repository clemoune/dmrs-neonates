function [E J] = SynthMeasNormalCylinders(params,protocol,params_deriv)
  
    if nargout == 1           
        E = NormalCylinders(params,protocol);
    else
        if nargin == 3
        [E Jinit] = NormalCylinders(params,protocol,params_deriv);
        J = zeros(length(E),nnz(params_deriv));
        it = 0;        
        for i = 1:length(params_deriv)
            if params_deriv(i)
                it = it+1;
                J(:,it) = Jinit(:,i);
            end
        end            
        else     
        [E J] = NormalCylinders(params,protocol);    
        end
    end
end
