function [E J] = SynthMeasNormalFiniteCylinders(params,protocol,params_deriv)
  
    if nargout == 1           
        E = NormalFiniteCylinders(params,protocol);
    else
        if nargin == 3
        [E Jinit] = NormalFiniteCylinders(params,protocol,params_deriv);
        J = zeros(length(E),nnz(params_deriv));
        it = 0;        
        for i = 1:length(params_deriv)
            if params_deriv(i)
                it = it+1;
                J(:,it) = Jinit(:,i);
            end
        end            
        else     
        [E J] = NormalFiniteCylinders(params,protocol);    
        end
    end
end