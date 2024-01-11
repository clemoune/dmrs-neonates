function [E J] = SynthMeasNormalFiniteAstroCylinders(params,protocol,params_deriv)
  
    if nargout == 1           
        E = NormalFiniteAstroCylinders(params,protocol);
    else
        if nargin == 3
        [E Jinit] = NormalFiniteAstroCylinders(params,protocol,params_deriv);
        J = zeros(length(E),nnz(params_deriv));
        it = 0;        
        for i = 1:length(params_deriv)
            if params_deriv(i)
                it = it+1;
                J(:,it) = Jinit(:,i);
            end
        end            
        else     
        [E J] = NormalFiniteAstroCylinders(params,protocol);    
        end
    end
end