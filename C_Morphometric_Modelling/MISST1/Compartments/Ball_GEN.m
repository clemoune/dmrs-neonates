function [E,J]=Ball_GEN(x, protocol, x_deriv)
% function [E J] = Ball_GEN(model, protocol, params_deriv)
% returns the diffusion signal and signal jacobian for the given protocol
% and Ball tissue model
%
% Substrate: free isotropic diffusion
% Pulse sequence: Generalized gradient spin echo.
% Signal approximation: -
%
% x is the list of model parameters in SI units:
% x(1) is the free diffusivity of the material inside and outside the cylinders.
%
% protocol is a structure containing information related to the diffusion
% sequence
% x_deriv is a vector of 0s and 1s indicating which derivatives are
% calculated
%
% author: Andrada Ianus (a.ianus.11@ucl.ac.uk), Ivana Drobnjak (i.drobnjak@ucl.ac.uk), Daniel
% C. Alexander (d.alexander@ucl.ac.uk)
%
% $Id$

% Model parameters
dIso = x(1);

if ~isfield(protocol,'GENj')
    protocol.GENj=1;
end

if ~isfield(protocol,'complex')
  protocol.complex='real';
end

tau=protocol.tau;
G=protocol.G;
[M totsizeG]=size(G);

GAMMA = 2.675987E8;
% Calculating the parallel & perpendicular signal is similar  (parallel plane & cylindrical restriction )

K=totsizeG/3;

% get Ndir directions on a spehere

Bval=zeros(M,1);
           



    for ind=1:K
     Fx=sum(G(:,1:3:3*(ind-1)+1)*tau,2);
     Fy=sum(G(:,2:3:3*(ind-1)+2)*tau,2);
     Fz=sum(G(:,3:3:3*ind)*tau,2);   
        
     Bval=Bval+(Fx.^2+Fy.^2+Fz.^2)*tau; 

    end



Bval=GAMMA^2*Bval;
E=exp(-Bval.*dIso);

%disp(ePerp)
if strcmp(protocol.complex,'complex')
  E=[real(E);imag(E)];
elseif strcmp(protocol.complex,'real')
  E=real(E);
elseif strcmp(protocol.complex,'abs')
  E=abs(E);
end





% Compute the Jacobian matrix
if(nargout>1)
     J = zeros(length(E),1);   
    if nargin < 3 
            
        dEtdD = -Bval.*exp(-Bval*dIso);           
        J(:,1) = dEtdD;
            
    else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
         
            if x_deriv(1) ~= 0  
                dEtdD = -Bval.*exp(-Bval*dIso);           
                J(:,1) = dEtdD;                  
            end             
        
    end   
   
end