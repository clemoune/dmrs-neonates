function [E J]=Stick_GEN(x,protocol,x_deriv)
% function [E J] = Stick_GEN(model, protocol, params_deriv)
% returns the diffusion signal and signal jacobian for the given protocol
% and Stick tissue model
%
% Substrate: unidirectional diffusion
% Pulse sequence: Generalized gradient spin echo.
% Signal approximation: -
%
% x is the list of model parameters in SI units:
% x(1) is the free diffusivity
% x(2) is the angle from the z direction of the sticks
% x(3) is the azymuthal angle from x direction
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
if ~isfield(protocol,'GENj')
    protocol.GENj=1;
end

if ~isfield(protocol,'complex')
  protocol.complex='real';
end
% Model parameters
dPar = x(1);
theta = x(2);
phi = x(3);

% calculate fibre direction from the specified angles
fibredir = [sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];

tau=protocol.tau;
G=protocol.G;
[M totsizeG]=size(G);

% Calculating the parallel & perpendicular signal is similar  (parallel plane & cylindrical restriction )

Bval_par = zeros(M,1);
GAMMA = 2.675987E8;

if isfield(protocol,'smalldel') && isfield(protocol,'delta')
     K = floor((protocol.smalldel+1E-10)./tau)+1;
     dstmp = floor(protocol.delta./protocol.tau)-floor((protocol.smalldel+1E-10)./protocol.tau);
     
     for m = 1:M
     G_dot_fibre = G(m,1:3:end)*fibredir(1)+ G(m,2:3:end)*fibredir(2) + G(m,3:3:end)*fibredir(3); 
    
     
       for k=1:(2*K(m)+dstmp(m))
            Fpar=sum(G_dot_fibre(1:k)*tau);            
            Bval_par(m)=Bval_par(m)+(Fpar.^2)*tau;
       end
      end
      Bval_par=GAMMA^2*Bval_par;     
      ePar=exp(-Bval_par.*dPar);    
      E = ePar;
else
     K=totsizeG/3;
     for m = 1:M
     G_dot_fibre = G(m,1:3:end)*fibredir(1)+ G(m,2:3:end)*fibredir(2) + G(m,3:3:end)*fibredir(3); 
    
       for k=1:K
            Fpar=sum(G_dot_fibre(1:k)*tau);          
            Bval_par(m)=Bval_par(m)+(Fpar.^2)*tau;
            
       end
      end
      Bval_par=GAMMA^2.*Bval_par;
      ePar=exp(-Bval_par.*dPar);
      E = ePar;
end
     
     
   
  


% Compute the Jacobian matrix
if(nargout>1)
   dx = 0.00001;
   J = zeros(length(E), 3);
    if nargin < 3 
         
         dEddPar = -Bval_par.*ePar;
         J(:,1) = dEddPar;         
         for i = 2:3
            xpert = x;
            xpert(i) = xpert(i)*(1+dx);
            Epert = Stick_GEN(xpert, protocol);
            dEdx = (Epert - E)/(xpert(i)*dx);
            J(:,i) = dEdx;
         end         
        
    else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
       
        for i = 1:length(x_deriv);
            if x_deriv(i) ~= 0  
                
                xpert = x;
                xpert(i) = xpert(i)*(1+dx);
                Epert = Stick_GEN(xpert, protocol);
                dEdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEdx;
            end
        end
        
       
        
    end   
   
end