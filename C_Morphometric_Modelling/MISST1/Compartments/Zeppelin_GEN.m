function [E J]=Zeppelin_GEN(x,protocol,x_deriv)
% function [E J] = Zeppelin_GEN(model, protocol, params_deriv)
% returns the diffusion signal and signal jacobian for the given protocol
% and Zeppelin tissue model
%
% Substrate: hindered diffusion
% Pulse sequence: Generalized gradient spin echo.
% Signal approximation: -
%
% x is the list of model parameters in SI units:
% x(1) is the free diffusivity of the material 
% x(2) is the hindered diffusivity 
% x(3) is the angle from the z direction of the axis of free diffusivity
% x(4) is the azymuthal angle from x direction
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


if ~isfield(protocol,'complex')
  protocol.complex='real';
end
% Model parameters
dPar = x(1);
dPerp = x(2);
theta = x(3);
phi = x(4);

GAMMA = 2.675987E8; % This is what is used throughout Wuzi.
%GAMMA = 2.6751525E8; % This is the latest best estimate of GAMMA (used in Camino)


% calculate fibre direction from the specified angles
fibredir = [sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];

tau=protocol.tau;
G=protocol.G;
[M totsizeG]=size(G);

% Calculating the parallel & perpendicular signal is similar  (parallel plane & cylindrical restriction )

Bval_par = zeros(M,1);
Bval_perp = zeros(M,1);

% x and y directions in he system of the principal axes
 v = [cos(theta)*cos(phi)^2+sin(phi)^2; -(1-cos(theta))*sin(phi)*cos(phi); -sin(theta)*cos(phi)]; % vectors for the new x and y directions; see Ozarslan 2010
 w = [ -(1-cos(theta))*sin(phi)*cos(phi); cos(theta)*sin(phi)^2+cos(phi)^2; -sin(theta)*sin(phi)];


     K=totsizeG/3;
     G_dot_fibre = G(:,1:3:end)*fibredir(1)+ G(:,2:3:end)*fibredir(2) + G(:,3:3:end)*fibredir(3); 
     G_dot_x =  G(:,1:3:end)*v(1)+ G(:,2:3:end)*v(2) +G(:,3:3:end)*v(3) ;  
     G_dot_y =  G(:,1:3:end)*w(1)+ G(:,2:3:end)*w(2) +G(:,3:3:end)*w(3) ;  

   for k=1:K
        Fpar=sum(G_dot_fibre(:,1:k)*tau,2);
        Fx = sum(G_dot_x(:,1:k)*tau,2);
        Fy = sum(G_dot_y(:,1:k)*tau,2);
        Bval_par=Bval_par+(Fpar.^2)*tau;
        Bval_perp=Bval_perp+(Fx.^2+Fy.^2)*tau;
   end

  Bval_par=GAMMA^2*Bval_par;
  Bval_perp=GAMMA^2*Bval_perp;
  ePar=exp(-Bval_par.*dPar);
  ePerp=exp(-Bval_perp.*dPerp);
  E = ePar.*ePerp;
   
  


% Compute the Jacobian matrix
if(nargout>1)
   dx = 0.00001;
   J = zeros(length(E), 4);
    if nargin < 3 
         
         dEddPar = -Bval_par.*E;
         dEddPerp = -Bval_perp.*E;
         J(:,1) = dEddPar;
         J(:,2) = dEddPerp;
         for i = 3:4
            xpert = x;
            xpert(i) = xpert(i)*(1+dx);
            Epert = Zeppelin_GEN(xpert, protocol);
            dEdx = (Epert - E)/(xpert(i)*dx);
            J(:,i) = dEdx;
         end         
        
    else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
    
        for i = 1:length(x_deriv);
            if x_deriv(i) ~= 0                  
                xpert = x;
                xpert(i) = xpert(i)*(1+dx);
                Epert = Zeppelin_GEN(xpert, protocol);
                dEdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEdx;
            end
        end
        
       
        
    end   
   
end
    