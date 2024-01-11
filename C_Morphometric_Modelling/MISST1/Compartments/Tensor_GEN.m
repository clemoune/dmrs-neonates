function [E J]=Tensor_GEN(x,protocol,x_deriv)
% function [E J] = Tensor_GEN(model, protocol, params_deriv)
% returns the diffusion signal and signal jacobian for the given protocol
% and Tensor tissue model
%
% Substrate: hindered diffusion
% Pulse sequence: Generalized gradient spin echo.
% Signal approximation: -
%
% x is the list of model parameters in SI units:
% x(1) is the diffusivity of the material along the first direction
% x(2) is the diffusivity of the material along the second direction
% x(3) is the diffusivity of the material along the third direction
% x(4) is the angle from the z direction of the main direction of the
% tensor
% x(5) is the azymuthal angle of the main direction of the
% tensor
% x(6) is the angle of rotation of the second direction 
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
d1 = x(1);
d2 = x(2);
d3 = x(3);
theta = x(4);
phi = x(5);
psi = x(6);


GAMMA = 2.675987E8; % This is what is used throughout Wuzi.
%GAMMA = 2.6751525E8; % This is the latest best estimate of GAMMA (used in Camino)


% calculate main direction from the specified angles
n1 = [sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];

% rotate by 90 degrees in xz plane
n2_zrot = [sin(theta+pi/2)*cos(phi); sin(theta+pi/2)*sin(phi); cos(theta+pi/2)];

n2 = RotateVec(n2_zrot,n1,psi);

n3 = cross(n1,n2);

tau=protocol.tau;
G=protocol.G;

% dot product between gradient direction and main direction of the tensor
F1 = cumsum((G(:,1:3:end).*n1(1)+G(:,2:3:end).*n1(2)+G(:,3:3:end).*n1(3)).*tau,2); 
% dot product between gradient direction and second direction of the tensor
F2 = cumsum((G(:,1:3:end).*n2(1)+G(:,2:3:end).*n2(2)+G(:,3:3:end).*n2(3)).*tau,2); 
% dot product between gradient direction and third direction of the tensor
F3 = cumsum((G(:,1:3:end).*n3(1)+G(:,2:3:end).*n3(2)+G(:,3:3:end).*n3(3)).*tau,2);

B1=GAMMA^2*sum((F1.^2)*tau,2);   
B2=GAMMA^2*sum((F2.^2)*tau,2); 
B3=GAMMA^2*sum((F3.^2)*tau,2); 


E=exp(-(B1.*d1 + B2*d2 +B3*d3) );


% Compute the Jacobian matrix
if(nargout>1)
   dx = 0.00001;
   J = zeros(length(E), 4);
    if nargin < 3 
         
         dEdd1 = -B1.*E;
         dEdd2 = -B2.*E;
         dEdd3 = -B3.*E;
         J(:,1) = dEdd1;
         J(:,2) = dEdd2;
         J(:,3) = dEdd3;
         for i = 4:6
            xpert = x;
            xpert(i) = xpert(i)*(1+dx);
            Epert = Tensor_GEN(xpert, protocol);
            dEdx = (Epert - E)/(xpert(i)*dx);
            J(:,i) = dEdx;
         end         
        
    else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
    
        for i = 1:length(x_deriv);
            if x_deriv(i) ~= 0                  
                xpert = x;
                xpert(i) = xpert(i)*(1+dx);
                Epert = Tensor_GEN(xpert, protocol);
                dEdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEdx;
            end
        end
        
       
        
    end   
   
end
end
function vec2 = RotateVec(vec1,axis,alpha) % rotates vector vec1 around axis by alpha
    axis = axis./norm(axis);
    % Rodriguez formula
    vec2 = cos(alpha)*vec1+sin(alpha)*cross(vec1,axis) + (1-cos(alpha))*dot(vec1,axis)*axis;
   
end