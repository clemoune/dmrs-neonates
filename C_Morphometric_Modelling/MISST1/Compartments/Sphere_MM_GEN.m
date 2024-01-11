function [E J]=Sphere_MM_GEN(x,protocol,x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a sphere compartment.
% 
% [E,J]=Sphere_MM_GEN(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate consisting of impermeable 
% spheres with a single radius and a diffusion protocol specified in the input
% Substrate: impermeable spheres with a single radius
% Diffusion pulse sequence: Generalized waveform (GEN)
% Signal approximation: Matrix method (MM) - propagator expressed via 
%       eigenmode expansion
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% x - size 2 vector of model parameters in SI units for Sphere:
%       x(1) - free diffusivity of the material inside the sphere.
%       x(2) - sphere radius   
% protocol - structure which includes all the information related to the 
%        diffusion protocol and required to generate the signal. 
%       protocol.G - the discretized gradient waveform, size [N, 3*K], 
%       where N is the number of measurements and K is the number of
%       gradient points in each direction.
%       [G11x G11y G11z G12x ... G1Kx G1Ky G1Kz]
%       ......
%       [GN1x GN1y GN1z GN2x ... GNKx GNKy GNKz]
%       protocol.tau - sampling interval of the gradient waveform, required 
%       for MM, size [1 1] 
%       protocol.cube_rotation - cell array which stores rotation matrices
%       describing the orientation of the cuboids. If it has more than one
%       matrix, than the average is computed; The default is the identity
%       matrix.
%       optional: if the total gradient waveform is formed by two repeated 
%       or reflected waveforms, specifying the duration of each waveform 
%       (protocol.smalldel, size [1 N]) and the separtation of the waveforms
%       (protocol.delta, size [1 N]) saves computation time.
%       all other fields are assigned by the MMConstants.m function 
% x_deriv - a vector of 0s and 1s with the same size as x, indicating
%       which parameters are considered in the Jacobian;
%
%------------------------------------------------------------------------
% This file is part of the camino.m toolbox.
% Copyright (c) 2015, UCL Microstructure Imaging Group (MIG), All rights reserved.
% Distributed under the Modified BSD Licence (see: LICENSE.pdf).
%
% Authors:
%   Andrada Ianus (a.ianus.11@ucl.ac.uk)
%   Ivana Drobnjak (i.drobnjak@ucl.ac.uk)
%   Daniel C. Alexander (d.alexander@ucl.ac.uk)
%
% 

if ~isfield(protocol,'GENj')
    protocol.GENj=1;
end

indj=protocol.GENj; % this is the index of the radii. 1 if only one.

tau=protocol.tau;
G=protocol.G;
[M totsizeG]=size(G);

if protocol.diff==0
    Amat0=protocol.A_sph{indj,1};
    Smat0=protocol.S_sph{indj,1};
    Rmat=protocol.R_sph{indj};
    Amat90=protocol.A_sph{indj,2};
    Smat90=protocol.S_sph{indj,2};
    
elseif protocol.diff==1 % matrices for cumputing derivatives with respect to diffusivity
    Amat0=protocol.A_sph{indj,1};
    Smat0=protocol.S_sph{indj,1};
    Rmat=protocol.RpertD_sph{indj};
    Amat90=protocol.A_sph{indj,2};
    Smat90=protocol.S_sph{indj,2};
    
elseif protocol.diff==2 % matrices for cumputing derivatives with respect to radius
    Amat0=protocol.Aperta_sph{indj,1};
    Smat0=protocol.Sperta_sph{indj,1};
    Rmat=protocol.Rperta_sph{indj};
    Amat90=protocol.Aperta_sph{indj,2};
    Smat90=protocol.Sperta_sph{indj,2};
    
else
    error('protocol.diff is not suitable')
end
kDn=protocol.kDn_sph;

E=zeros(M,1);


if isfield(protocol,'smalldel') && isfield(protocol,'delta')
    K = floor((protocol.smalldel+1E-10)./tau)+1;
   dstmp = floor((protocol.delta-1E-10)./protocol.tau)-floor((protocol.smalldel+1E-10)./protocol.tau);
   
    
     %disp('Angle method using exp(i(k-n)theta)')
    G_mag=sqrt(G(:,1:3:end).^2+G(:,2:3:end).^2+G(:,3:3:end).^2);
    Gx = G(:,1:3:end);
    Gy = G(:,2:3:end);
    Gz = G(:,3:3:end);
     cos_theta = zeros(size(Gx));
    sin_phi = zeros(size(cos_theta));
    cos_phi = zeros(size(cos_theta));
    cos_theta(G_mag >0) = Gz(G_mag >0)./G_mag(G_mag >0);
    sin_theta = sqrt(1-cos_theta.^2);    
    
    sin_phi(sin_theta ~= 0 & G_mag >0) = Gy(sin_theta~=0 & G_mag >0)./G_mag(sin_theta~=0 & G_mag >0)./sin_theta(sin_theta~=0 & G_mag >0);
    cos_phi(sin_theta ~= 0 & G_mag >0) = Gx(sin_theta~=0 & G_mag >0)./G_mag(sin_theta~=0 & G_mag >0)./sin_theta(sin_theta~=0 & G_mag >0);  
    
    Gind_mat=round(G_mag./protocol.gstep);   
    if protocol.mirror == 1
        for m=1:M
                  Smat = real(Smat0) + imag(Smat90)*sin_phi(m,1)+1i*(imag(Smat0)*cos_theta(m,1)+imag(Smat90)*cos_phi(m,1));  
                  Ecal=Smat';
                  for ind=1:K(m)
                    Gind=Gind_mat(m,ind);
                    if Gind~=0
                      AmatU=real(Amat0)+(imag(Amat90).*kDn).*sin_phi(m,ind).*sin_theta(m,ind)+1i*(imag(Amat0)*cos_theta(m,ind)+imag(Amat90)*cos_phi(m,ind)*sin_theta(m,ind));
                      AmatUD=AmatU^Gind;
                      Ecal=Ecal*Rmat*AmatUD;
                    elseif Gind==0
                      Ecal=Ecal*Rmat;
                    end
                  end
                  Rmid=Rmat^dstmp(m);
                  E(m)=Ecal*Rmid*Ecal';                  
         end
    elseif protocol.mirror == 0
        for m=1:M
                  Smat = real(Smat0) + imag(Smat90)*sin_phi(m,1)+1i*(imag(Smat0)*cos_theta(m,1)+imag(Smat90)*cos_phi(m,1));    
                  Ecal=Smat';
                  ProdMat2 = eye(size(Rmat));
                  for ind=1:K(m)
                    Gind=Gind_mat(m,ind);
                    if Gind~=0
                    AmatU=real(Amat0)+(imag(Amat90).*kDn).*sin_phi(m,ind).*sin_theta(m,ind)+1i*(imag(Amat0)*cos_theta(m,ind)+imag(Amat90)*cos_phi(m,ind)*sin_theta(m,ind));
                      AmatUD=AmatU^Gind;
                      AmatUDT = AmatUD';
                      ProdMat2 = ProdMat2*Rmat*AmatUDT; % the matrix product for the second gradient
                      Ecal=Ecal*Rmat*AmatUD;
                    elseif Gind==0
                      ProdMat2 = ProdMat2*Rmat;  
                      Ecal=Ecal*Rmat;
                    end
                  end
                  Rmid=Rmat^dstmp(m);
                  E(m)=Ecal*Rmid*ProdMat2*Smat;                 
        end     
        
    end
else
    K=totsizeG/3;
    % Perpendicular
  

    %disp('Angle method using exp(i(k-n)theta)')
    G_mag=sqrt(G(:,1:3:end).^2+G(:,2:3:end).^2+G(:,3:3:end).^2);
    Gx = G(:,1:3:end);
    Gy = G(:,2:3:end); 
    Gz = G(:,3:3:end); 
    
    cos_theta = zeros(size(Gx));
    sin_phi = zeros(size(cos_theta));
    cos_phi = zeros(size(cos_theta));
    cos_theta(G_mag >0) = Gz(G_mag >0)./G_mag(G_mag >0);
    sin_theta = sqrt(1-cos_theta.^2);    
    
    sin_phi(sin_theta ~= 0 & G_mag >0) = Gy(sin_theta~=0 & G_mag >0)./G_mag(sin_theta~=0 & G_mag >0)./sin_theta(sin_theta~=0 & G_mag >0);
    cos_phi(sin_theta ~= 0 & G_mag >0) = Gx(sin_theta~=0 & G_mag >0)./G_mag(sin_theta~=0 & G_mag >0)./sin_theta(sin_theta~=0 & G_mag >0);  
     
    Gind_mat=round(G_mag./protocol.gstep);
    for m=1:M
              Smat = real(Smat0) + imag(Smat90)*sin_phi(m,1)+1i*(imag(Smat0)*cos_theta(m,1)+imag(Smat90)*cos_phi(m,1));      
              Ecal=Smat';
              for ind=1:K
                Gind=Gind_mat(m,ind);
                if Gind~=0
                 AmatU=real(Amat0)+(imag(Amat90).*kDn).*sin_phi(m,ind).*sin_theta(m,ind)+1i*(imag(Amat0)*cos_theta(m,ind)+imag(Amat90)*cos_phi(m,ind)*sin_theta(m,ind));
                  AmatUD=AmatU^Gind;
                  Ecal=Ecal*Rmat*AmatUD;
                elseif Gind==0
                  Ecal=Ecal*Rmat;
                end
              end
              E(m)=Ecal*Rmat*Smat;
    end
end
    if strcmp(protocol.complex,'complex')
      E=[real(mean(E,2));imag(mean(E,2))];
    elseif strcmp(protocol.complex,'real')
      E=real(mean(E,2));
    elseif strcmp(protocol.complex,'abs')
      E=abs(mean(E,2));
    end

% Compute the Jacobian matrix
if(nargout>1)
    J = zeros(length(E),2);
    dx = protocol.pert;
    if nargin < 3 
         
        protocol.diff=1;   
        xpert = x;
        xpert(1) = xpert(1)*(1+dx);    
        Epert = Sphere_MM_GEN(xpert,protocol);
        dEtdD = (Epert - E)/(xpert(1)*dx);

        protocol.diff=2;
        xpert = x;
        xpert(2) = xpert(2)*(1+dx);    
        Epert = Sphere_MM_GEN(xpert,protocol);
        dEtda = (Epert - E)/(xpert(2)*dx);
        
        J(:,1) = dEtdD;
        J(:,2) = dEtda;
        
    else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
       
        
            if x_deriv(1) ~= 0  
               
                 xpert = x;
                protocol.diff=1;   
                xpert(1) = xpert(1)*(1+dx);    
                Epert = Sphere_MM_GEN(xpert,protocol,x_deriv);
                dEtdx = (Epert - E)/(xpert(1)*dx);
                J(:,1) = dEtdx;
            elseif  x_deriv(2) ~= 0  
                
                 xpert = x;
                protocol.diff=2;   
                xpert(2) = xpert(2)*(1+dx);    
                Epert = Sphere_MM_GEN(xpert,protocol,x_deriv);
                dEtdx = (Epert - E)/(xpert(2)*dx);
                J(:,2) = dEtdx;
                
            end             

    end  
    protocol.diff=0;
  
   
end