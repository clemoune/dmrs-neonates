function [E J]=Sphere_MM_SWOGSE_3D(x,protocol,x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a sphere compartment.
% 
% [E,J]=Sphere_MM_SWOGSE_3D(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate consisting of impermeable 
% spheres with a single radius and a diffusion protocol specified in the input
% Substrate: impermeable spheres with a single radius
% Diffusion pulse sequence: Square wave oscillating gradients with varying
% gradient orientation (SWOGSE_3D)
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
%       protocol.delta - pulse separation, size [1 N]
%       protocol.smalldel - pulse duration, size [1 N]
%       protocol.tau - sampling interval of the gradient waveform, 
%       required for MM, size [1 1] 
%       protocol.Gx - gradient strength in x direction, size [1 N];
%       protocol.Gy - gradient strength in y direction, size [1 N];
%       protocol.Gz - gradient strength in z direction, size [1 N];
%       protocol.omegax - gradient angular frequency in x direction, size [1 N];
%       protocol.omegay - gradient angular frequency in y direction, size [1 N];
%       protocol.omegaz - gradient angular frequency in z direction, size [1 N];
%       protocol.phix - phase of the gradient waveform in x direction, size [1 N];
%       protocol.phiy - phase of the gradient waveform in y direction, size [1 N];
%       protocol.phiz - phase of the gradient waveform in z direction, size [1 N];
%       protocol.angle = 1 - gradient in x direction only
%                      = 2 - gradient in x and y directions
%                      = 3 - gradient in x and z directions
%                      = 4 - gradient in x, y and z directions
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

if ~isfield(protocol,'GENj')
    protocol.GENj=1;
end

indj=protocol.GENj; % this is the index of the radii. 1 if only one.

tau=protocol.tau;


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
K = floor((protocol.smalldel+1E-10)./tau)+1;
dstmp = floor((protocol.delta-1E-10)./protocol.tau)-floor((protocol.smalldel+1E-10)./protocol.tau);

if protocol.angle == 1 % gradient only along x
    Gx = protocol.Gx';
    M = size(Gx,1);
    E=zeros(M,1); % cylinder restriction           
    G_mag = abs(Gx);
    cos_theta = zeros(size(Gx)); %Gz/G_mag, Gz = 0;
    sin_theta = sqrt(1-cos_theta.^2);
    sin_phi = zeros(size(cos_theta)); cos_phi_p = zeros(size(Gx));
    cos_phi_p(G_mag>0) = Gx(G_mag>0)./G_mag(G_mag>0)./sin_theta(G_mag>0);    
    G_ind=round(G_mag./protocol.gstep);   
    if protocol.mirror == 1
        for m=1:M
            
            time_vec = tau*(0:K(m)-1);
            sign_vec1 = (-1).^floor((protocol.omegax(m)*time_vec-protocol.phix(m))./pi-1E-10);
            Smat = real(Smat0) ;%+ imag(Smat90)*sin_phi(m)+1i*(imag(Smat0)*cos_theta(m)+imag(Smat90)*cos_phi_p(m));  
            Ecal=Smat'; 
            
            Amat_p=real(Amat0)+(imag(Amat90).*kDn).*sin_phi(m).*sin_theta(m)+1i*(imag(Amat0)*cos_theta(m)+imag(Amat90)*cos_phi_p(m)*sin_theta(m));
            Amat_m = Amat_p'; 
            RAmat_p = Rmat*Amat_p^G_ind(m);
            RAmat_m = Rmat*Amat_m^G_ind(m); 
            
            Ecal = Ecal*Rmat;  % first gradient point is 0
          
            for i = 2:K(m)-1
                 if sign_vec1(i) >0 
                     Ecal = Ecal*RAmat_p;  
                     
                 else
                     Ecal = Ecal*RAmat_m;                      
                 end                        
            end
             Ecal = Ecal*Rmat;  % last gradient point is 0
          Rmid=Rmat^dstmp(m);
          E(m)=Ecal*Rmid*Ecal';     
        end       
    elseif protocol.mirror == 0
        for m=1:M
                            
            time_vec = tau*(0:K(m)-1);
            sign_vec1 = (-1).^floor((protocol.omegax(m)*time_vec-protocol.phix(m))./pi-1E-10);
            Smat = real(Smat0); %+ imag(Smat90)*sin_phi(m)+1i*(imag(Smat0)*cos_theta(m)+imag(Smat90)*cos_phi_p(m));  
            Ecal=Smat'; 
            
            Amat_p=real(Amat0)+(imag(Amat90).*kDn).*sin_phi(m).*sin_theta(m)+1i*(imag(Amat0)*cos_theta(m)+imag(Amat90)*cos_phi_p(m)*sin_theta(m));
            Amat_m = Amat_p'; 
            RAmat_p = Rmat*Amat_p^G_ind(m);
            RAmat_m = Rmat*Amat_m^G_ind(m); 
            RAmatT_p =  Rmat*(Amat_p')^G_ind(m);
            RAmatT_m =  Rmat*(Amat_m')^G_ind(m);
          
            Ecal = Ecal*Rmat;  % first gradient point is 0
            ProdMat2 = Rmat; 
           
            for i=2:K(m)-1
                if sign_vec1(i) > 0  
                     Ecal = Ecal*RAmat_p;
                     ProdMat2 = ProdMat2*RAmatT_p;                                           
                 else
                     Ecal = Ecal*RAmat_m;
                     ProdMat2 = ProdMat2*RAmatT_m;
                 end
            end
            Ecal = Ecal*Rmat;  % first gradient point is 0
            ProdMat2 = ProdMat2*Rmat;            
          Rmid=Rmat^dstmp(m);
          E(m)=Ecal*Rmid*ProdMat2*Smat;
          
                  
        end       
    end
   
    if strcmp(protocol.complex,'complex')
      E=[real(mean(E,2));imag(mean(E,2))];
    elseif strcmp(protocol.complex,'real')
      E=real(mean(E,2));
    elseif strcmp(protocol.complex,'abs')
      E=abs(mean(E,2));
    end
elseif protocol.angle == 2
    Gx = protocol.Gx';
    Gy = protocol.Gy';
    M = size(Gx,1);
   
    E=zeros(M,1); % cylinder restriction           
    G_mag = sqrt(Gx.^2+Gy.^2);
    cos_theta = zeros(size(Gx)); %Gz/G_mag, Gz = 0;
    sin_theta = sqrt(1-cos_theta.^2);
    sin_phi_p = zeros(size(Gx)); sin_phi_m = zeros(size(Gx));
    cos_phi_p = zeros(size(Gx)); cos_phi_m = zeros(size(Gx));
    sin_phi_p(G_mag>0)  = Gy(G_mag>0) ./G_mag(G_mag>0) ./sin_theta(G_mag>0) ;
    sin_phi_m(G_mag>0)  = -Gy(G_mag>0) ./G_mag(G_mag>0) ./sin_theta(G_mag>0) ;
    cos_phi_p(G_mag>0)  = Gx(G_mag>0) ./G_mag(G_mag>0) ./sin_theta(G_mag>0) ;
    cos_phi_m(G_mag>0)  = -Gx(G_mag>0) ./G_mag(G_mag>0) ./sin_theta(G_mag>0) ;
    G_ind=round(G_mag./protocol.gstep);   
    if protocol.mirror == 1
        for m=1:M
                     
            time_vec = tau*(0:K(m)-1);
            sign_vec1 = (-1).^floor((protocol.omegax(m)*time_vec-protocol.phix(m))./pi-1E-10);
            sign_vec2 = (-1).^floor((protocol.omegay(m)*time_vec-protocol.phiy(m))./pi-1E-10);
            Smat = real(Smat0) ; %+ imag(Smat90)*sin_phi_p(m)+1i*(imag(Smat0)*cos_theta(m)+imag(Smat90)*cos_phi_p(m));  %first gradient point is 0;
            Ecal=Smat'; 
            
            Amat_pp=real(Amat0)+(imag(Amat90).*kDn).*sin_phi_p(m).*sin_theta(m)+1i*(imag(Amat0)*cos_theta(m)+imag(Amat90)*cos_phi_p(m)*sin_theta(m));
            Amat_pm=real(Amat0)+(imag(Amat90).*kDn).*sin_phi_m(m).*sin_theta(m)+1i*(imag(Amat0)*cos_theta(m)+imag(Amat90)*cos_phi_p(m)*sin_theta(m));
            Amat_mp=real(Amat0)+(imag(Amat90).*kDn).*sin_phi_p(m).*sin_theta(m)+1i*(imag(Amat0)*cos_theta(m)+imag(Amat90)*cos_phi_m(m)*sin_theta(m));
            Amat_mm=real(Amat0)+(imag(Amat90).*kDn).*sin_phi_m(m).*sin_theta(m)+1i*(imag(Amat0)*cos_theta(m)+imag(Amat90)*cos_phi_m(m)*sin_theta(m));
            
            RAmat_pp = Rmat*Amat_pp^G_ind(m);
            RAmat_pm = Rmat*Amat_pm^G_ind(m); 
            RAmat_mp = Rmat*Amat_mp^G_ind(m);
            RAmat_mm = Rmat*Amat_mm^G_ind(m);
            
           

            Ecal = Ecal*Rmat;  % first gradient point is 0
           
            for i = 2:K(m)-1
                 if sign_vec1(i)>0 && sign_vec2(i)>0
                    Ecal = Ecal*RAmat_pp;
                 elseif sign_vec1(i)>0 && sign_vec2(i)<0 
                     Ecal = Ecal*RAmat_pm;
                 elseif sign_vec1(i)<0 && sign_vec2(i)>0
                     Ecal = Ecal*RAmat_mp;
                 else
                     Ecal = Ecal*RAmat_mm;
                end               
            end
            Ecal = Ecal*Rmat;
          Rmid=Rmat^dstmp(m);
          E(m)=Ecal*Rmid*Ecal';           
       
        end      
    elseif protocol.mirror == 0
        for m=1:M


             time_vec = tau*(0:K(m)-1);
            sign_vec1 = (-1).^floor((protocol.omegax(m)*time_vec-protocol.phix(m))./pi-1E-10);
            sign_vec2 = (-1).^floor((protocol.omegay(m)*time_vec-protocol.phiy(m))./pi-1E-10);
            Smat = real(Smat0) ;%+ imag(Smat90)*sin_phi_p(m)+1i*(imag(Smat0)*cos_theta(m)+imag(Smat90)*cos_phi_p(m));  
            Ecal=Smat'; 
            
            Amat_pp=real(Amat0)+(imag(Amat90).*kDn).*sin_phi_p(m).*sin_theta(m)+1i*(imag(Amat0)*cos_theta(m)+imag(Amat90)*cos_phi_p(m)*sin_theta(m));
            Amat_pm=real(Amat0)+(imag(Amat90).*kDn).*sin_phi_m(m).*sin_theta(m)+1i*(imag(Amat0)*cos_theta(m)+imag(Amat90)*cos_phi_p(m)*sin_theta(m));
            Amat_mp=real(Amat0)+(imag(Amat90).*kDn).*sin_phi_p(m).*sin_theta(m)+1i*(imag(Amat0)*cos_theta(m)+imag(Amat90)*cos_phi_m(m)*sin_theta(m));
            Amat_mm=real(Amat0)+(imag(Amat90).*kDn).*sin_phi_m(m).*sin_theta(m)+1i*(imag(Amat0)*cos_theta(m)+imag(Amat90)*cos_phi_m(m)*sin_theta(m));
            
            RAmat_pp = Rmat*Amat_pp^G_ind(m);
            RAmat_pm = Rmat*Amat_pm^G_ind(m); 
            RAmat_mp = Rmat*Amat_mp^G_ind(m);
            RAmat_mm = Rmat*Amat_mm^G_ind(m);
            RAmatT_pp = Rmat*(Amat_pp')^G_ind(m);
            RAmatT_pm = Rmat*(Amat_pm')^G_ind(m);
            RAmatT_mp = Rmat*(Amat_mp')^G_ind(m);
            RAmatT_mm = Rmat*(Amat_mm')^G_ind(m);

           

            Ecal = Ecal*Rmat;  % first gradient point is 0
            ProdMat2 = Rmat; 
          
            for i=2:K(m)-1                  
                 if sign_vec1(i)>0 && sign_vec2(i)>0
                    Ecal = Ecal*RAmat_pp;
                    ProdMat2 = ProdMat2*RAmatT_pp;
                 elseif sign_vec1(i)>0 && sign_vec2(i)<0 
                     Ecal = Ecal*RAmat_pm;
                     ProdMat2 = ProdMat2*RAmatT_pm;
                 elseif sign_vec1(i)<0 && sign_vec2(i)>0
                     Ecal = Ecal*RAmat_mp;
                     ProdMat2 = ProdMat2*RAmatT_mp;
                 else
                     Ecal = Ecal*RAmat_mm;
                     ProdMat2 = ProdMat2*RAmatT_mm;
                end                               
            end
            Ecal = Ecal*Rmat;  % last gradient point is 0
            ProdMat2 = ProdMat2*Rmat; 
          Rmid=Rmat^dstmp(m);
          E(m)=Ecal*Rmid*ProdMat2*Smat;                        
        end
       
    end

    if strcmp(protocol.complex,'complex')
      E=[real(mean(E,2));imag(mean(E,2))];
    elseif strcmp(protocol.complex,'real')
      E=real(mean(E,2));
    elseif strcmp(protocol.complex,'abs')
      E=abs(mean(E,2));
    end
elseif protocol.angle == 3
    Gx = protocol.Gx';
    Gz = protocol.Gz';
    M = size(Gx,1);   

    E=zeros(M,1); % cylinder restriction           
    G_mag = sqrt(Gx.^2+Gz.^2);
    cos_theta_p = zeros(size(Gx)); cos_theta_m = zeros(size(Gx));
    cos_phi_p = zeros(size(Gx)); cos_phi_m = zeros(size(Gx));
    cos_theta_p(G_mag>0) = Gz(G_mag>0)./G_mag(G_mag>0); 
    cos_theta_m(G_mag>0) = -Gz(G_mag>0)./G_mag(G_mag>0);
    sin_theta = sqrt(1-cos_theta_p.^2); % the same 
    sin_phi = zeros(size(cos_theta_p));    
    cos_phi_p(G_mag>0 & sin_theta ~=0) = Gx(G_mag>0 & sin_theta ~=0)./...
        G_mag(G_mag>0 & sin_theta ~=0)./sin_theta(G_mag>0 & sin_theta ~=0);
    cos_phi_m(G_mag>0 & sin_theta ~=0) = -Gx(G_mag>0 & sin_theta ~=0)./...
        G_mag(G_mag>0 & sin_theta ~=0)./sin_theta(G_mag>0 & sin_theta ~=0);
    G_ind=round(G_mag./protocol.gstep);   
    if protocol.mirror == 1
        for m=1:M
                     
            time_vec = tau*(0:K(m)-1);
            sign_vec1 = (-1).^floor((protocol.omegax(m)*time_vec-protocol.phix(m))./pi-1E-10);
            sign_vec2 = (-1).^floor((protocol.omegaz(m)*time_vec-protocol.phiz(m))./pi-1E-10);
            Smat = real(Smat0) ;%+ imag(Smat90)*sin_phi(m)+1i*(imag(Smat0)*cos_theta_p(m)+imag(Smat90)*cos_phi_p(m));  
            Ecal=Smat'; 
            
            Amat_pp=real(Amat0)+(imag(Amat90).*kDn).*sin_phi(m).*sin_theta(m)+1i*(imag(Amat0)*cos_theta_p(m)+imag(Amat90)*cos_phi_p(m)*sin_theta(m));
            Amat_pm=real(Amat0)+(imag(Amat90).*kDn).*sin_phi(m).*sin_theta(m)+1i*(imag(Amat0)*cos_theta_m(m)+imag(Amat90)*cos_phi_p(m)*sin_theta(m));
            Amat_mp=real(Amat0)+(imag(Amat90).*kDn).*sin_phi(m).*sin_theta(m)+1i*(imag(Amat0)*cos_theta_p(m)+imag(Amat90)*cos_phi_m(m)*sin_theta(m));
            Amat_mm=real(Amat0)+(imag(Amat90).*kDn).*sin_phi(m).*sin_theta(m)+1i*(imag(Amat0)*cos_theta_m(m)+imag(Amat90)*cos_phi_m(m)*sin_theta(m));
            
            RAmat_pp = Rmat*Amat_pp^G_ind(m);
            RAmat_pm = Rmat*Amat_pm^G_ind(m); 
            RAmat_mp = Rmat*Amat_mp^G_ind(m);
            RAmat_mm = Rmat*Amat_mm^G_ind(m);            
           

            Ecal = Ecal*Rmat;  % first gradient point is 0
           
            for i = 2:K(m)-1
                 if sign_vec1(i)>0 && sign_vec2(i)>0
                    Ecal = Ecal*RAmat_pp;
                 elseif sign_vec1(i)>0 && sign_vec2(i)<0 
                     Ecal = Ecal*RAmat_pm;
                 elseif sign_vec1(i)<0 && sign_vec2(i)>0
                     Ecal = Ecal*RAmat_mp;
                 else
                     Ecal = Ecal*RAmat_mm;
                end              
            end
             Ecal = Ecal*Rmat;  % last gradient point is 0 
          Rmid=Rmat^dstmp(m);
          E(m)=Ecal*Rmid*Ecal';           
       
        end      
    elseif protocol.mirror == 0
        for m=1:M

            time_vec = tau*(0:K(m)-1);
            sign_vec1 = (-1).^floor((protocol.omegax(m)*time_vec-protocol.phix(m))./pi-1E-10);
            sign_vec2 = (-1).^floor((protocol.omegaz(m)*time_vec-protocol.phiz(m))./pi-1E-10);
            Smat = real(Smat0) ;%+ imag(Smat90)*sin_phi(m)+1i*(imag(Smat0)*cos_theta_p(m)+imag(Smat90)*cos_phi_p(m));  
            Ecal=Smat'; 
            
            Amat_pp=real(Amat0)+(imag(Amat90).*kDn).*sin_phi(m).*sin_theta(m)+1i*(imag(Amat0)*cos_theta_p(m)+imag(Amat90)*cos_phi_p(m)*sin_theta(m));
            Amat_pm=real(Amat0)+(imag(Amat90).*kDn).*sin_phi(m).*sin_theta(m)+1i*(imag(Amat0)*cos_theta_m(m)+imag(Amat90)*cos_phi_p(m)*sin_theta(m));
            Amat_mp=real(Amat0)+(imag(Amat90).*kDn).*sin_phi(m).*sin_theta(m)+1i*(imag(Amat0)*cos_theta_p(m)+imag(Amat90)*cos_phi_m(m)*sin_theta(m));
            Amat_mm=real(Amat0)+(imag(Amat90).*kDn).*sin_phi(m).*sin_theta(m)+1i*(imag(Amat0)*cos_theta_m(m)+imag(Amat90)*cos_phi_m(m)*sin_theta(m));
            
            RAmat_pp = Rmat*Amat_pp^G_ind(m);
            RAmat_pm = Rmat*Amat_pm^G_ind(m); 
            RAmat_mp = Rmat*Amat_mp^G_ind(m);
            RAmat_mm = Rmat*Amat_mm^G_ind(m);
             RAmatT_pp = Rmat*(Amat_pp')^G_ind(m);
            RAmatT_pm = Rmat*(Amat_pm')^G_ind(m);
            RAmatT_mp = Rmat*(Amat_mp')^G_ind(m);
            RAmatT_mm = Rmat*(Amat_mm')^G_ind(m);

           

            Ecal = Ecal*Rmat;  % first gradient point is 0
            ProdMat2 = Rmat; 
          
            for i=2:K(m)-1                   
                 if sign_vec1(i)>0 && sign_vec2(i)>0
                    Ecal = Ecal*RAmat_pp;
                    ProdMat2 = ProdMat2*RAmatT_pp;
                 elseif sign_vec1(i)>0 && sign_vec2(i)<0 
                     Ecal = Ecal*RAmat_pm;
                     ProdMat2 = ProdMat2*RAmatT_pm;
                 elseif sign_vec1(i)<0 && sign_vec2(i)>0
                     Ecal = Ecal*RAmat_mp;
                     ProdMat2 = ProdMat2*RAmatT_mp;
                 else
                     Ecal = Ecal*RAmat_mm;
                     ProdMat2 = ProdMat2*RAmatT_mm;
                end                    
            end
            Ecal = Ecal*Rmat;  % last gradient point is 0
            ProdMat2 = ProdMat2*Rmat;  
          Rmid=Rmat^dstmp(m);
          E(m)=Ecal*Rmid*ProdMat2*Smat;                        
        end
       
    end

    if strcmp(protocol.complex,'complex')
      E=[real(mean(E,2));imag(mean(E,2))];
    elseif strcmp(protocol.complex,'real')
      E=real(mean(E,2));
    elseif strcmp(protocol.complex,'abs')
      E=abs(mean(E,2));
    end
elseif protocol.angle == 4
     Gx = protocol.Gx';
     Gy = protocol.Gy';
     Gz = protocol.Gz';
    M = size(Gx,1);
    
    E=zeros(M,1); % cylinder restriction           
    G_mag = sqrt(Gx.^2+Gy.^2+Gz.^2);
    cos_theta_p = zeros(size(Gx)); cos_theta_m = zeros(size(Gx));
    sin_phi_p = zeros(size(Gx)); sin_phi_m = zeros(size(Gx));
    cos_phi_p = zeros(size(Gx)); cos_phi_m = zeros(size(Gx));
    
    cos_theta_p(G_mag>0 ) = Gz(G_mag~=0)./G_mag(G_mag~=0); 
    cos_theta_m(G_mag>0 ) = -Gz(G_mag~=0)./G_mag(G_mag~=0);
    sin_theta = sqrt(1-cos_theta_p.^2); % the same 
    sin_phi_p(G_mag>0 & sin_theta ~=0) = Gy(G_mag>0 & sin_theta ~=0)./...
        G_mag(G_mag>0 & sin_theta ~=0)./sin_theta(G_mag>0 & sin_theta ~=0);  
    sin_phi_m(G_mag>0 & sin_theta ~=0) = -Gy(G_mag>0 & sin_theta ~=0)./...
        G_mag(G_mag>0 & sin_theta ~=0)./sin_theta(G_mag>0 & sin_theta ~=0);    
    cos_phi_p(G_mag>0 & sin_theta ~=0) = Gx(G_mag>0 & sin_theta ~=0)./...
        G_mag(G_mag>0 & sin_theta ~=0)./sin_theta(G_mag>0 & sin_theta ~=0);
    cos_phi_m(G_mag>0 & sin_theta ~=0) = -Gx(G_mag>0 & sin_theta ~=0)./...
        G_mag(G_mag>0 & sin_theta ~=0)./sin_theta(G_mag>0 & sin_theta ~=0);
    G_ind=round(G_mag./protocol.gstep);   
    if protocol.mirror == 1
        for m=1:M
                     
            time_vec = tau*(0:K(m)-1);
            sign_vec1 = (-1).^floor((protocol.omegax(m)*time_vec-protocol.phix(m))./pi-1E-10);
            sign_vec2 = (-1).^floor((protocol.omegay(m)*time_vec-protocol.phiy(m))./pi-1E-10);
            sign_vec3 = (-1).^floor((protocol.omegaz(m)*time_vec-protocol.phiz(m))./pi-1E-10);
            Smat = real(Smat0) ;%+ imag(Smat90)*sin_phi_p(m)+1i*(imag(Smat0)*cos_theta_p(m)+imag(Smat90)*cos_phi_p(m));  
            Ecal=Smat'; 
            
            Amat_ppp=real(Amat0)+(imag(Amat90).*kDn).*sin_phi_p(m).*sin_theta(m)+1i*(imag(Amat0)*cos_theta_p(m)+imag(Amat90)*cos_phi_p(m)*sin_theta(m));
            Amat_ppm=real(Amat0)+(imag(Amat90).*kDn).*sin_phi_p(m).*sin_theta(m)+1i*(imag(Amat0)*cos_theta_m(m)+imag(Amat90)*cos_phi_p(m)*sin_theta(m));
            Amat_pmp=real(Amat0)+(imag(Amat90).*kDn).*sin_phi_m(m).*sin_theta(m)+1i*(imag(Amat0)*cos_theta_p(m)+imag(Amat90)*cos_phi_p(m)*sin_theta(m));
            Amat_pmm=real(Amat0)+(imag(Amat90).*kDn).*sin_phi_m(m).*sin_theta(m)+1i*(imag(Amat0)*cos_theta_m(m)+imag(Amat90)*cos_phi_p(m)*sin_theta(m));
            Amat_mpp=real(Amat0)+(imag(Amat90).*kDn).*sin_phi_p(m).*sin_theta(m)+1i*(imag(Amat0)*cos_theta_p(m)+imag(Amat90)*cos_phi_m(m)*sin_theta(m));
            Amat_mpm=real(Amat0)+(imag(Amat90).*kDn).*sin_phi_p(m).*sin_theta(m)+1i*(imag(Amat0)*cos_theta_m(m)+imag(Amat90)*cos_phi_m(m)*sin_theta(m));
            Amat_mmp=real(Amat0)+(imag(Amat90).*kDn).*sin_phi_m(m).*sin_theta(m)+1i*(imag(Amat0)*cos_theta_p(m)+imag(Amat90)*cos_phi_m(m)*sin_theta(m));
            Amat_mmm=real(Amat0)+(imag(Amat90).*kDn).*sin_phi_m(m).*sin_theta(m)+1i*(imag(Amat0)*cos_theta_m(m)+imag(Amat90)*cos_phi_m(m)*sin_theta(m));
            
            RAmat_ppp = Rmat*Amat_ppp^G_ind(m);
            RAmat_ppm = Rmat*Amat_ppm^G_ind(m); 
            RAmat_pmp = Rmat*Amat_pmp^G_ind(m);
            RAmat_pmm = Rmat*Amat_pmm^G_ind(m);
            RAmat_mpp = Rmat*Amat_mpp^G_ind(m);
            RAmat_mpm = Rmat*Amat_mpm^G_ind(m); 
            RAmat_mmp = Rmat*Amat_mmp^G_ind(m);
            RAmat_mmm = Rmat*Amat_mmm^G_ind(m);
            
           

            Ecal = Ecal*Rmat;  % first gradient point is 0
           
             for i = 2:K(m)-1
                 if sign_vec1(i) >0  && sign_vec2(i) >0  && sign_vec3(i) >0 
                    Ecal = Ecal*RAmat_ppp;

                 elseif sign_vec1(i) >0  && sign_vec2(i) >0 && sign_vec3(i) <0 
                     Ecal = Ecal*RAmat_ppm;

                 elseif sign_vec1(i) >0  && sign_vec2(i) <0 && sign_vec3(i) >0 
                     Ecal = Ecal*RAmat_pmp;

                 elseif sign_vec1(i) >0  && sign_vec2(i) <0 && sign_vec3(i) <0
                     Ecal = Ecal*RAmat_pmm;

                 elseif sign_vec1(i) <0  && sign_vec2(i) >0 && sign_vec3(i) >0
                    Ecal = Ecal*RAmat_mpp;

                 elseif sign_vec1(i) <0  && sign_vec2(i) >0 && sign_vec3(i) <0
                     Ecal = Ecal*RAmat_mpm;

                 elseif sign_vec1(i) <0  && sign_vec2(i) <0 && sign_vec3(i) >0
                     Ecal = Ecal*RAmat_mmp;

                 else 
                     Ecal = Ecal*RAmat_mmm;  

                 end                     
             end
              Ecal = Ecal*Rmat;  % last gradient point is 0
          Rmid=Rmat^dstmp(m);
          E(m)=Ecal*Rmid*Ecal';           
       
        end      
    elseif protocol.mirror == 0
        for m=1:M


            time_vec = tau*(0:K(m)-1);
            sign_vec1 = (-1).^floor((protocol.omegax(m)*time_vec-protocol.phix(m))./pi-1E-10);
            sign_vec2 = (-1).^floor((protocol.omegay(m)*time_vec-protocol.phiy(m))./pi-1E-10);
            sign_vec3 = (-1).^floor((protocol.omegaz(m)*time_vec-protocol.phiz(m))./pi-1E-10);
            Smat = real(Smat0) ;% + imag(Smat90)*sin_phi_p(m)+1i*(imag(Smat0)*cos_theta_p(m)+imag(Smat90)*cos_phi_p(m));  
            Ecal=Smat'; 
            
            Amat_ppp=real(Amat0)+(imag(Amat90).*kDn).*sin_phi_p(m).*sin_theta(m)+1i*(imag(Amat0)*cos_theta_p(m)+imag(Amat90)*cos_phi_p(m)*sin_theta(m));
            Amat_ppm=real(Amat0)+(imag(Amat90).*kDn).*sin_phi_p(m).*sin_theta(m)+1i*(imag(Amat0)*cos_theta_m(m)+imag(Amat90)*cos_phi_p(m)*sin_theta(m));
            Amat_pmp=real(Amat0)+(imag(Amat90).*kDn).*sin_phi_m(m).*sin_theta(m)+1i*(imag(Amat0)*cos_theta_p(m)+imag(Amat90)*cos_phi_p(m)*sin_theta(m));
            Amat_pmm=real(Amat0)+(imag(Amat90).*kDn).*sin_phi_m(m).*sin_theta(m)+1i*(imag(Amat0)*cos_theta_m(m)+imag(Amat90)*cos_phi_p(m)*sin_theta(m));
            Amat_mpp=real(Amat0)+(imag(Amat90).*kDn).*sin_phi_p(m).*sin_theta(m)+1i*(imag(Amat0)*cos_theta_p(m)+imag(Amat90)*cos_phi_m(m)*sin_theta(m));
            Amat_mpm=real(Amat0)+(imag(Amat90).*kDn).*sin_phi_p(m).*sin_theta(m)+1i*(imag(Amat0)*cos_theta_m(m)+imag(Amat90)*cos_phi_m(m)*sin_theta(m));
            Amat_mmp=real(Amat0)+(imag(Amat90).*kDn).*sin_phi_m(m).*sin_theta(m)+1i*(imag(Amat0)*cos_theta_p(m)+imag(Amat90)*cos_phi_m(m)*sin_theta(m));
            Amat_mmm=real(Amat0)+(imag(Amat90).*kDn).*sin_phi_m(m).*sin_theta(m)+1i*(imag(Amat0)*cos_theta_m(m)+imag(Amat90)*cos_phi_m(m)*sin_theta(m));
            
            RAmat_ppp = Rmat*Amat_ppp^G_ind(m);
            RAmat_ppm = Rmat*Amat_ppm^G_ind(m); 
            RAmat_pmp = Rmat*Amat_pmp^G_ind(m);
            RAmat_pmm = Rmat*Amat_pmm^G_ind(m);
            RAmat_mpp = Rmat*Amat_mpp^G_ind(m);
            RAmat_mpm = Rmat*Amat_mpm^G_ind(m); 
            RAmat_mmp = Rmat*Amat_mmp^G_ind(m);
            RAmat_mmm = Rmat*Amat_mmm^G_ind(m);
            RAmatT_ppp = Rmat*(Amat_ppp')^G_ind(m);
            RAmatT_ppm = Rmat*(Amat_ppm')^G_ind(m);
            RAmatT_pmp = Rmat*(Amat_pmp')^G_ind(m);
            RAmatT_pmm = Rmat*(Amat_pmm')^G_ind(m);
            RAmatT_mpp = Rmat*(Amat_mpp')^G_ind(m);
            RAmatT_mpm = Rmat*(Amat_mpm')^G_ind(m);
            RAmatT_mmp = Rmat*(Amat_mmp')^G_ind(m);
            RAmatT_mmm = Rmat*(Amat_mmm')^G_ind(m);           

           

            Ecal = Ecal*Rmat;  % first gradient point is 0
            ProdMat2 = Rmat; 
                     
            for i = 2:K(m)-1
                if sign_vec1(i) >0  && sign_vec2(i) >0  && sign_vec3(i) >0 
                    Ecal = Ecal*RAmat_ppp;
                    ProdMat2 = ProdMat2*RAmatT_ppp;
                 elseif sign_vec1(i) >0  && sign_vec2(i) >0 && sign_vec3(i) <0 
                     Ecal = Ecal*RAmat_ppm;
                      ProdMat2 = ProdMat2*RAmatT_ppm; 
                 elseif sign_vec1(i) >0  && sign_vec2(i) <0 && sign_vec3(i) >0 
                     Ecal = Ecal*RAmat_pmp;
                     ProdMat2 = ProdMat2*RAmatT_pmp;
                 elseif sign_vec1(i) >0  && sign_vec2(i) <0 && sign_vec3(i) <0
                     Ecal = Ecal*RAmat_pmm;
                     ProdMat2 = ProdMat2*RAmatT_pmm;
                 elseif sign_vec1(i) <0  && sign_vec2(i) >0 && sign_vec3(i) >0
                    Ecal = Ecal*RAmat_mpp;
                    ProdMat2 = ProdMat2*RAmatT_mpp;
                 elseif sign_vec1(i) <0  && sign_vec2(i) >0 && sign_vec3(i) <0
                     Ecal = Ecal*RAmat_mpm;
                     ProdMat2 = ProdMat2*RAmatT_mpm;
                 elseif sign_vec1(i) <0  && sign_vec2(i) <0 && sign_vec3(i) >0
                     Ecal = Ecal*RAmat_mmp;
                     ProdMat2 = ProdMat2*RAmatT_mmp;
                 else 
                     Ecal = Ecal*RAmat_mmm;  
                     ProdMat2 = ProdMat2*RAmatT_mmm;
                 end                                     
            end
           Ecal = Ecal*Rmat;  % last gradient point is 0
           ProdMat2 = ProdMat2*Rmat; 
          Rmid=Rmat^dstmp(m);
          E(m)=Ecal*Rmid*ProdMat2*Smat;                        
        end
       
    end

    if strcmp(protocol.complex,'complex')
      E=[real(mean(E,2));imag(mean(E,2))];
    elseif strcmp(protocol.complex,'real')
      E=real(mean(E,2));
    elseif strcmp(protocol.complex,'abs')
      E=abs(mean(E,2));
    end
end




% Compute the Jacobian matrix
if(nargout>1)
   
    dx = protocol.pert;
    J = zeros(length(E),2);
    if nargin < 3 
         
        protocol.diff=1;   
        xpert = x;
        xpert(1) = xpert(1)*(1+dx);    
        Epert = Sphere_MM_SWOGSE_3D(xpert,protocol);
        dEtdD = (Epert - E)/(xpert(1)*dx);

        protocol.diff=2;
        xpert = x;
        xpert(2) = xpert(2)*(1+dx);    
        Epert = Sphere_MM_SWOGSE_3D(xpert,protocol);
        dEtda = (Epert - E)/(xpert(2)*dx);
        
        J(:,1) = dEtdD;
        J(:,2) = dEtda;
        
    else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
        
        
            if x_deriv(1) ~= 0  
              
                 xpert = x;
                protocol.diff=1;   
                xpert(1) = xpert(1)*(1+dx);    
                Epert = Sphere_MM_SWOGSE_3D(xpert,protocol,x_deriv);
                dEtdx = (Epert - E)/(xpert(1)*dx);
                J(:,1) = dEtdx;
            elseif  x_deriv(2) ~= 0  
                
                 xpert = x;
                protocol.diff=2;   
                xpert(2) = xpert(2)*(1+dx);    
                Epert = Sphere_MM_SWOGSE_3D(xpert,protocol,x_deriv);
                dEtdx = (Epert - E)/(xpert(2)*dx);
                J(:,2) = dEtdx;
                
            end             

    end  
    protocol.diff=0;
  
   
end