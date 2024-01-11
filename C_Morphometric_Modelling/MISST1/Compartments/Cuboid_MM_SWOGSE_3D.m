function [E J]=Cuboid_MM_SWOGSE_3D(x,protocol,x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a cuboid/ensemble of cuboids.
% 
% [E,J]=Cuboid_MM_SWOGSE_3D(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate consisting of a cuboid/ 
% an ensemble of cuboids and a diffusion protocol specified in the input
% Substrate: cuboid / ensemble of cuboids
% Diffusion pulse sequence: Square wave oscillating gradients with varying
% gradient orientation (SWOGSE_3D)
% Signal approximation: Matrix method (MM) - propagator expressed via 
%       eigenmode expansion
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% x - size 4 vector of model parameters in SI units for Cuboid:
%       x(1) - free diffusivity of the material inside the cuboids.
%       x(2) - length of the cuboid (x dir)
%       x(3) - width of the cuboid (y dir)
%       x(4) - height of the cuboid (z dir) 
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
% 
if ~isfield(protocol,'GENj')
    protocol.GENj=1;
end
indj=protocol.GENj; % this is the index of the radii. 1 if only one.
if ~isfield(protocol,'complex')
  protocol.complex='real';
end
if ~isfield(protocol,'cube_rotation') % a cell containing rotation matrices. the signal is computed for the ensemble of rotations
    cube_rotation{1} = {eye(3)};
else
    cube_rotation = protocol.cube_rotation;
end
% Model parameters
tau=protocol.tau;
RotN = length(cube_rotation);
K = floor((protocol.smalldel+1E-10)./tau)+1;
dstmp = floor((protocol.delta-1E-10)./protocol.tau)-floor((protocol.smalldel+1E-10)./protocol.tau);
% Calculating the parallel & perpendicular signal is similar  (parallel plane & cylindrical restriction )


% Perpendicular
if protocol.diff==0
    Amat_x=protocol.A_x{indj};
    Smat_x=protocol.S_x{indj};
    Rmat_x=protocol.R_x{indj};
    Amat_y=protocol.A_y{indj};
    Smat_y=protocol.S_y{indj};
    Rmat_y=protocol.R_y{indj};
    Amat_z=protocol.A_z{indj};
    Smat_z=protocol.S_z{indj};
    Rmat_z=protocol.R_z{indj};    
 
elseif protocol.diff==1
    Amat_x=protocol.A_x{indj};
    Smat_x=protocol.S_x{indj};
    Rmat_x=protocol.RpertD_x{indj};
    Amat_y=protocol.A_y{indj};
    Smat_y=protocol.S_y{indj};
    Rmat_y=protocol.RpertD_y{indj};
    Amat_z=protocol.A_z{indj};
    Smat_z=protocol.S_z{indj};
    Rmat_z=protocol.RpertD_z{indj};
elseif protocol.diff==2
    Amat_x=protocol.Aperta_x{indj};
    Smat_x=protocol.Sperta_x{indj};
    Rmat_x=protocol.Rperta_x{indj};
    Amat_y=protocol.A_y{indj};
    Smat_y=protocol.S_y{indj};
    Rmat_y=protocol.R_y{indj};
    Amat_z=protocol.A_z{indj};
    Smat_z=protocol.S_z{indj};
    Rmat_z=protocol.R_z{indj}; 
 
elseif protocol.diff==3
    Amat_x=protocol.A_x{indj};
    Smat_x=protocol.S_x{indj};
    Rmat_x=protocol.R_x{indj};
    Amat_y=protocol.Apertb_y{indj};
    Smat_y=protocol.Spertb_y{indj};
    Rmat_y=protocol.Rpertb_y{indj};
    Amat_z=protocol.A_z{indj};
    Smat_z=protocol.S_z{indj};
    Rmat_z=protocol.R_z{indj};
elseif protocol.diff == 4
    Amat_x=protocol.A_x{indj};
    Smat_x=protocol.S_x{indj};
    Rmat_x=protocol.R_x{indj};
    Amat_y=protocol.A_y{indj};
    Smat_y=protocol.S_y{indj};
    Rmat_y=protocol.R_y{indj};
    Amat_z=protocol.Apertc_z{indj};
    Smat_z=protocol.Spertc_z{indj};
    Rmat_z=protocol.Rpertc_z{indj};
else
    error('protocol.diff is not suitable')
end

if protocol.angle == 1 % gradient only along x
    Gx = protocol.Gx';
    M = size(Gx,1);
   
      % get Ndir directions on a spehere
    eX = zeros(M,RotN); % x direction
    eY = zeros(M,RotN); % y direction
    eZ = zeros(M,RotN); % z direction  
    G_mag = Gx;
    G_ind=round(G_mag./protocol.gstep);
    qhat_X=Gx./G_mag; % qhatx  
    if protocol.mirror == 1
        for m=1:M
            time_vec = tau*(0:K(m)-1);
            sign_vec1 = (-1).^floor((protocol.omegax(m)*time_vec-protocol.phix(m))./pi-1E-10);
            for n = 1:RotN
                % Calculating parallel & perpendicular signals
                inv_rot = cube_rotation{n}';                               
                  
                Ecal_x=Smat_x'; 
                Ecal_y=Smat_y'; 
                Ecal_z=Smat_z';                 

                grad_dir_p = inv_rot*[qhat_X(m); 0; 0];                  

                Amat_x_p=real(Amat_x)+ 1i* imag(Amat_x)*grad_dir_p(1);
                Amat_x_m=Amat_x_p'; 
                RAmat_x_p = Rmat_x*Amat_x_p^G_ind(m);
                RAmat_x_m = Rmat_x*Amat_x_m^G_ind(m); 
 
                Amat_y_p=real(Amat_y)+ 1i* imag(Amat_y)*grad_dir_p(2);
                Amat_y_m=Amat_y_p'; 
                RAmat_y_p = Rmat_y*Amat_y_p^G_ind(m);
                RAmat_y_m = Rmat_y*Amat_y_m^G_ind(m);

                Amat_z_p=real(Amat_z)+ 1i* imag(Amat_z)*grad_dir_p(3);
                Amat_z_m=Amat_z_p'; 
                RAmat_z_p = Rmat_z*Amat_z_p^G_ind(m);
                RAmat_z_m = Rmat_z*Amat_z_m^G_ind(m);
               
                
                Ecal_x = Ecal_x*Rmat_x;  % first gradient point is 0
                Ecal_y = Ecal_y*Rmat_y;  % first gradient point is 0
                Ecal_z = Ecal_z*Rmat_z;  % first gradient point is 0
                for i = 2:K(m)-1
                     if sign_vec1(i) >0 
                         Ecal_x = Ecal_x*RAmat_x_p;  
                         Ecal_y = Ecal_y*RAmat_y_p;   
                         Ecal_z = Ecal_z*RAmat_z_p;   
                     else
                         Ecal_x = Ecal_x*RAmat_x_m;
                         Ecal_y = Ecal_y*RAmat_y_m;
                         Ecal_z = Ecal_z*RAmat_z_m;
                     end                        
                end
                 Ecal_x = Ecal_x*Rmat_x;
                 Ecal_y = Ecal_y*Rmat_y;
                 Ecal_z = Ecal_z*Rmat_z;
              Rmid_x=Rmat_x^dstmp(m);
              eX(m,n)=Ecal_x*Rmid_x*Ecal_x';  
              Rmid_y=Rmat_y^dstmp(m);
              eY(m,n)=Ecal_y*Rmid_y*Ecal_y'; 
              Rmid_z=Rmat_z^dstmp(m);
              eZ(m,n)=Ecal_z*Rmid_z*Ecal_z'; 
                 
            end
        end       
    elseif protocol.mirror == 0
        for m=1:M
            time_vec = tau*(0:K(m)-1);
            sign_vec1 = (-1).^floor((protocol.omegax(m)*time_vec-protocol.phix(m))./pi-1E-10);
             for n = 1:RotN
                % Calculating parallel & perpendicular signals
                inv_rot = cube_rotation{n}';                    
                  
                Ecal_x=Smat_x'; 
                Ecal_y=Smat_y'; 
                Ecal_z=Smat_z'; 
                
                grad_dir_p = inv_rot*[qhat_X(m); 0; 0];                  

                Amat_x_p=real(Amat_x)+ 1i* imag(Amat_x)*grad_dir_p(1);
                Amat_x_m=Amat_x_p'; 
                RAmat_x_p = Rmat_x*Amat_x_p^G_ind(m);
                RAmat_x_m = Rmat_x*Amat_x_m^G_ind(m); 
                RAmatT_x_p = Rmat_x*(Amat_x_p')^G_ind(m);
                RAmatT_x_m = Rmat_x*(Amat_x_m')^G_ind(m);
 
                Amat_y_p=real(Amat_y)+ 1i* imag(Amat_y)*grad_dir_p(2);
                Amat_y_m=Amat_y_p'; 
                RAmat_y_p = Rmat_y*Amat_y_p^G_ind(m);
                RAmat_y_m = Rmat_y*Amat_y_m^G_ind(m);
                RAmatT_y_p = Rmat_y*(Amat_y_p')^G_ind(m);
                RAmatT_y_m = Rmat_y*(Amat_y_m')^G_ind(m);

                Amat_z_p=real(Amat_z)+ 1i* imag(Amat_z)*grad_dir_p(3);
                Amat_z_m=Amat_z_p'; 
                RAmat_z_p = Rmat_z*Amat_z_p^G_ind(m);
                RAmat_z_m = Rmat_z*Amat_z_m^G_ind(m);   
                RAmatT_z_p = Rmat_z*(Amat_z_p')^G_ind(m);
                RAmatT_z_m = Rmat_z*(Amat_z_m')^G_ind(m);
                
                Ecal_x = Ecal_x*Rmat_x;  % first gradient point is 0
                Ecal_y = Ecal_y*Rmat_y;  % first gradient point is 0
                Ecal_z = Ecal_z*Rmat_z;  % first gradient point is 0
                ProdMat2_x = Rmat_x; 
                ProdMat2_y = Rmat_y; 
                ProdMat2_z = Rmat_z; 
                for i = 2:K(m)-1
                     if sign_vec1(i) >0 
                         Ecal_x = Ecal_x*RAmat_x_p;  
                         Ecal_y = Ecal_y*RAmat_y_p;   
                         Ecal_z = Ecal_z*RAmat_z_p;  
                         ProdMat2_x = ProdMat2_x*RAmatT_x_p;
                         ProdMat2_y = ProdMat2_y*RAmatT_y_p;
                         ProdMat2_z = ProdMat2_z*RAmatT_z_p;
                     else
                         Ecal_x = Ecal_x*RAmat_x_m;
                          Ecal_y = Ecal_y*RAmat_y_m;
                         Ecal_z = Ecal_z*RAmat_z_m;
                         ProdMat2_x = ProdMat2_x*RAmatT_x_m;
                         ProdMat2_y = ProdMat2_y*RAmatT_y_m;
                         ProdMat2_z = ProdMat2_z*RAmatT_z_m;
                     end                        
                end
                 Ecal_x = Ecal_x*Rmat_x;
                 Ecal_y = Ecal_y*Rmat_y;
                 Ecal_z = Ecal_z*Rmat_z;
                 ProdMat2_x = ProdMat2_x*Rmat_x;
                 ProdMat2_y = ProdMat2_y*Rmat_y;
                 ProdMat2_z = ProdMat2_z*Rmat_z;
              Rmid_x=Rmat_x^dstmp(m);
              eX(m,n)=Ecal_x*Rmid_x*ProdMat2_x*Smat_x;
              Rmid_y=Rmat_y^dstmp(m);
              eY(m,n)=Ecal_y*Rmid_y*ProdMat2_y*Smat_y;
              Rmid_z=Rmat_z^dstmp(m);
              eZ(m,n)=Ecal_z*Rmid_z*ProdMat2_z*Smat_z;
                 
            end
        end       
    end
    E = eX.*eY.*eZ;
    %disp(ePerp)
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
   
      % get Ndir directions on a spehere
    eX=zeros(M,RotN); % x direction
    eY = zeros(M,RotN); % y direction
    eZ = zeros(M,RotN); % z direction  
    G_mag = sqrt(Gx.^2+Gy.^2);
    G_ind=round(G_mag./protocol.gstep);
    qhat_X=Gx./G_mag; % qhatx 
    qhat_Y=Gy./G_mag; % qhaty 
    if protocol.mirror == 1
        for m=1:M
            time_vec =tau*(0:K(m)-1);
            sign_vec1 = (-1).^floor((protocol.omegax(m)*time_vec-protocol.phix(m))./pi-1E-10);
            sign_vec2 = (-1).^floor((protocol.omegay(m)*time_vec-protocol.phiy(m))./pi-1E-10);
            for n = 1:RotN
                % Calculating parallel & perpendicular signals
                inv_rot = cube_rotation{n}';
                             
                  
                Ecal_x=Smat_x'; 
                Ecal_y=Smat_y'; 
                Ecal_z=Smat_z'; 
                
                grad_dir_pp = inv_rot*[qhat_X(m); qhat_Y(m); 0];
                grad_dir_pm = inv_rot*[qhat_X(m); -qhat_Y(m); 0];
                grad_dir_mp = inv_rot*[-qhat_X(m); qhat_Y(m); 0];                 

                Amat_x_pp=real(Amat_x)+ 1i* imag(Amat_x)*grad_dir_pp(1);
                Amat_x_pm=real(Amat_x)+ 1i* imag(Amat_x)*grad_dir_pm(1);
                Amat_x_mp=real(Amat_x)+ 1i* imag(Amat_x)*grad_dir_mp(1);
                Amat_x_mm=Amat_x_pp'; 
                RAmat_x_pp = Rmat_x*Amat_x_pp^G_ind(m);
                RAmat_x_pm = Rmat_x*Amat_x_pm^G_ind(m);
                RAmat_x_mp = Rmat_x*Amat_x_mp^G_ind(m);
                RAmat_x_mm = Rmat_x*Amat_x_mm^G_ind(m); 
 
                Amat_y_pp=real(Amat_y)+ 1i* imag(Amat_y)*grad_dir_pp(2);
                Amat_y_pm=real(Amat_y)+ 1i* imag(Amat_y)*grad_dir_pm(2);
                Amat_y_mp=real(Amat_y)+ 1i* imag(Amat_y)*grad_dir_mp(2);
                Amat_y_mm=Amat_y_pp'; 
                RAmat_y_pp = Rmat_y*Amat_y_pp^G_ind(m);
                RAmat_y_pm = Rmat_y*Amat_y_pm^G_ind(m);
                RAmat_y_mp = Rmat_y*Amat_y_mp^G_ind(m);
                RAmat_y_mm = Rmat_y*Amat_y_mm^G_ind(m);

                Amat_z_pp=real(Amat_z)+ 1i* imag(Amat_z)*grad_dir_pp(3);
                Amat_z_pm=real(Amat_z)+ 1i* imag(Amat_z)*grad_dir_pm(3);
                Amat_z_mp=real(Amat_z)+ 1i* imag(Amat_z)*grad_dir_mp(3);
                Amat_z_mm=Amat_z_pp'; 
                RAmat_z_pp = Rmat_z*Amat_z_pp^G_ind(m);
                RAmat_z_pm = Rmat_z*Amat_z_pm^G_ind(m);
                RAmat_z_mp = Rmat_z*Amat_z_mp^G_ind(m);
                RAmat_z_mm = Rmat_z*Amat_z_mm^G_ind(m);
               
                
                Ecal_x = Ecal_x*Rmat_x;  % first gradient point is 0
                Ecal_y = Ecal_y*Rmat_y;  % first gradient point is 0
                Ecal_z = Ecal_z*Rmat_z;  % first gradient point is 0
                for i = 2:K(m)-1
                     if sign_vec1(i) >0 && sign_vec2(i) >0 
                        Ecal_x = Ecal_x*RAmat_x_pp;
                        Ecal_y = Ecal_y*RAmat_y_pp;
                        Ecal_z = Ecal_z*RAmat_z_pp;
                     elseif sign_vec1(i) >0  && sign_vec2(i) <0  
                        Ecal_x = Ecal_x*RAmat_x_pm;
                        Ecal_y = Ecal_y*RAmat_y_pm;
                        Ecal_z = Ecal_z*RAmat_z_pm;
                     elseif sign_vec1(i) <0  && sign_vec2(i) >0 
                         Ecal_x = Ecal_x*RAmat_x_mp;
                         Ecal_y = Ecal_y*RAmat_y_mp;
                         Ecal_z = Ecal_z*RAmat_z_mp;
                     else
                         Ecal_x = Ecal_x*RAmat_x_mm;
                        Ecal_y = Ecal_y*RAmat_y_mm;
                        Ecal_z = Ecal_z*RAmat_z_mm;
                     end                    
                end
                Ecal_x = Ecal_x*Rmat_x;
                 Ecal_y = Ecal_y*Rmat_y;
                 Ecal_z = Ecal_z*Rmat_z;
               
              Rmid_x=Rmat_x^dstmp(m);
              eX(m,n)=Ecal_x*Rmid_x*Ecal_x';  
              Rmid_y=Rmat_y^dstmp(m);
              eY(m,n)=Ecal_y*Rmid_y*Ecal_y'; 
              Rmid_z=Rmat_z^dstmp(m);
              eZ(m,n)=Ecal_z*Rmid_z*Ecal_z'; 
                 
            end
        end       
    elseif protocol.mirror == 0
        for m=1:M
            time_vec = tau*(0:K(m)-1);
            sign_vec1 = (-1).^floor((protocol.omegax(m)*time_vec-protocol.phix(m))./pi-1E-10);
            sign_vec2 = (-1).^floor((protocol.omegay(m)*time_vec-protocol.phiy(m))./pi-1E-10);
             for n = 1:RotN
                % Calculating parallel & perpendicular signals
                inv_rot = cube_rotation{n}';                 
                 Ecal_x=Smat_x'; 
                Ecal_y=Smat_y'; 
                Ecal_z=Smat_z'; 
                  
                grad_dir_pp = inv_rot*[qhat_X(m); qhat_Y(m); 0];
                grad_dir_pm = inv_rot*[qhat_X(m); -qhat_Y(m); 0];
                grad_dir_mp = inv_rot*[-qhat_X(m); qhat_Y(m); 0];                 

                Amat_x_pp=real(Amat_x)+ 1i* imag(Amat_x)*grad_dir_pp(1);
                Amat_x_pm=real(Amat_x)+ 1i* imag(Amat_x)*grad_dir_pm(1);
                Amat_x_mp=real(Amat_x)+ 1i* imag(Amat_x)*grad_dir_mp(1);
                Amat_x_mm=Amat_x_pp'; 
                RAmat_x_pp = Rmat_x*Amat_x_pp^G_ind(m);
                RAmat_x_pm = Rmat_x*Amat_x_pm^G_ind(m);
                RAmat_x_mp = Rmat_x*Amat_x_mp^G_ind(m);
                RAmat_x_mm = Rmat_x*Amat_x_mm^G_ind(m); 
                RAmatT_x_pp = Rmat_x*(Amat_x_pp')^G_ind(m);
                RAmatT_x_pm = Rmat_x*(Amat_x_pm')^G_ind(m);
                RAmatT_x_mp = Rmat_x*(Amat_x_mp')^G_ind(m);
                RAmatT_x_mm = Rmat_x*(Amat_x_mm')^G_ind(m);
 
                Amat_y_pp=real(Amat_y)+ 1i* imag(Amat_y)*grad_dir_pp(2);
                Amat_y_pm=real(Amat_y)+ 1i* imag(Amat_y)*grad_dir_pm(2);
                Amat_y_mp=real(Amat_y)+ 1i* imag(Amat_y)*grad_dir_mp(2);
                Amat_y_mm=Amat_y_pp'; 
                RAmat_y_pp = Rmat_y*Amat_y_pp^G_ind(m);
                RAmat_y_pm = Rmat_y*Amat_y_pm^G_ind(m);
                RAmat_y_mp = Rmat_y*Amat_y_mp^G_ind(m);
                RAmat_y_mm = Rmat_y*Amat_y_mm^G_ind(m);
                RAmatT_y_pp = Rmat_y*(Amat_y_pp')^G_ind(m);
                RAmatT_y_pm = Rmat_y*(Amat_y_pm')^G_ind(m);
                RAmatT_y_mp = Rmat_y*(Amat_y_mp')^G_ind(m);
                RAmatT_y_mm = Rmat_y*(Amat_y_mm')^G_ind(m);

                Amat_z_pp=real(Amat_z)+ 1i* imag(Amat_z)*grad_dir_pp(3);
                Amat_z_pm=real(Amat_z)+ 1i* imag(Amat_z)*grad_dir_pm(3);
                Amat_z_mp=real(Amat_z)+ 1i* imag(Amat_z)*grad_dir_mp(3);
                Amat_z_mm=Amat_z_pp'; 
                RAmat_z_pp = Rmat_z*Amat_z_pp^G_ind(m);
                RAmat_z_pm = Rmat_z*Amat_z_pm^G_ind(m);
                RAmat_z_mp = Rmat_z*Amat_z_mp^G_ind(m);
                RAmat_z_mm = Rmat_z*Amat_z_mm^G_ind(m);
                RAmatT_z_pp = Rmat_z*(Amat_z_pp')^G_ind(m);
                RAmatT_z_pm = Rmat_z*(Amat_z_pm')^G_ind(m);
                RAmatT_z_mp = Rmat_z*(Amat_z_mp')^G_ind(m);
                RAmatT_z_mm = Rmat_z*(Amat_z_mm')^G_ind(m);

                
                Ecal_x = Ecal_x*Rmat_x;  % first gradient point is 0
                Ecal_y = Ecal_y*Rmat_y;  % first gradient point is 0
                Ecal_z = Ecal_z*Rmat_z;  % first gradient point is 0
                ProdMat2_x = Rmat_x; 
                ProdMat2_y = Rmat_y; 
                ProdMat2_z = Rmat_z; 
                for i = 2:K(m)-1
                     if sign_vec1(i) >0 && sign_vec2(i) >0
                        Ecal_x = Ecal_x*RAmat_x_pp;
                        Ecal_y = Ecal_y*RAmat_y_pp;
                        Ecal_z = Ecal_z*RAmat_z_pp;
                         ProdMat2_x = ProdMat2_x*RAmatT_x_pp;
                         ProdMat2_y = ProdMat2_y*RAmatT_y_pp;
                         ProdMat2_z = ProdMat2_z*RAmatT_z_pp;
                     elseif sign_vec1(i) >0 && sign_vec2(i) <0
                        Ecal_x = Ecal_x*RAmat_x_pm;
                        Ecal_y = Ecal_y*RAmat_y_pm;
                        Ecal_z = Ecal_z*RAmat_z_pm;
                         ProdMat2_x = ProdMat2_x*RAmatT_x_pm;
                         ProdMat2_y = ProdMat2_y*RAmatT_y_pm;
                         ProdMat2_z = ProdMat2_z*RAmatT_z_pm;
                     elseif sign_vec1(i) <0 && sign_vec2(i) >0
                         Ecal_x = Ecal_x*RAmat_x_mp;
                         Ecal_y = Ecal_y*RAmat_y_mp;
                         Ecal_z = Ecal_z*RAmat_z_mp;
                         ProdMat2_x = ProdMat2_x*RAmatT_x_mp;
                         ProdMat2_y = ProdMat2_y*RAmatT_y_mp;
                         ProdMat2_z = ProdMat2_z*RAmatT_z_mp;
                     else
                         Ecal_x = Ecal_x*RAmat_x_mm;
                        Ecal_y = Ecal_y*RAmat_y_mm;
                        Ecal_z = Ecal_z*RAmat_z_mm;
                         ProdMat2_x = ProdMat2_x*RAmatT_x_mm;
                         ProdMat2_y = ProdMat2_y*RAmatT_y_mm;
                         ProdMat2_z = ProdMat2_z*RAmatT_z_mm;
                     end                    
                end
                Ecal_x = Ecal_x*Rmat_x;
                 Ecal_y = Ecal_y*Rmat_y;
                 Ecal_z = Ecal_z*Rmat_z;
                 ProdMat2_x = ProdMat2_x*Rmat_x;
                 ProdMat2_y = ProdMat2_y*Rmat_y;
                 ProdMat2_z = ProdMat2_z*Rmat_z;
              Rmid_x=Rmat_x^dstmp(m);
              eX(m,n)=Ecal_x*Rmid_x*ProdMat2_x*Smat_x;
              Rmid_y=Rmat_y^dstmp(m);
              eY(m,n)=Ecal_y*Rmid_y*ProdMat2_y*Smat_y;
              Rmid_z=Rmat_z^dstmp(m);
              eZ(m,n)=Ecal_z*Rmid_z*ProdMat2_z*Smat_z;
                 
            end
        end       
    end
    E = eX.*eY.*eZ;
    %disp(ePerp)
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
   
      % get Ndir directions on a spehere
    eX=zeros(M,RotN); % x direction
    eY = zeros(M,RotN); % y direction
    eZ = zeros(M,RotN); % z direction  
    G_mag = sqrt(Gx.^2+Gz.^2);
    G_ind=round(G_mag./protocol.gstep);
    qhat_X=Gx./G_mag; % qhatx 
    qhat_Z=Gz./G_mag; % qhatz 
    if protocol.mirror == 1
        for m=1:M
            time_vec = tau*(0:K(m)-1);
            sign_vec1 = (-1).^floor((protocol.omegax(m)*time_vec-protocol.phix(m))./pi-1E-10);
            sign_vec2 = (-1).^floor((protocol.omegaz(m)*time_vec-protocol.phiz(m))./pi-1E-10);
            for n = 1:RotN
                % Calculating parallel & perpendicular signals
                inv_rot = cube_rotation{n}';
                           
                  
                Ecal_x=Smat_x'; 
                Ecal_y=Smat_y'; 
                Ecal_z=Smat_z'; 
                
                grad_dir_pp = inv_rot*[qhat_X(m); 0; qhat_Z(m); ];
                grad_dir_pm = inv_rot*[qhat_X(m); 0; -qhat_Z(m); ];
                grad_dir_mp = inv_rot*[-qhat_X(m); 0; qhat_Z(m); ];                 

                Amat_x_pp=real(Amat_x)+ 1i* imag(Amat_x)*grad_dir_pp(1);
                Amat_x_pm=real(Amat_x)+ 1i* imag(Amat_x)*grad_dir_pm(1);
                Amat_x_mp=real(Amat_x)+ 1i* imag(Amat_x)*grad_dir_mp(1);
                Amat_x_mm=Amat_x_pp'; 
                RAmat_x_pp = Rmat_x*Amat_x_pp^G_ind(m);
                RAmat_x_pm = Rmat_x*Amat_x_pm^G_ind(m);
                RAmat_x_mp = Rmat_x*Amat_x_mp^G_ind(m);
                RAmat_x_mm = Rmat_x*Amat_x_mm^G_ind(m); 
 
                Amat_y_pp=real(Amat_y)+ 1i* imag(Amat_y)*grad_dir_pp(2);
                Amat_y_pm=real(Amat_y)+ 1i* imag(Amat_y)*grad_dir_pm(2);
                Amat_y_mp=real(Amat_y)+ 1i* imag(Amat_y)*grad_dir_mp(2);
                Amat_y_mm=Amat_y_pp'; 
                RAmat_y_pp = Rmat_y*Amat_y_pp^G_ind(m);
                RAmat_y_pm = Rmat_y*Amat_y_pm^G_ind(m);
                RAmat_y_mp = Rmat_y*Amat_y_mp^G_ind(m);
                RAmat_y_mm = Rmat_y*Amat_y_mm^G_ind(m);

                Amat_z_pp=real(Amat_z)+ 1i* imag(Amat_z)*grad_dir_pp(3);
                Amat_z_pm=real(Amat_z)+ 1i* imag(Amat_z)*grad_dir_pm(3);
                Amat_z_mp=real(Amat_z)+ 1i* imag(Amat_z)*grad_dir_mp(3);
                Amat_z_mm=Amat_z_pp'; 
                RAmat_z_pp = Rmat_z*Amat_z_pp^G_ind(m);
                RAmat_z_pm = Rmat_z*Amat_z_pm^G_ind(m);
                RAmat_z_mp = Rmat_z*Amat_z_mp^G_ind(m);
                RAmat_z_mm = Rmat_z*Amat_z_mm^G_ind(m);
               
                
                Ecal_x = Ecal_x*Rmat_x;  % first gradient point is 0
                Ecal_y = Ecal_y*Rmat_y;  % first gradient point is 0
                Ecal_z = Ecal_z*Rmat_z;  % first gradient point is 0
                for i = 2:K(m)-1
                     if sign_vec1(i) >0 && sign_vec2(i) >0
                        Ecal_x = Ecal_x*RAmat_x_pp;
                        Ecal_y = Ecal_y*RAmat_y_pp;
                        Ecal_z = Ecal_z*RAmat_z_pp;
                     elseif sign_vec1(i) >0 && sign_vec2(i) <0
                        Ecal_x = Ecal_x*RAmat_x_pm;
                        Ecal_y = Ecal_y*RAmat_y_pm;
                        Ecal_z = Ecal_z*RAmat_z_pm;
                     elseif sign_vec1(i) <0 && sign_vec2(i) >0
                         Ecal_x = Ecal_x*RAmat_x_mp;
                         Ecal_y = Ecal_y*RAmat_y_mp;
                         Ecal_z = Ecal_z*RAmat_z_mp;
                     else
                         Ecal_x = Ecal_x*RAmat_x_mm;
                        Ecal_y = Ecal_y*RAmat_y_mm;
                        Ecal_z = Ecal_z*RAmat_z_mm;
                     end                    
                end
                Ecal_x = Ecal_x*Rmat_x;
                 Ecal_y = Ecal_y*Rmat_y;
                 Ecal_z = Ecal_z*Rmat_z;
              
              Rmid_x=Rmat_x^dstmp(m);
              eX(m,n)=Ecal_x*Rmid_x*Ecal_x';  
              Rmid_y=Rmat_y^dstmp(m);
              eY(m,n)=Ecal_y*Rmid_y*Ecal_y'; 
              Rmid_z=Rmat_z^dstmp(m);
              eZ(m,n)=Ecal_z*Rmid_z*Ecal_z'; 
                 
            end
        end       
    elseif protocol.mirror == 0
        for m=1:M
            time_vec = tau*(0:K(m)-1);
            sign_vec1 = (-1).^floor((protocol.omegax(m)*time_vec-protocol.phix(m))./pi-1E-10);
            sign_vec2 = (-1).^floor((protocol.omegaz(m)*time_vec-protocol.phiz(m))./pi-1E-10);
             for n = 1:RotN
                % Calculating parallel & perpendicular signals
                inv_rot = cube_rotation{n}';
               
                 Ecal_x=Smat_x'; 
                Ecal_y=Smat_y'; 
                Ecal_z=Smat_z'; 
                  
                grad_dir_pp = inv_rot*[qhat_X(m); 0; qhat_Z(m); ];
                grad_dir_pm = inv_rot*[qhat_X(m); 0; -qhat_Z(m); ];
                grad_dir_mp = inv_rot*[-qhat_X(m); 0; qhat_Z(m); ];                

                Amat_x_pp=real(Amat_x)+ 1i* imag(Amat_x)*grad_dir_pp(1);
                Amat_x_pm=real(Amat_x)+ 1i* imag(Amat_x)*grad_dir_pm(1);
                Amat_x_mp=real(Amat_x)+ 1i* imag(Amat_x)*grad_dir_mp(1);
                Amat_x_mm=Amat_x_pp'; 
                RAmat_x_pp = Rmat_x*Amat_x_pp^G_ind(m);
                RAmat_x_pm = Rmat_x*Amat_x_pm^G_ind(m);
                RAmat_x_mp = Rmat_x*Amat_x_mp^G_ind(m);
                RAmat_x_mm = Rmat_x*Amat_x_mm^G_ind(m); 
                RAmatT_x_pp = Rmat_x*(Amat_x_pp')^G_ind(m);
                RAmatT_x_pm = Rmat_x*(Amat_x_pm')^G_ind(m);
                RAmatT_x_mp = Rmat_x*(Amat_x_mp')^G_ind(m);
                RAmatT_x_mm = Rmat_x*(Amat_x_mm')^G_ind(m);
 
                Amat_y_pp=real(Amat_y)+ 1i* imag(Amat_y)*grad_dir_pp(2);
                Amat_y_pm=real(Amat_y)+ 1i* imag(Amat_y)*grad_dir_pm(2);
                Amat_y_mp=real(Amat_y)+ 1i* imag(Amat_y)*grad_dir_mp(2);
                Amat_y_mm=Amat_y_pp'; 
                RAmat_y_pp = Rmat_y*Amat_y_pp^G_ind(m);
                RAmat_y_pm = Rmat_y*Amat_y_pm^G_ind(m);
                RAmat_y_mp = Rmat_y*Amat_y_mp^G_ind(m);
                RAmat_y_mm = Rmat_y*Amat_y_mm^G_ind(m);
                RAmatT_y_pp = Rmat_y*(Amat_y_pp')^G_ind(m);
                RAmatT_y_pm = Rmat_y*(Amat_y_pm')^G_ind(m);
                RAmatT_y_mp = Rmat_y*(Amat_y_mp')^G_ind(m);
                RAmatT_y_mm = Rmat_y*(Amat_y_mm')^G_ind(m);

                Amat_z_pp=real(Amat_z)+ 1i* imag(Amat_z)*grad_dir_pp(3);
                Amat_z_pm=real(Amat_z)+ 1i* imag(Amat_z)*grad_dir_pm(3);
                Amat_z_mp=real(Amat_z)+ 1i* imag(Amat_z)*grad_dir_mp(3);
                Amat_z_mm=Amat_z_pp'; 
                RAmat_z_pp = Rmat_z*Amat_z_pp^G_ind(m);
                RAmat_z_pm = Rmat_z*Amat_z_pm^G_ind(m);
                RAmat_z_mp = Rmat_z*Amat_z_mp^G_ind(m);
                RAmat_z_mm = Rmat_z*Amat_z_mm^G_ind(m);
                RAmatT_z_pp = Rmat_z*(Amat_z_pp')^G_ind(m);
                RAmatT_z_pm = Rmat_z*(Amat_z_pm')^G_ind(m);
                RAmatT_z_mp = Rmat_z*(Amat_z_mp')^G_ind(m);
                RAmatT_z_mm = Rmat_z*(Amat_z_mm')^G_ind(m);

                
                Ecal_x = Ecal_x*Rmat_x;  % first gradient point is 0
                Ecal_y = Ecal_y*Rmat_y;  % first gradient point is 0
                Ecal_z = Ecal_z*Rmat_z;  % first gradient point is 0
                ProdMat2_x = Rmat_x; 
                ProdMat2_y = Rmat_y; 
                ProdMat2_z = Rmat_z; 
                for i = 2:K(m)-1
                     if sign_vec1(i) >0 && sign_vec2(i) >0
                        Ecal_x = Ecal_x*RAmat_x_pp;
                        Ecal_y = Ecal_y*RAmat_y_pp;
                        Ecal_z = Ecal_z*RAmat_z_pp;
                         ProdMat2_x = ProdMat2_x*RAmatT_x_pp;
                         ProdMat2_y = ProdMat2_y*RAmatT_y_pp;
                         ProdMat2_z = ProdMat2_z*RAmatT_z_pp;
                     elseif sign_vec1(i) >0 && sign_vec2(i) <0
                        Ecal_x = Ecal_x*RAmat_x_pm;
                        Ecal_y = Ecal_y*RAmat_y_pm;
                        Ecal_z = Ecal_z*RAmat_z_pm;
                         ProdMat2_x = ProdMat2_x*RAmatT_x_pm;
                         ProdMat2_y = ProdMat2_y*RAmatT_y_pm;
                         ProdMat2_z = ProdMat2_z*RAmatT_z_pm;
                     elseif sign_vec1(i) <0 && sign_vec2(i) >0
                         Ecal_x = Ecal_x*RAmat_x_mp;
                         Ecal_y = Ecal_y*RAmat_y_mp;
                         Ecal_z = Ecal_z*RAmat_z_mp;
                         ProdMat2_x = ProdMat2_x*RAmatT_x_mp;
                         ProdMat2_y = ProdMat2_y*RAmatT_y_mp;
                         ProdMat2_z = ProdMat2_z*RAmatT_z_mp;
                     else
                         Ecal_x = Ecal_x*RAmat_x_mm;
                        Ecal_y = Ecal_y*RAmat_y_mm;
                        Ecal_z = Ecal_z*RAmat_z_mm;
                         ProdMat2_x = ProdMat2_x*RAmatT_x_mm;
                         ProdMat2_y = ProdMat2_y*RAmatT_y_mm;
                         ProdMat2_z = ProdMat2_z*RAmatT_z_mm;
                     end                    
                end
                Ecal_x = Ecal_x*Rmat_x;
                 Ecal_y = Ecal_y*Rmat_y;
                 Ecal_z = Ecal_z*Rmat_z;
                 ProdMat2_x = ProdMat2_x*Rmat_x;
                 ProdMat2_y = ProdMat2_y*Rmat_y;
                 ProdMat2_z = ProdMat2_z*Rmat_z;
             Rmid_x=Rmat_x^dstmp(m);
              eX(m,n)=Ecal_x*Rmid_x*ProdMat2_x*Smat_x;
              Rmid_y=Rmat_y^dstmp(m);
              eY(m,n)=Ecal_y*Rmid_y*ProdMat2_y*Smat_y;
              Rmid_z=Rmat_z^dstmp(m);
              eZ(m,n)=Ecal_z*Rmid_z*ProdMat2_z*Smat_z;
                 
            end
        end       
    end
    E = eX.*eY.*eZ;
    %disp(ePerp)
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
   
      % get Ndir directions on a spehere
    eX=zeros(M,RotN); % x direction
    eY = zeros(M,RotN); % y direction
    eZ = zeros(M,RotN); % z direction  
    G_mag = sqrt(Gx.^2+Gy.^2+Gz.^2);
    G_ind=round(G_mag./protocol.gstep);
    qhat_X=Gx./G_mag; % qhatx 
    qhat_Y=Gy./G_mag; % qhaty 
    qhat_Z=Gz./G_mag; % qhatz 
    if protocol.mirror == 1
        for m=1:M
            time_vec = tau*(0:K(m)-1);
            sign_vec1 = (-1).^floor((protocol.omegax(m)*time_vec-protocol.phix(m))./pi-1E-10);
            sign_vec2 = (-1).^floor((protocol.omegay(m)*time_vec-protocol.phiy(m))./pi-1E-10);
            sign_vec3 = (-1).^floor((protocol.omegaz(m)*time_vec-protocol.phiz(m))./pi-1E-10);
            for n = 1:RotN
                % Calculating parallel & perpendicular signals
                inv_rot = cube_rotation{n}';
                             
                  
                Ecal_x=Smat_x'; 
                Ecal_y=Smat_y'; 
                Ecal_z=Smat_z'; 
                
                grad_dir_ppp = inv_rot*[qhat_X(m); qhat_Y(m); qhat_Z(m); ];
                grad_dir_ppm = inv_rot*[qhat_X(m); qhat_Y(m); -qhat_Z(m); ];
                grad_dir_pmp = inv_rot*[qhat_X(m); -qhat_Y(m); qhat_Z(m); ];
                grad_dir_pmm = inv_rot*[qhat_X(m); -qhat_Y(m); -qhat_Z(m); ];
                grad_dir_mpp = inv_rot*[-qhat_X(m); qhat_Y(m); qhat_Z(m); ];
                grad_dir_mpm = inv_rot*[-qhat_X(m); qhat_Y(m); -qhat_Z(m); ];
                grad_dir_mmp = inv_rot*[-qhat_X(m); -qhat_Y(m); qhat_Z(m); ];
                grad_dir_mmm = inv_rot*[-qhat_X(m); -qhat_Y(m); -qhat_Z(m); ];
               

                Amat_x_ppp=real(Amat_x)+ 1i* imag(Amat_x)*grad_dir_ppp(1);
                Amat_x_ppm=real(Amat_x)+ 1i* imag(Amat_x)*grad_dir_ppm(1);
                Amat_x_pmp=real(Amat_x)+ 1i* imag(Amat_x)*grad_dir_pmp(1);
                Amat_x_pmm=real(Amat_x)+ 1i* imag(Amat_x)*grad_dir_pmm(1); 
                Amat_x_mpp=real(Amat_x)+ 1i* imag(Amat_x)*grad_dir_mpp(1);
                Amat_x_mpm=real(Amat_x)+ 1i* imag(Amat_x)*grad_dir_mpm(1);
                Amat_x_mmp=real(Amat_x)+ 1i* imag(Amat_x)*grad_dir_mmp(1);
                Amat_x_mmm=real(Amat_x)+ 1i* imag(Amat_x)*grad_dir_mmm(1);
                RAmat_x_ppp = Rmat_x*Amat_x_ppp^G_ind(m);
                RAmat_x_ppm = Rmat_x*Amat_x_ppm^G_ind(m);
                RAmat_x_pmp = Rmat_x*Amat_x_pmp^G_ind(m);
                RAmat_x_pmm = Rmat_x*Amat_x_pmm^G_ind(m); 
                RAmat_x_mpp = Rmat_x*Amat_x_mpp^G_ind(m);
                RAmat_x_mpm = Rmat_x*Amat_x_mpm^G_ind(m);
                RAmat_x_mmp = Rmat_x*Amat_x_mmp^G_ind(m);
                RAmat_x_mmm = Rmat_x*Amat_x_mmm^G_ind(m);
 
                Amat_y_ppp=real(Amat_y)+ 1i* imag(Amat_y)*grad_dir_ppp(2);
                Amat_y_ppm=real(Amat_y)+ 1i* imag(Amat_y)*grad_dir_ppm(2);
                Amat_y_pmp=real(Amat_y)+ 1i* imag(Amat_y)*grad_dir_pmp(2);
                Amat_y_pmm=real(Amat_y)+ 1i* imag(Amat_y)*grad_dir_pmm(2); 
                Amat_y_mpp=real(Amat_y)+ 1i* imag(Amat_y)*grad_dir_mpp(2);
                Amat_y_mpm=real(Amat_y)+ 1i* imag(Amat_y)*grad_dir_mpm(2);
                Amat_y_mmp=real(Amat_y)+ 1i* imag(Amat_y)*grad_dir_mmp(2);
                Amat_y_mmm=real(Amat_y)+ 1i* imag(Amat_y)*grad_dir_mmm(2);
                RAmat_y_ppp = Rmat_y*Amat_y_ppp^G_ind(m);
                RAmat_y_ppm = Rmat_y*Amat_y_ppm^G_ind(m);
                RAmat_y_pmp = Rmat_y*Amat_y_pmp^G_ind(m);
                RAmat_y_pmm = Rmat_y*Amat_y_pmm^G_ind(m); 
                RAmat_y_mpp = Rmat_y*Amat_y_mpp^G_ind(m);
                RAmat_y_mpm = Rmat_y*Amat_y_mpm^G_ind(m);
                RAmat_y_mmp = Rmat_y*Amat_y_mmp^G_ind(m);
                RAmat_y_mmm = Rmat_y*Amat_y_mmm^G_ind(m);

                Amat_z_ppp=real(Amat_z)+ 1i* imag(Amat_z)*grad_dir_ppp(3);
                Amat_z_ppm=real(Amat_z)+ 1i* imag(Amat_z)*grad_dir_ppm(3);
                Amat_z_pmp=real(Amat_z)+ 1i* imag(Amat_z)*grad_dir_pmp(3);
                Amat_z_pmm=real(Amat_z)+ 1i* imag(Amat_z)*grad_dir_pmm(3); 
                Amat_z_mpp=real(Amat_z)+ 1i* imag(Amat_z)*grad_dir_mpp(3);
                Amat_z_mpm=real(Amat_z)+ 1i* imag(Amat_z)*grad_dir_mpm(3);
                Amat_z_mmp=real(Amat_z)+ 1i* imag(Amat_z)*grad_dir_mmp(3);
                Amat_z_mmm=real(Amat_z)+ 1i* imag(Amat_z)*grad_dir_mmm(3);
                RAmat_z_ppp = Rmat_z*Amat_z_ppp^G_ind(m);
                RAmat_z_ppm = Rmat_z*Amat_z_ppm^G_ind(m);
                RAmat_z_pmp = Rmat_z*Amat_z_pmp^G_ind(m);
                RAmat_z_pmm = Rmat_z*Amat_z_pmm^G_ind(m); 
                RAmat_z_mpp = Rmat_z*Amat_z_mpp^G_ind(m);
                RAmat_z_mpm = Rmat_z*Amat_z_mpm^G_ind(m);
                RAmat_z_mmp = Rmat_z*Amat_z_mmp^G_ind(m);
                RAmat_z_mmm = Rmat_z*Amat_z_mmm^G_ind(m);
               
                
                Ecal_x = Ecal_x*Rmat_x;  % first gradient point is 0
                Ecal_y = Ecal_y*Rmat_y;  % first gradient point is 0
                Ecal_z = Ecal_z*Rmat_z;  % first gradient point is 0
                for i = 2:K(m)-1
                     if sign_vec1(i) >0 && sign_vec2(i) >0 && sign_vec3(i) >0
                        Ecal_x = Ecal_x*RAmat_x_ppp;
                        Ecal_y = Ecal_y*RAmat_y_ppp;
                        Ecal_z = Ecal_z*RAmat_z_ppp;
                     elseif sign_vec1(i) >0 && sign_vec2(i) >0 && sign_vec3(i) <0
                        Ecal_x = Ecal_x*RAmat_x_ppm;
                        Ecal_y = Ecal_y*RAmat_y_ppm;
                        Ecal_z = Ecal_z*RAmat_z_ppm;
                     elseif sign_vec1(i) >0 && sign_vec2(i) <0 && sign_vec3(i) >0
                        Ecal_x = Ecal_x*RAmat_x_pmp;
                        Ecal_y = Ecal_y*RAmat_y_pmp;
                        Ecal_z = Ecal_z*RAmat_z_pmp;
                     elseif sign_vec1(i) >0 && sign_vec2(i) <0 && sign_vec3(i) <0
                        Ecal_x = Ecal_x*RAmat_x_pmm;
                        Ecal_y = Ecal_y*RAmat_y_pmm;
                        Ecal_z = Ecal_z*RAmat_z_pmm;
                     elseif sign_vec1(i) <0 && sign_vec2(i) >0 && sign_vec3(i) >0
                        Ecal_x = Ecal_x*RAmat_x_mpp;
                        Ecal_y = Ecal_y*RAmat_y_mpp;
                        Ecal_z = Ecal_z*RAmat_z_mpp;
                     elseif sign_vec1(i) <0 && sign_vec2(i) >0 && sign_vec3(i) <0
                        Ecal_x = Ecal_x*RAmat_x_mpm;
                        Ecal_y = Ecal_y*RAmat_y_mpm;
                        Ecal_z = Ecal_z*RAmat_z_mpm;
                     elseif sign_vec1(i) <0 && sign_vec2(i) <0 && sign_vec3(i) >0
                        Ecal_x = Ecal_x*RAmat_x_mmp;
                        Ecal_y = Ecal_y*RAmat_y_mmp;
                        Ecal_z = Ecal_z*RAmat_z_mmp;
                     else 
                        Ecal_x = Ecal_x*RAmat_x_mmm; 
                        Ecal_y = Ecal_y*RAmat_y_mmm; 
                        Ecal_z = Ecal_z*RAmat_z_mmm;
                     end   
                end
                Ecal_x = Ecal_x*Rmat_x;
                 Ecal_y = Ecal_y*Rmat_y;
                 Ecal_z = Ecal_z*Rmat_z;                
              Rmid_x=Rmat_x^dstmp(m);
              eX(m,n)=Ecal_x*Rmid_x*Ecal_x';  
              Rmid_y=Rmat_y^dstmp(m);
              eY(m,n)=Ecal_y*Rmid_y*Ecal_y'; 
              Rmid_z=Rmat_z^dstmp(m);
              eZ(m,n)=Ecal_z*Rmid_z*Ecal_z'; 
                 
            end
        end       
    elseif protocol.mirror == 0
        for m=1:M
            time_vec = tau*(0:K(m)-1);
            sign_vec1 = (-1).^floor((protocol.omegax(m)*time_vec-protocol.phix(m))./pi-1E-10);
            sign_vec2 = (-1).^floor((protocol.omegay(m)*time_vec-protocol.phiy(m))./pi-1E-10);
            sign_vec3 = (-1).^floor((protocol.omegaz(m)*time_vec-protocol.phiz(m))./pi-1E-10);
             for n = 1:RotN
                % Calculating parallel & perpendicular signals
                inv_rot = cube_rotation{n}';
                
                Ecal_x=Smat_x'; 
                Ecal_y=Smat_y'; 
                Ecal_z=Smat_z'; 
                
                grad_dir_ppp = inv_rot*[qhat_X(m); qhat_Y(m); qhat_Z(m); ];
                grad_dir_ppm = inv_rot*[qhat_X(m); qhat_Y(m); -qhat_Z(m); ];
                grad_dir_pmp = inv_rot*[qhat_X(m); -qhat_Y(m); qhat_Z(m); ];
                grad_dir_pmm = inv_rot*[qhat_X(m); -qhat_Y(m); -qhat_Z(m); ];
                grad_dir_mpp = inv_rot*[-qhat_X(m); qhat_Y(m); qhat_Z(m); ];
                grad_dir_mpm = inv_rot*[-qhat_X(m); qhat_Y(m); -qhat_Z(m); ];
                grad_dir_mmp = inv_rot*[-qhat_X(m); -qhat_Y(m); qhat_Z(m); ];
                grad_dir_mmm = inv_rot*[-qhat_X(m); -qhat_Y(m); -qhat_Z(m); ];
               

                Amat_x_ppp=real(Amat_x)+ 1i* imag(Amat_x)*grad_dir_ppp(1);
                Amat_x_ppm=real(Amat_x)+ 1i* imag(Amat_x)*grad_dir_ppm(1);
                Amat_x_pmp=real(Amat_x)+ 1i* imag(Amat_x)*grad_dir_pmp(1);
                Amat_x_pmm=real(Amat_x)+ 1i* imag(Amat_x)*grad_dir_pmm(1); 
                Amat_x_mpp=real(Amat_x)+ 1i* imag(Amat_x)*grad_dir_mpp(1);
                Amat_x_mpm=real(Amat_x)+ 1i* imag(Amat_x)*grad_dir_mpm(1);
                Amat_x_mmp=real(Amat_x)+ 1i* imag(Amat_x)*grad_dir_mmp(1);
                Amat_x_mmm=real(Amat_x)+ 1i* imag(Amat_x)*grad_dir_mmm(1);
                RAmat_x_ppp = Rmat_x*Amat_x_ppp^G_ind(m);
                RAmat_x_ppm = Rmat_x*Amat_x_ppm^G_ind(m);
                RAmat_x_pmp = Rmat_x*Amat_x_pmp^G_ind(m);
                RAmat_x_pmm = Rmat_x*Amat_x_pmm^G_ind(m); 
                RAmat_x_mpp = Rmat_x*Amat_x_mpp^G_ind(m);
                RAmat_x_mpm = Rmat_x*Amat_x_mpm^G_ind(m);
                RAmat_x_mmp = Rmat_x*Amat_x_mmp^G_ind(m);
                RAmat_x_mmm = Rmat_x*Amat_x_mmm^G_ind(m);
                RAmatT_x_ppp = Rmat_x*(Amat_x_ppp')^G_ind(m);
                RAmatT_x_ppm = Rmat_x*(Amat_x_ppm')^G_ind(m);
                RAmatT_x_pmp = Rmat_x*(Amat_x_pmp')^G_ind(m);
                RAmatT_x_pmm = Rmat_x*(Amat_x_pmm')^G_ind(m);
                RAmatT_x_mpp = Rmat_x*(Amat_x_mpp')^G_ind(m);
                RAmatT_x_mpm = Rmat_x*(Amat_x_mpm')^G_ind(m);
                RAmatT_x_mmp = Rmat_x*(Amat_x_mmp')^G_ind(m);
                RAmatT_x_mmm = Rmat_x*(Amat_x_mmm')^G_ind(m);
                
                
 
                Amat_y_ppp=real(Amat_y)+ 1i* imag(Amat_y)*grad_dir_ppp(2);
                Amat_y_ppm=real(Amat_y)+ 1i* imag(Amat_y)*grad_dir_ppm(2);
                Amat_y_pmp=real(Amat_y)+ 1i* imag(Amat_y)*grad_dir_pmp(2);
                Amat_y_pmm=real(Amat_y)+ 1i* imag(Amat_y)*grad_dir_pmm(2); 
                Amat_y_mpp=real(Amat_y)+ 1i* imag(Amat_y)*grad_dir_mpp(2);
                Amat_y_mpm=real(Amat_y)+ 1i* imag(Amat_y)*grad_dir_mpm(2);
                Amat_y_mmp=real(Amat_y)+ 1i* imag(Amat_y)*grad_dir_mmp(2);
                Amat_y_mmm=real(Amat_y)+ 1i* imag(Amat_y)*grad_dir_mmm(2);
                RAmat_y_ppp = Rmat_y*Amat_y_ppp^G_ind(m);
                RAmat_y_ppm = Rmat_y*Amat_y_ppm^G_ind(m);
                RAmat_y_pmp = Rmat_y*Amat_y_pmp^G_ind(m);
                RAmat_y_pmm = Rmat_y*Amat_y_pmm^G_ind(m); 
                RAmat_y_mpp = Rmat_y*Amat_y_mpp^G_ind(m);
                RAmat_y_mpm = Rmat_y*Amat_y_mpm^G_ind(m);
                RAmat_y_mmp = Rmat_y*Amat_y_mmp^G_ind(m);
                RAmat_y_mmm = Rmat_y*Amat_y_mmm^G_ind(m);
                RAmatT_y_ppp = Rmat_y*(Amat_y_ppp')^G_ind(m);
                RAmatT_y_ppm = Rmat_y*(Amat_y_ppm')^G_ind(m);
                RAmatT_y_pmp = Rmat_y*(Amat_y_pmp')^G_ind(m);
                RAmatT_y_pmm = Rmat_y*(Amat_y_pmm')^G_ind(m);
                RAmatT_y_mpp = Rmat_y*(Amat_y_mpp')^G_ind(m);
                RAmatT_y_mpm = Rmat_y*(Amat_y_mpm')^G_ind(m);
                RAmatT_y_mmp = Rmat_y*(Amat_y_mmp')^G_ind(m);
                RAmatT_y_mmm = Rmat_y*(Amat_y_mmm')^G_ind(m);

                Amat_z_ppp=real(Amat_z)+ 1i* imag(Amat_z)*grad_dir_ppp(3);
                Amat_z_ppm=real(Amat_z)+ 1i* imag(Amat_z)*grad_dir_ppm(3);
                Amat_z_pmp=real(Amat_z)+ 1i* imag(Amat_z)*grad_dir_pmp(3);
                Amat_z_pmm=real(Amat_z)+ 1i* imag(Amat_z)*grad_dir_pmm(3); 
                Amat_z_mpp=real(Amat_z)+ 1i* imag(Amat_z)*grad_dir_mpp(3);
                Amat_z_mpm=real(Amat_z)+ 1i* imag(Amat_z)*grad_dir_mpm(3);
                Amat_z_mmp=real(Amat_z)+ 1i* imag(Amat_z)*grad_dir_mmp(3);
                Amat_z_mmm=real(Amat_z)+ 1i* imag(Amat_z)*grad_dir_mmm(3);
                RAmat_z_ppp = Rmat_z*Amat_z_ppp^G_ind(m);
                RAmat_z_ppm = Rmat_z*Amat_z_ppm^G_ind(m);
                RAmat_z_pmp = Rmat_z*Amat_z_pmp^G_ind(m);
                RAmat_z_pmm = Rmat_z*Amat_z_pmm^G_ind(m); 
                RAmat_z_mpp = Rmat_z*Amat_z_mpp^G_ind(m);
                RAmat_z_mpm = Rmat_z*Amat_z_mpm^G_ind(m);
                RAmat_z_mmp = Rmat_z*Amat_z_mmp^G_ind(m);
                RAmat_z_mmm = Rmat_z*Amat_z_mmm^G_ind(m);
                RAmatT_z_ppp = Rmat_z*(Amat_z_ppp')^G_ind(m);
                RAmatT_z_ppm = Rmat_z*(Amat_z_ppm')^G_ind(m);
                RAmatT_z_pmp = Rmat_z*(Amat_z_pmp')^G_ind(m);
                RAmatT_z_pmm = Rmat_z*(Amat_z_pmm')^G_ind(m);
                RAmatT_z_mpp = Rmat_z*(Amat_z_mpp')^G_ind(m);
                RAmatT_z_mpm = Rmat_z*(Amat_z_mpm')^G_ind(m);
                RAmatT_z_mmp = Rmat_z*(Amat_z_mmp')^G_ind(m);
                RAmatT_z_mmm = Rmat_z*(Amat_z_mmm')^G_ind(m);
               
                
                Ecal_x = Ecal_x*Rmat_x;  % first gradient point is 0
                Ecal_y = Ecal_y*Rmat_y;  % first gradient point is 0
                Ecal_z = Ecal_z*Rmat_z;  % first gradient point is 0
                ProdMat2_x = Rmat_x; 
                ProdMat2_y = Rmat_y; 
                ProdMat2_z = Rmat_z; 
                for i = 2:K(m)-1
                     if sign_vec1(i) >0 && sign_vec2(i) >0 && sign_vec3(i) >0
                        Ecal_x = Ecal_x*RAmat_x_ppp;
                        Ecal_y = Ecal_y*RAmat_y_ppp;
                        Ecal_z = Ecal_z*RAmat_z_ppp;
                        ProdMat2_x = ProdMat2_x*RAmatT_x_ppp;
                        ProdMat2_y = ProdMat2_y*RAmatT_y_ppp;
                        ProdMat2_z = ProdMat2_z*RAmatT_z_ppp;
                     elseif sign_vec1(i) >0 && sign_vec2(i) >0 && sign_vec3(i) <0
                        Ecal_x = Ecal_x*RAmat_x_ppm;
                        Ecal_y = Ecal_y*RAmat_y_ppm;
                        Ecal_z = Ecal_z*RAmat_z_ppm;
                        ProdMat2_x = ProdMat2_x*RAmatT_x_ppm;
                        ProdMat2_y = ProdMat2_y*RAmatT_y_ppm;
                        ProdMat2_z = ProdMat2_z*RAmatT_z_ppm;
                     elseif sign_vec1(i) >0 && sign_vec2(i) <0 && sign_vec3(i) >0
                        Ecal_x = Ecal_x*RAmat_x_pmp;
                        Ecal_y = Ecal_y*RAmat_y_pmp;
                        Ecal_z = Ecal_z*RAmat_z_pmp;
                        ProdMat2_x = ProdMat2_x*RAmatT_x_pmp;
                        ProdMat2_y = ProdMat2_y*RAmatT_y_pmp;
                        ProdMat2_z = ProdMat2_z*RAmatT_z_pmp;
                     elseif sign_vec1(i) >0 && sign_vec2(i) <0 && sign_vec3(i) <0
                        Ecal_x = Ecal_x*RAmat_x_pmm;
                        Ecal_y = Ecal_y*RAmat_y_pmm;
                        Ecal_z = Ecal_z*RAmat_z_pmm;
                        ProdMat2_x = ProdMat2_x*RAmatT_x_pmm;
                        ProdMat2_y = ProdMat2_y*RAmatT_y_pmm;
                        ProdMat2_z = ProdMat2_z*RAmatT_z_pmm;
                     elseif sign_vec1(i) <0 && sign_vec2(i) >0 && sign_vec3(i) >0
                        Ecal_x = Ecal_x*RAmat_x_mpp;
                        Ecal_y = Ecal_y*RAmat_y_mpp;
                        Ecal_z = Ecal_z*RAmat_z_mpp;
                        ProdMat2_x = ProdMat2_x*RAmatT_x_mpp;
                        ProdMat2_y = ProdMat2_y*RAmatT_y_mpp;
                        ProdMat2_z = ProdMat2_z*RAmatT_z_mpp;
                     elseif sign_vec1(i) <0 && sign_vec2(i) >0 && sign_vec3(i) <0
                        Ecal_x = Ecal_x*RAmat_x_mpm;
                        Ecal_y = Ecal_y*RAmat_y_mpm;
                        Ecal_z = Ecal_z*RAmat_z_mpm;
                        ProdMat2_x = ProdMat2_x*RAmatT_x_mpm;
                        ProdMat2_y = ProdMat2_y*RAmatT_y_mpm;
                        ProdMat2_z = ProdMat2_z*RAmatT_z_mpm;
                     elseif sign_vec1(i) <0 && sign_vec2(i) <0 && sign_vec3(i) >0
                        Ecal_x = Ecal_x*RAmat_x_mmp;
                        Ecal_y = Ecal_y*RAmat_y_mmp;
                        Ecal_z = Ecal_z*RAmat_z_mmp;
                        ProdMat2_x = ProdMat2_x*RAmatT_x_mmp;
                        ProdMat2_y = ProdMat2_y*RAmatT_y_mmp;
                        ProdMat2_z = ProdMat2_z*RAmatT_z_mmp;
                     else 
                        Ecal_x = Ecal_x*RAmat_x_mmm; 
                        Ecal_y = Ecal_y*RAmat_y_mmm; 
                        Ecal_z = Ecal_z*RAmat_z_mmm;
                        ProdMat2_x = ProdMat2_x*RAmatT_x_mmm;
                        ProdMat2_y = ProdMat2_y*RAmatT_y_mmm;
                        ProdMat2_z = ProdMat2_z*RAmatT_z_mmm;
                     end   
                end
                Ecal_x = Ecal_x*Rmat_x;
                 Ecal_y = Ecal_y*Rmat_y;
                 Ecal_z = Ecal_z*Rmat_z;
                 ProdMat2_x = ProdMat2_x*Rmat_x;
                 ProdMat2_y = ProdMat2_y*Rmat_y;
                 ProdMat2_z = ProdMat2_z*Rmat_z;
              Rmid_x=Rmat_x^dstmp(m);
              eX(m,n)=Ecal_x*Rmid_x*ProdMat2_x*Smat_x;
              Rmid_y=Rmat_y^dstmp(m);
              eY(m,n)=Ecal_y*Rmid_y*ProdMat2_y*Smat_y;
              Rmid_z=Rmat_z^dstmp(m);
              eZ(m,n)=Ecal_z*Rmid_z*ProdMat2_z*Smat_z;
                 
            end
        end  
            
    end
    E = eX.*eY.*eZ;
    %disp(ePerp)
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
    J = zeros(length(E),4);
    if nargin < 3 
         
        for i = 1:4
        xpert = x;
        protocol.diff=i;   
        xpert(i) = xpert(i)*(1+dx);    
        Epert = Cuboid_MM_SWOGSE_3D(xpert,protocol);
        dEtdx = (Epert - E)/(xpert(i)*dx);
        J(:,i) = dEtdx;
        end       
      
    else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
       
        
        for i = 1:4            
            if x_deriv(i) ~= 0  
              
                xpert = x;
                protocol.diff=i;   
                xpert(i) = xpert(i)*(1+dx);    
                Epert = Cuboid_MM_SWOGSE_3D(xpert,protocol,x_deriv);
                dEtdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEtdx;
            end
        end
         
    end   
    protocol.diff=0;
   
end
    


    
    