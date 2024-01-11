function [E J]=Cuboid_MM_GEN(x,protocol,x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a cuboid/ensemble of cuboids.
% 
% [E,J]=Cuboid_MM_GEN(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate consisting of a cuboid/ 
% an ensemble of cuboids and a diffusion protocol specified in the input
% Substrate: cuboid / ensemble of cuboids
% Diffusion pulse sequence: Generalized waveform (GEN)
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
G=protocol.G;
[M totsizeG]=size(G);
RotN = length(cube_rotation);

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

if isfield(protocol,'smalldel') && isfield(protocol,'delta')
    K = floor((protocol.smalldel+1E-10)./tau)+1;
    dstmp = floor((protocol.delta-1E-10)./protocol.tau)-floor((protocol.smalldel+1E-10)./protocol.tau);
     %disp('Angle method using exp(i(k-n)theta)')
    G_mag=sqrt(G(:,1:3:end).^2+G(:,2:3:end).^2+ G(:,3:3:end).^2);
    indG=G_mag>1E-6;
    G_tmpX=G(:,1:3:end);
    G_tmpY=G(:,2:3:end);
    G_tmpZ=G(:,3:3:end);
    Gx_hat0 = zeros(size(G_mag));
    Gy_hat0 = zeros(size(G_mag));
    Gz_hat0 = zeros(size(G_mag));

    Gx_hat0(indG)=G_tmpX(indG)./G_mag(indG);% unitless projection on x axis
    Gy_hat0(indG)=G_tmpY(indG)./G_mag(indG);% unitless projection on y axis
    Gz_hat0(indG)=G_tmpZ(indG)./G_mag(indG);% unitless projection on z axis
    
       
    eX=zeros(M,RotN); % x direction
    eY = zeros(M,RotN); % y direction
    eZ = zeros(M,RotN); % z direction
   
    Gind_mat = round(G_mag./protocol.gstep);
    if protocol.mirror == 1        
        for n = 1:RotN
            inv_rot = cube_rotation{n}';
            for m=1:M             
               Ecal_x = Smat_x';
               Ecal_y = Smat_y';
               Ecal_z = Smat_z';
                for ind=1:K(m)                
                    Gind = Gind_mat(m,ind);
                    if Gind ~=0
                        rot_gradient = inv_rot*[Gx_hat0(m,ind); Gy_hat0(m,ind); Gz_hat0(m,ind)];
                        Gx_hat  =rot_gradient(1);
                        Gy_hat  =rot_gradient(2);
                        Gz_hat  =rot_gradient(3);

                        AmatU_x = real(Amat_x)+1i*Gx_hat.*imag(Amat_x);
                        AmatUD_x = AmatU_x^Gind;                       
                        Ecal_x=Ecal_x*Rmat_x*AmatUD_x;

                        AmatU_y = real(Amat_y)+1i*Gy_hat.*imag(Amat_y);
                        AmatUD_y = AmatU_y^Gind;                         
                        Ecal_y=Ecal_y*Rmat_y*AmatUD_y;

                        AmatU_z = real(Amat_z)+1i*Gz_hat.*imag(Amat_z);
                        AmatUD_z = AmatU_z^Gind;                         
                        Ecal_z=Ecal_z*Rmat_z*AmatUD_z;                        
                    else
                        Ecal_x=Ecal_x*Rmat_x;                        
                        Ecal_y=Ecal_y*Rmat_y;                       
                        Ecal_z=Ecal_z*Rmat_z;                           
                    end
                end
                Rmid_x=Rmat_x^dstmp(m);
                Rmid_y=Rmat_y^dstmp(m);
                Rmid_z=Rmat_z^dstmp(m);
                eX(m,n)=Ecal_x*Rmid_x*Ecal_x';
                eY(m,n)=Ecal_y*Rmid_y*Ecal_y';
                eZ(m,n)=Ecal_z*Rmid_z*Ecal_z';
            end
                          
        end
    elseif protocol.mirror == 0
        for n = 1:RotN
            inv_rot = cube_rotation{n}';
            for m=1:M             
               Ecal_x = Smat_x';
               Ecal_y = Smat_y';
               Ecal_z = Smat_z';
               ProdMat2_x = eye(size(Rmat_x));
               ProdMat2_y = eye(size(Rmat_x));
               ProdMat2_z = eye(size(Rmat_x));
                for ind=1:K(m)
                          
                    Gind = Gind_mat(m,ind);
                    if Gind
                         rot_gradient = inv_rot*[Gx_hat0(m,ind); Gy_hat0(m,ind); Gz_hat0(m,ind)];
                        Gx_hat  =rot_gradient(1);
                        Gy_hat  =rot_gradient(2);
                        Gz_hat  =rot_gradient(3);         
                        AmatU_x = real(Amat_x)+1i*Gx_hat.*imag(Amat_x);
                        AmatUD_x = AmatU_x^Gind;     
                        AmatUDT_x = AmatUD_x';    
                        ProdMat2_x = ProdMat2_x*Rmat_x*AmatUDT_x;
                        Ecal_x=Ecal_x*Rmat_x*AmatUD_x;

                        AmatU_y = real(Amat_y)+1i*Gy_hat.*imag(Amat_y);
                        AmatUD_y = AmatU_y^Gind;     
                        AmatUDT_y = AmatUD_y';    
                        ProdMat2_y = ProdMat2_y*Rmat_y*AmatUDT_y;
                        Ecal_y=Ecal_y*Rmat_y*AmatUD_y;

                        AmatU_z = real(Amat_z)+1i*Gz_hat.*imag(Amat_z);
                        AmatUD_z = AmatU_z^Gind;     
                        AmatUDT_z = AmatUD_z';    
                        ProdMat2_z = ProdMat2_z*Rmat_z*AmatUDT_z;
                        Ecal_z=Ecal_z*Rmat_z*AmatUD_z;                        
                    else
                        Ecal_x=Ecal_x*Rmat_x;
                        ProdMat2_x = ProdMat2_x*Rmat_x;  
                  
                        Ecal_y=Ecal_y*Rmat_y;
                        ProdMat2_y = ProdMat2_y*Rmat_y;

                        Ecal_z=Ecal_z*Rmat_z;
                        ProdMat2_z = ProdMat2_z*Rmat_z;                    
                    end
                end
                Rmid_x=Rmat_x^dstmp(m);
                Rmid_y=Rmat_y^dstmp(m);
                Rmid_z=Rmat_z^dstmp(m);
                eX(m,n)=Ecal_x*Rmid_x*ProdMat2_x*Smat_x;
                eY(m,n)=Ecal_y*Rmid_y*ProdMat2_y*Smat_y;
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
else
    K=totsizeG/3;
   %disp('Angle method using exp(i(k-n)theta)')
   G_mag=sqrt(G(:,1:3:end).^2+G(:,2:3:end).^2+ G(:,3:3:end).^2);
    indG=G_mag>1E-6;
    G_tmpX=G(:,1:3:end);
    G_tmpY=G(:,2:3:end);
    G_tmpZ=G(:,3:3:end);
    Gx_hat0 = zeros(size(G_mag));
    Gy_hat0 = zeros(size(G_mag));
    Gz_hat0 = zeros(size(G_mag));
    Gx_hat0(indG)=G_tmpX(indG)./G_mag(indG);% unitless projection on x axis
    Gy_hat0(indG)=G_tmpY(indG)./G_mag(indG);% unitless projection on y axis
    Gz_hat0(indG)=G_tmpZ(indG)./G_mag(indG);% unitless projection on z axis
    
       
    eX=zeros(M,RotN); % x direction
    eY = zeros(M,RotN); % y direction
    eZ = zeros(M,RotN); % z direction

    Gind_mat = round(G_mag./protocol.gstep);
    for n = 1:RotN
        inv_rot = cube_rotation{n}';
        for m=1:M             
           Ecal_x = Smat_x';
           Ecal_y = Smat_y';
           Ecal_z = Smat_z';
            for ind=1:K
               
                Gind = Gind_mat(m,ind);
                if Gind ~=0

                    rot_gradient = inv_rot*[Gx_hat0(m,ind); Gy_hat0(m,ind); Gz_hat0(m,ind)];
                    Gx_hat  =rot_gradient(1);
                    Gy_hat  =rot_gradient(2);
                    Gz_hat  =rot_gradient(3);
                    AmatU_x = real(Amat_x)+1i*Gx_hat.*imag(Amat_x);
                    AmatUD_x = AmatU_x^Gind;                        
                    Ecal_x=Ecal_x*Rmat_x*AmatUD_x;

                    AmatU_y = real(Amat_y)+1i*Gy_hat.*imag(Amat_y);
                    AmatUD_y = AmatU_y^Gind;                       
                    Ecal_y=Ecal_y*Rmat_y*AmatUD_y;

                    AmatU_z = real(Amat_z)+1i*Gz_hat.*imag(Amat_z);
                    AmatUD_z = AmatU_z^Gind;                  
                    Ecal_z=Ecal_z*Rmat_z*AmatUD_z;                        
                else
                    Ecal_x=Ecal_x*Rmat_x;
                  
                    Ecal_y=Ecal_y*Rmat_y;
                   
                    Ecal_z=Ecal_z*Rmat_z;                                 
                end
            end
            eX(m,n)=Ecal_x*Rmat_x*Smat_x;
            eY(m,n)=Ecal_y*Rmat_y*Smat_y;
            eZ(m,n)=Ecal_z*Rmat_z*Smat_z;   
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
        Epert = Cuboid_MM_GEN(xpert,protocol);
        dEtdx = (Epert - E)/(xpert(i)*dx);
        J(:,i) = dEtdx;
        end       
      
    else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
     
        
        for i = 1:4            
            if x_deriv(i) ~= 0               
                xpert = x;
                protocol.diff=i;   
                xpert(i) = xpert(i)*(1+dx);    
                Epert = Cuboid_MM_GEN(xpert,protocol,x_deriv);
                dEtdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEtdx;
            end
        end
         
    end   
    protocol.diff=0;
   
end
    


    
    