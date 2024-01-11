function [E J]=GammaFiniteCylinders_MM_SWOGSE_3D(x,protocol,x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a GammaFiniteCylinders compartment.
% 
% [E,J]=GammaFiniteCylinders_MM_SWOGSE_3D(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate consisting of parallel
% finite cylinders with a Gamma distribution of radii and a diffusion 
% protocol specified in the input
% Substrate: Parallel, impermeable finite cylinders with a Gamma
% distribution of radii
% Diffusion pulse sequence: Square wave oscillating gradients with varying
% gradient orientation (SWOGSE_3D)
% Signal approximation: Matrix method (MM) - propagator expressed via 
%       eigenmode expansion
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% x - size 6 vector of model parameters in SI units for GammaFiniteCylinders:
%       x(1) - free diffusivity of the material inside the cylinders.
%       x(2) - mean radius of the distribution.
%       x(3) - shape parameter of the Gamma distribution
%       x(4) - eccentricity (ratio between length and diameter, the same for all sizes)
%       x(5) - polar angle theta in spherical coordinates desbribing the fibre
% direction
%       x(6) - azimuthal angle phi in spherical coordinates describing the
% fibre direction
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

if ~isfield(protocol,'complex')
  protocol.complex='real';
end
% Model parameters
theta = x(4);
phi = x(5);
fibredir = [sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];
tau=protocol.tau;

dstmp = floor((protocol.delta-1E-10)./protocol.tau)-floor(protocol.smalldel./protocol.tau);

cos_theta = cos(theta);
sin_theta = sqrt(1-cos_theta^2);
cos_phi = cos(phi);
sin_phi = sqrt(1-cos_phi^2);


 v = [cos_theta*cos_phi^2+sin_phi^2; -(1-cos_theta)*sin_phi*cos_phi; -sin_theta*cos_phi]; % vectors for the new x and y directions; see Ozarslan 2010
 w = [ -(1-cos_theta)*sin_phi*cos_phi; cos_theta*sin_phi^2+cos_phi^2; -sin_theta*sin_phi];

% Calculating the parallel & perpendicular signal is similar  (parallel plane & cylindrical restriction )
M = size(protocol.Gx',1);
ePerp=zeros(M,length(protocol.A_cyl)); % cylinder restriction
ePar = zeros(M,length(protocol.A_cyl)); % parallel planes 

  for indj = 1:length(protocol.A_cyl)
        % Perpendicular
    if protocol.diff==0
        Amat_cyl=protocol.A_cyl{indj};
        Smat_cyl=protocol.S_cyl{indj};
        Rmat_cyl=protocol.R_cyl{indj};
        Amat_plane=protocol.A_plane{indj};
        Smat_plane=protocol.S_plane{indj};
        Rmat_plane=protocol.R_plane{indj};    
        Rweight = repmat(protocol.Rweight,[M,1]);
    elseif protocol.diff==1
        Amat_cyl=protocol.A_cyl{indj};
        Smat_cyl=protocol.S_cyl{indj};
        Rmat_cyl=protocol.RpertD_cyl{indj};
        Amat_plane=protocol.A_plane{indj};
        Smat_plane=protocol.S_plane{indj};
        Rmat_plane=protocol.RpertD_plane{indj};
        Rweight = repmat(protocol.Rweight,[M,1]);
    elseif protocol.diff==2
        Amat_cyl=protocol.Aperta_cyl{indj};
        Smat_cyl=protocol.Sperta_cyl{indj};
        Rmat_cyl=protocol.Rperta_cyl{indj};
        Amat_plane=protocol.Aperta_plane{indj};
        Smat_plane=protocol.Sperta_plane{indj};
        Rmat_plane=protocol.Rperta_plane{indj};   
        Rweight = repmat(protocol.Rweight_perta,[M,1]);
    elseif protocol.diff==3
        Amat_cyl=protocol.Apertsh_cyl{indj};
        Smat_cyl=protocol.Spertsh_cyl{indj};
        Rmat_cyl=protocol.Rpertsh_cyl{indj};
        Amat_plane=protocol.Apertsh_plane{indj};
        Smat_plane=protocol.Spertsh_plane{indj};
        Rmat_plane=protocol.Rpertsh_plane{indj};     
        Rweight = repmat(protocol.Rweight_pertsh,[M,1]);
    elseif protocol.diff==4
        Amat_cyl=protocol.A_cyl{indj};
        Smat_cyl=protocol.S_cyl{indj};
        Rmat_cyl=protocol.R_cyl{indj};
        Amat_plane=protocol.ApertEcc_plane{indj};
        Smat_plane=protocol.SpertEcc_plane{indj};
        Rmat_plane=protocol.RpertEcc_plane{indj};   
        Rweight = repmat(protocol.Rweight,[M,1]);

    else
        error('protocol.diff is not suitable')
    end
    kDn_cyl=protocol.kDn_cyl{indj};K = floor((protocol.smalldel+1E-10)./tau)+1;

    if protocol.angle == 1 % gradient only along x
        Gx = protocol.Gx';       
        G_mag = Gx;
        G_ind=round(G_mag./protocol.gstep);
        qhat_X=Gx./G_mag; % qhatx  
        if protocol.mirror == 1
            for m=1:M

                time_vec = tau*(0:K(m)-1);
                sign_vec1 = (-1).^floor((protocol.omegax(m)*time_vec-protocol.phix(m))./pi-1E-10);
                Ecal_cyl=Smat_cyl';     
                Ecal_plane=Smat_plane';  

                grad_dir_p = [qhat_X(m); 0; 0];
                Amat_cyl_p=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dir_p,w))+1i*imag(Amat_cyl)*dot(grad_dir_p,v); 
                Amat_cyl_m=Amat_cyl_p'; 
                RAmat_cyl_p = Rmat_cyl*Amat_cyl_p^G_ind(m);
                RAmat_cyl_m = Rmat_cyl*Amat_cyl_m^G_ind(m);  

                Amat_plane_p=real(Amat_plane)+ 1i* imag(Amat_plane)*dot(grad_dir_p,fibredir);
                Amat_plane_m=Amat_plane_p'; 
                RAmat_plane_p = Rmat_plane*Amat_plane_p^G_ind(m);
                RAmat_plane_m = Rmat_plane*Amat_plane_m^G_ind(m);  


                Ecal_cyl = Ecal_cyl*Rmat_cyl;  % first gradient point is 0
                Ecal_plane = Ecal_plane*Rmat_plane;  % first gradient point is 0
                for i = 2:K(m)-1
                     if sign_vec1(i) >0 
                         Ecal_cyl = Ecal_cyl*RAmat_cyl_p;  
                         Ecal_plane = Ecal_plane*RAmat_plane_p;   
                     else
                         Ecal_cyl = Ecal_cyl*RAmat_cyl_m;
                          Ecal_plane = Ecal_plane*RAmat_plane_m;
                     end                        
                end
                Ecal_cyl = Ecal_cyl*Rmat_cyl;            
                Ecal_plane = Ecal_plane*Rmat_plane;

              Rmid_cyl=Rmat_cyl^dstmp(m);
              ePerp(m,indj)=Ecal_cyl*Rmid_cyl*Ecal_cyl';  
              Rmid_plane=Rmat_plane^dstmp(m);
              ePar(m,indj)=Ecal_plane*Rmid_plane*Ecal_plane'; 


            end       
        elseif protocol.mirror == 0
            for m=1:M

                time_vec = tau*(0:K(m)-1);
                sign_vec1 = (-1).^floor((protocol.omegax(m)*time_vec-protocol.phix(m))./pi-1E-10);

                Ecal_plane=Smat_plane';

                Ecal_cyl=Smat_cyl';                   

                grad_dir_p = [qhat_X(m); 0; 0];
                Amat_cyl_p=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dir_p,w))+1i*imag(Amat_cyl)*dot(grad_dir_p,v); 
                Amat_cyl_m=Amat_cyl_p'; 
                RAmat_cyl_p = Rmat_cyl*Amat_cyl_p^G_ind(m);
                RAmat_cyl_m = Rmat_cyl*Amat_cyl_m^G_ind(m);  
                RAmatT_cyl_p = Rmat_cyl*(Amat_cyl_p')^G_ind(m);
                RAmatT_cyl_m = Rmat_cyl*(Amat_cyl_m')^G_ind(m);

                Amat_plane_p= real(Amat_plane)+ 1i* imag(Amat_plane)*dot(grad_dir_p,fibredir);
                Amat_plane_m=Amat_plane_p'; 
                RAmat_plane_p = Rmat_plane*Amat_plane_p^G_ind(m);
                RAmat_plane_m = Rmat_plane*Amat_plane_m^G_ind(m);  
                RAmatT_plane_p = Rmat_plane*(Amat_plane_p')^G_ind(m);
                RAmatT_plane_m = Rmat_plane*(Amat_plane_m')^G_ind(m);

                Ecal_cyl = Ecal_cyl*Rmat_cyl;  % first gradient point is 0
                ProdMat2_cyl = Rmat_cyl; 
                Ecal_plane = Ecal_plane*Rmat_plane;  % first gradient point is 0
                ProdMat2_plane = Rmat_plane; 
                for i=2:K(m)-1
                    if sign_vec1(i) >0 
                         Ecal_cyl = Ecal_cyl*RAmat_cyl_p;
                         ProdMat2_cyl = ProdMat2_cyl*RAmatT_cyl_p; 
                         Ecal_plane = Ecal_plane*RAmat_plane_p;
                        ProdMat2_plane = ProdMat2_plane*RAmatT_plane_p;
                     else
                         Ecal_cyl = Ecal_cyl*RAmat_cyl_m;
                         ProdMat2_cyl = ProdMat2_cyl*RAmatT_cyl_m;
                        Ecal_plane = Ecal_plane*RAmat_plane_m;
                         ProdMat2_plane = ProdMat2_plane*RAmatT_plane_m;
                    end
                end
                Ecal_cyl = Ecal_cyl*Rmat_cyl;            
                Ecal_plane = Ecal_plane*Rmat_plane;
                ProdMat2_cyl = ProdMat2_cyl*Rmat_cyl;
                ProdMat2_plane = ProdMat2_plane*Rmat_plane;
              Rmid_cyl=Rmat_cyl^dstmp(m);
              ePerp(m,indj)=Ecal_cyl*Rmid_cyl*ProdMat2_cyl*Smat_cyl;
              Rmid_plane=Rmat_plane^dstmp(m);
              ePar(m,indj)=Ecal_plane*Rmid_plane*ProdMat2_plane*Smat_plane;  

            end       
        end
    
    elseif protocol.angle == 2
        Gx = protocol.Gx';
        Gy = protocol.Gy';
      
        G_mag = sqrt(Gx.^2+Gy.^2);
        G_ind=round(G_mag./protocol.gstep);
        qhat_X=Gx./G_mag; % qhatx  
        qhat_Y=Gy./G_mag; % qhaty 
        if protocol.mirror == 1
            for m=1:M

               time_vec = tau*(0:K(m)-1);
                sign_vec1 = (-1).^floor((protocol.omegax(m)*time_vec-protocol.phix(m))./pi-1E-10);
                sign_vec2 = (-1).^floor((protocol.omegay(m)*time_vec-protocol.phiy(m))./pi-1E-10);               

                % Calculating perpendicular
                Ecal_cyl=Smat_cyl'; 
                Ecal_plane=Smat_plane';           

                grad_dir_pp = [qhat_X(m); qhat_Y(m); 0];
                grad_dir_pm = [qhat_X(m); -qhat_Y(m); 0];
                grad_dir_mp = [-qhat_X(m); qhat_Y(m); 0];

                Amat_cyl_pp=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dir_pp,w))+1i*imag(Amat_cyl)*dot(grad_dir_pp,v);
                Amat_cyl_pm=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dir_pm,w))+1i*imag(Amat_cyl)*dot(grad_dir_pm,v); 
                Amat_cyl_mp=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dir_mp,w))+1i*imag(Amat_cyl)*dot(grad_dir_mp,v); 
                Amat_cyl_mm=Amat_cyl_pp';

                Amat_plane_pp= real(Amat_plane)+ 1i* imag(Amat_plane)*dot(grad_dir_pp,fibredir);
                Amat_plane_pm= real(Amat_plane)+ 1i* imag(Amat_plane)*dot(grad_dir_pm,fibredir);
                Amat_plane_mp= real(Amat_plane)+ 1i* imag(Amat_plane)*dot(grad_dir_mp,fibredir);
                Amat_plane_mm = Amat_plane_pp';

                 RAmat_cyl_pp = Rmat_cyl*Amat_cyl_pp^G_ind(m);
                 RAmat_cyl_pm = Rmat_cyl*Amat_cyl_pm^G_ind(m);
                 RAmat_cyl_mp = Rmat_cyl*Amat_cyl_mp^G_ind(m);
                 RAmat_cyl_mm = Rmat_cyl*Amat_cyl_mm^G_ind(m);  

                 RAmat_plane_pp = Rmat_plane*Amat_plane_pp^G_ind(m);
                 RAmat_plane_pm = Rmat_plane*Amat_plane_pm^G_ind(m);
                 RAmat_plane_mp = Rmat_plane*Amat_plane_mp^G_ind(m);
                 RAmat_plane_mm = Rmat_plane*Amat_plane_mm^G_ind(m);

                Ecal_cyl = Ecal_cyl*Rmat_cyl;  % first gradient point is 0
                Ecal_plane = Ecal_plane*Rmat_plane; 
                for i = 2:K(m)-1
                     if sign_vec1(i) > 0 && sign_vec2(i) > 0
                        Ecal_cyl = Ecal_cyl*RAmat_cyl_pp;
                        Ecal_plane = Ecal_plane*RAmat_plane_pp;
                     elseif sign_vec1(i) > 0 && sign_vec2(i) < 0
                         Ecal_cyl = Ecal_cyl*RAmat_cyl_pm;
                        Ecal_plane = Ecal_plane*RAmat_plane_pm;
                     elseif sign_vec1(i) < 0 && sign_vec2(i) > 0
                         Ecal_cyl = Ecal_cyl*RAmat_cyl_mp;
                         Ecal_plane = Ecal_plane*RAmat_plane_mp;
                     else
                         Ecal_cyl = Ecal_cyl*RAmat_cyl_mm;
                        Ecal_plane = Ecal_plane*RAmat_plane_mm;
                     end                                 
                end
                Ecal_cyl = Ecal_cyl*Rmat_cyl;            
                Ecal_plane = Ecal_plane*Rmat_plane;

              Rmid_cyl=Rmat_cyl^dstmp(m);
              ePerp(m,indj)=Ecal_cyl*Rmid_cyl*Ecal_cyl';  
              Rmid_plane=Rmat_plane^dstmp(m);
              ePar(m,indj)=Ecal_plane*Rmid_plane*Ecal_plane'; 

            end      
        elseif protocol.mirror == 0
            for m=1:M
                time_vec = tau*(0:K(m)-1);
                sign_vec1 = (-1).^floor((protocol.omegax(m)*time_vec-protocol.phix(m))./pi-1E-10);
                sign_vec2 = (-1).^floor((protocol.omegay(m)*time_vec-protocol.phiy(m))./pi-1E-10);            

                % Calculating perpendicular
                Ecal_cyl=Smat_cyl';       
                Ecal_plane = Smat_plane';           

                grad_dir_pp = [qhat_X(m); qhat_Y(m); 0];
                grad_dir_pm = [qhat_X(m); -qhat_Y(m); 0];
                grad_dir_mp = [-qhat_X(m); qhat_Y(m); 0];

                 Amat_cyl_pp=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dir_pp,w))+1i*imag(Amat_cyl)*dot(grad_dir_pp,v);
                Amat_cyl_pm=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dir_pm,w))+1i*imag(Amat_cyl)*dot(grad_dir_pm,v); 
                Amat_cyl_mp=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dir_mp,w))+1i*imag(Amat_cyl)*dot(grad_dir_mp,v); 
                Amat_cyl_mm=Amat_cyl_pp';

                 Amat_plane_pp= real(Amat_plane)+ 1i* imag(Amat_plane)*dot(grad_dir_pp,fibredir);
                Amat_plane_pm= real(Amat_plane)+ 1i* imag(Amat_plane)*dot(grad_dir_pm,fibredir);
                Amat_plane_mp= real(Amat_plane)+ 1i* imag(Amat_plane)*dot(grad_dir_mp,fibredir);
                Amat_plane_mm = Amat_plane_pp';

                 RAmat_cyl_pp = Rmat_cyl*Amat_cyl_pp^G_ind(m);
                 RAmat_cyl_pm = Rmat_cyl*Amat_cyl_pm^G_ind(m);
                 RAmat_cyl_mp = Rmat_cyl*Amat_cyl_mp^G_ind(m);
                 RAmat_cyl_mm = Rmat_cyl*Amat_cyl_mm^G_ind(m);

                RAmat_plane_pp = Rmat_plane*Amat_plane_pp^G_ind(m);
                 RAmat_plane_pm = Rmat_plane*Amat_plane_pm^G_ind(m);
                 RAmat_plane_mp = Rmat_plane*Amat_plane_mp^G_ind(m);
                 RAmat_plane_mm = Rmat_plane*Amat_plane_mm^G_ind(m);

                 RAmatT_cyl_pp = Rmat_cyl*(Amat_cyl_pp')^G_ind(m);
                RAmatT_cyl_pm = Rmat_cyl*(Amat_cyl_pm')^G_ind(m);
                RAmatT_cyl_mp = Rmat_cyl*(Amat_cyl_mp')^G_ind(m);
                RAmatT_cyl_mm = Rmat_cyl*(Amat_cyl_mm')^G_ind(m);

                RAmatT_plane_pp = Rmat_plane*(Amat_plane_pp')^G_ind(m);
                RAmatT_plane_pm = Rmat_plane*(Amat_plane_pm')^G_ind(m);
                RAmatT_plane_mp = Rmat_plane*(Amat_plane_mp')^G_ind(m);
                RAmatT_plane_mm = Rmat_plane*(Amat_plane_mm')^G_ind(m);

                Ecal_cyl = Ecal_cyl*Rmat_cyl;  % first gradient point is 0
                ProdMat2_cyl = Rmat_cyl; 
                Ecal_plane = Ecal_plane*Rmat_plane;  % first gradient point is 0
                ProdMat2_plane = Rmat_plane; 
                for i=2:K(m)-1                    
                    if sign_vec1(i) >0 && sign_vec2(i) >0
                         Ecal_cyl = Ecal_cyl*RAmat_cyl_pp;
                          ProdMat2_cyl = ProdMat2_cyl*RAmatT_cyl_pp;
                          Ecal_plane = Ecal_plane*RAmat_plane_pp;
                          ProdMat2_plane = ProdMat2_plane*RAmatT_plane_pp;
                     elseif sign_vec1(i) >0 && sign_vec2(i) <0
                         Ecal_cyl = Ecal_cyl*RAmat_cyl_pm;
                         ProdMat2_cyl = ProdMat2_cyl*RAmatT_cyl_pm;   
                         Ecal_plane = Ecal_plane*RAmat_plane_pm;
                         ProdMat2_plane = ProdMat2_plane*RAmatT_plane_pm; 
                     elseif sign_vec1(i) <0 && sign_vec2(i) >0
                         Ecal_cyl = Ecal_cyl*RAmat_cyl_mp;
                         ProdMat2_cyl = ProdMat2_cyl*RAmatT_cyl_mp;   
                         Ecal_plane = Ecal_plane*RAmat_plane_mp;
                         ProdMat2_plane = ProdMat2_plane*RAmatT_plane_mp;   
                     else
                         Ecal_cyl = Ecal_cyl*RAmat_cyl_mm;
                         ProdMat2_cyl = ProdMat2_cyl*RAmatT_cyl_mm;
                        Ecal_plane = Ecal_plane*RAmat_plane_mm;
                         ProdMat2_plane = ProdMat2_plane*RAmatT_plane_mm;
                     end                       
                end
                Ecal_cyl = Ecal_cyl*Rmat_cyl;            
                Ecal_plane = Ecal_plane*Rmat_plane;
                ProdMat2_cyl = ProdMat2_cyl*Rmat_cyl;
                ProdMat2_plane = ProdMat2_plane*Rmat_plane;
              Rmid_cyl=Rmat_cyl^dstmp(m);
              ePerp(m,indj)=Ecal_cyl*Rmid_cyl*ProdMat2_cyl*Smat_cyl;  
              Rmid_plane=Rmat_plane^dstmp(m);
              ePar(m,indj)=Ecal_plane*Rmid_plane*ProdMat2_plane*Smat_plane;                 
            end

        end
     
    elseif protocol.angle == 3
        Gx = protocol.Gx';
        Gz = protocol.Gz';
      
        G_mag = sqrt(Gx.^2+Gz.^2);
        G_ind=round(G_mag./protocol.gstep);
        qhat_X=Gx./G_mag; % qhatx  
        qhat_Z=Gz./G_mag; % qhatz
        if protocol.mirror == 1
            for m=1:M

                time_vec = tau*(0:K(m)-1);
                sign_vec1 = (-1).^floor((protocol.omegax(m)*time_vec-protocol.phix(m))./pi-1E-10);
                sign_vec2 = (-1).^floor((protocol.omegaz(m)*time_vec-protocol.phiz(m))./pi-1E-10);      

                % Calculating perpendicular
                Ecal_cyl=Smat_cyl'; 
                Ecal_plane=Smat_plane';  

                grad_dir_pp = [qhat_X(m); 0; qhat_Z(m); ];
                grad_dir_pm = [qhat_X(m); 0;-qhat_Z(m); ];
                grad_dir_mp = [-qhat_X(m); 0; qhat_Z(m);];

                Amat_cyl_pp=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dir_pp,w))+1i*imag(Amat_cyl)*dot(grad_dir_pp,v);
                Amat_cyl_pm=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dir_pm,w))+1i*imag(Amat_cyl)*dot(grad_dir_pm,v); 
                Amat_cyl_mp=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dir_mp,w))+1i*imag(Amat_cyl)*dot(grad_dir_mp,v); 
                Amat_cyl_mm=Amat_cyl_pp';

                Amat_plane_pp= real(Amat_plane)+ 1i* imag(Amat_plane)*dot(grad_dir_pp,fibredir);
                Amat_plane_pm= real(Amat_plane)+ 1i* imag(Amat_plane)*dot(grad_dir_pm,fibredir);
                Amat_plane_mp= real(Amat_plane)+ 1i* imag(Amat_plane)*dot(grad_dir_mp,fibredir);
                Amat_plane_mm = Amat_plane_pp';

                 RAmat_cyl_pp = Rmat_cyl*Amat_cyl_pp^G_ind(m);
                 RAmat_cyl_pm = Rmat_cyl*Amat_cyl_pm^G_ind(m);
                 RAmat_cyl_mp = Rmat_cyl*Amat_cyl_mp^G_ind(m);
                 RAmat_cyl_mm = Rmat_cyl*Amat_cyl_mm^G_ind(m);  

                 RAmat_plane_pp = Rmat_plane*Amat_plane_pp^G_ind(m);
                 RAmat_plane_pm = Rmat_plane*Amat_plane_pm^G_ind(m);
                 RAmat_plane_mp = Rmat_plane*Amat_plane_mp^G_ind(m);
                 RAmat_plane_mm = Rmat_plane*Amat_plane_mm^G_ind(m);

                Ecal_cyl = Ecal_cyl*Rmat_cyl;  % first gradient point is 0
                Ecal_plane = Ecal_plane*Rmat_plane; 
                for i = 2:K(m)-1
                     if sign_vec1(i) > 0 && sign_vec2(i) > 0
                        Ecal_cyl = Ecal_cyl*RAmat_cyl_pp;
                        Ecal_plane = Ecal_plane*RAmat_plane_pp;
                     elseif sign_vec1(i) > 0 && sign_vec2(i) < 0
                         Ecal_cyl = Ecal_cyl*RAmat_cyl_pm;
                        Ecal_plane = Ecal_plane*RAmat_plane_pm;
                     elseif sign_vec1(i) < 0 && sign_vec2(i) > 0
                         Ecal_cyl = Ecal_cyl*RAmat_cyl_mp;
                         Ecal_plane = Ecal_plane*RAmat_plane_mp;
                     else
                         Ecal_cyl = Ecal_cyl*RAmat_cyl_mm;
                        Ecal_plane = Ecal_plane*RAmat_plane_mm;
                     end                             
                end
                Ecal_cyl = Ecal_cyl*Rmat_cyl;            
                Ecal_plane = Ecal_plane*Rmat_plane;

              Rmid_cyl=Rmat_cyl^dstmp(m);
              ePerp(m,indj)=Ecal_cyl*Rmid_cyl*Ecal_cyl';  
              Rmid_plane=Rmat_plane^dstmp(m);
              ePar(m,indj)=Ecal_plane*Rmid_plane*Ecal_plane'; 

            end      
        elseif protocol.mirror == 0
            for m=1:M

               time_vec = tau*(0:K(m)-1);
                sign_vec1 = (-1).^floor((protocol.omegax(m)*time_vec-protocol.phix(m))./pi-1E-10);
                sign_vec2 = (-1).^floor((protocol.omegaz(m)*time_vec-protocol.phiz(m))./pi-1E-10);
                % Calculating perpendicular
                Ecal_cyl=Smat_cyl';       
                Ecal_plane = Smat_plane';

                 grad_dir_pp = [qhat_X(m); 0; qhat_Z(m); ];
                grad_dir_pm = [qhat_X(m); 0;-qhat_Z(m); ];
                grad_dir_mp = [-qhat_X(m); 0; qhat_Z(m);];

                 Amat_cyl_pp=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dir_pp,w))+1i*imag(Amat_cyl)*dot(grad_dir_pp,v);
                Amat_cyl_pm=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dir_pm,w))+1i*imag(Amat_cyl)*dot(grad_dir_pm,v); 
                Amat_cyl_mp=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dir_mp,w))+1i*imag(Amat_cyl)*dot(grad_dir_mp,v); 
                Amat_cyl_mm=Amat_cyl_pp';

                Amat_plane_pp= real(Amat_plane)+ 1i* imag(Amat_plane)*dot(grad_dir_pp,fibredir);
                Amat_plane_pm= real(Amat_plane)+ 1i* imag(Amat_plane)*dot(grad_dir_pm,fibredir);
                Amat_plane_mp= real(Amat_plane)+ 1i* imag(Amat_plane)*dot(grad_dir_mp,fibredir);
                Amat_plane_mm = Amat_plane_pp';

                 RAmat_cyl_pp = Rmat_cyl*Amat_cyl_pp^G_ind(m);
                 RAmat_cyl_pm = Rmat_cyl*Amat_cyl_pm^G_ind(m);
                 RAmat_cyl_mp = Rmat_cyl*Amat_cyl_mp^G_ind(m);
                 RAmat_cyl_mm = Rmat_cyl*Amat_cyl_mm^G_ind(m);

                RAmat_plane_pp = Rmat_plane*Amat_plane_pp^G_ind(m);
                 RAmat_plane_pm = Rmat_plane*Amat_plane_pm^G_ind(m);
                 RAmat_plane_mp = Rmat_plane*Amat_plane_mp^G_ind(m);
                 RAmat_plane_mm = Rmat_plane*Amat_plane_mm^G_ind(m);

                RAmatT_cyl_pp = Rmat_cyl*(Amat_cyl_pp')^G_ind(m);
                RAmatT_cyl_pm = Rmat_cyl*(Amat_cyl_pm')^G_ind(m);
                RAmatT_cyl_mp = Rmat_cyl*(Amat_cyl_mp')^G_ind(m);
                RAmatT_cyl_mm = Rmat_cyl*(Amat_cyl_mm')^G_ind(m);

                RAmatT_plane_pp = Rmat_plane*(Amat_plane_pp')^G_ind(m);
                RAmatT_plane_pm = Rmat_plane*(Amat_plane_pm')^G_ind(m);
                RAmatT_plane_mp = Rmat_plane*(Amat_plane_mp')^G_ind(m);
                RAmatT_plane_mm = Rmat_plane*(Amat_plane_mm')^G_ind(m);

                Ecal_cyl = Ecal_cyl*Rmat_cyl;  % first gradient point is 0
                ProdMat2_cyl = Rmat_cyl; 
                Ecal_plane = Ecal_plane*Rmat_plane;  % first gradient point is 0
                ProdMat2_plane = Rmat_plane; 
                for i=2:K(m)-1                    
                      if sign_vec1(i) >0 && sign_vec2(i) >0
                         Ecal_cyl = Ecal_cyl*RAmat_cyl_pp;
                          ProdMat2_cyl = ProdMat2_cyl*RAmatT_cyl_pp;
                          Ecal_plane = Ecal_plane*RAmat_plane_pp;
                          ProdMat2_plane = ProdMat2_plane*RAmatT_plane_pp;
                     elseif sign_vec1(i) >0 && sign_vec2(i) <0
                         Ecal_cyl = Ecal_cyl*RAmat_cyl_pm;
                         ProdMat2_cyl = ProdMat2_cyl*RAmatT_cyl_pm;   
                         Ecal_plane = Ecal_plane*RAmat_plane_pm;
                         ProdMat2_plane = ProdMat2_plane*RAmatT_plane_pm; 
                     elseif sign_vec1(i) <0 && sign_vec2(i) >0
                         Ecal_cyl = Ecal_cyl*RAmat_cyl_mp;
                         ProdMat2_cyl = ProdMat2_cyl*RAmatT_cyl_mp;   
                         Ecal_plane = Ecal_plane*RAmat_plane_mp;
                         ProdMat2_plane = ProdMat2_plane*RAmatT_plane_mp;   
                     else
                         Ecal_cyl = Ecal_cyl*RAmat_cyl_mm;
                         ProdMat2_cyl = ProdMat2_cyl*RAmatT_cyl_mm;
                        Ecal_plane = Ecal_plane*RAmat_plane_mm;
                         ProdMat2_plane = ProdMat2_plane*RAmatT_plane_mm;
                     end                               
                end
                Ecal_cyl = Ecal_cyl*Rmat_cyl;            
                Ecal_plane = Ecal_plane*Rmat_plane;
                ProdMat2_cyl = ProdMat2_cyl*Rmat_cyl;
                ProdMat2_plane = ProdMat2_plane*Rmat_plane;
              Rmid_cyl=Rmat_cyl^dstmp(m);
              ePerp(m,indj)=Ecal_cyl*Rmid_cyl*ProdMat2_cyl*Smat_cyl;  
              Rmid_plane=Rmat_plane^dstmp(m);
              ePar(m,indj)=Ecal_plane*Rmid_plane*ProdMat2_plane*Smat_plane; 

            end

        end
 
    elseif protocol.angle == 4
         Gx = protocol.Gx';
         Gy = protocol.Gy';
         Gz = protocol.Gz';
     
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

                Ecal_cyl=Smat_cyl';   
                Ecal_plane=Smat_plane'; 

                grad_dir_ppp = [qhat_X(m); qhat_Y(m); qhat_Z(m); ];
                grad_dir_ppm = [qhat_X(m); qhat_Y(m); -qhat_Z(m); ];
                grad_dir_pmp = [qhat_X(m); -qhat_Y(m); qhat_Z(m); ];
                grad_dir_pmm = [qhat_X(m); -qhat_Y(m); -qhat_Z(m); ];
                grad_dir_mpp = [-qhat_X(m); qhat_Y(m); qhat_Z(m); ];
                grad_dir_mpm = [-qhat_X(m); qhat_Y(m); -qhat_Z(m); ];
                grad_dir_mmp = [-qhat_X(m); -qhat_Y(m); qhat_Z(m); ];
                grad_dir_mmm = [-qhat_X(m); -qhat_Y(m); -qhat_Z(m); ];


                Amat_cyl_ppp=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dir_ppp,w))+1i*imag(Amat_cyl)*dot(grad_dir_ppp,v);
                Amat_cyl_ppm=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dir_ppm,w))+1i*imag(Amat_cyl)*dot(grad_dir_ppm,v); 
                Amat_cyl_pmp=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dir_pmp,w))+1i*imag(Amat_cyl)*dot(grad_dir_pmp,v); 
                Amat_cyl_pmm=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dir_pmm,w))+1i*imag(Amat_cyl)*dot(grad_dir_pmm,v);
                Amat_cyl_mpp=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dir_mpp,w))+1i*imag(Amat_cyl)*dot(grad_dir_mpp,v);
                Amat_cyl_mpm=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dir_mpm,w))+1i*imag(Amat_cyl)*dot(grad_dir_mpm,v); 
                Amat_cyl_mmp=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dir_mmp,w))+1i*imag(Amat_cyl)*dot(grad_dir_mmp,v); 
                Amat_cyl_mmm=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dir_mmm,w))+1i*imag(Amat_cyl)*dot(grad_dir_mmm,v);

                Amat_plane_ppp= real(Amat_plane)+ 1i* imag(Amat_plane)*dot(grad_dir_ppp,fibredir);
                Amat_plane_ppm= real(Amat_plane)+ 1i* imag(Amat_plane)*dot(grad_dir_ppm,fibredir);
                Amat_plane_pmp= real(Amat_plane)+ 1i* imag(Amat_plane)*dot(grad_dir_pmp,fibredir);
                Amat_plane_pmm= real(Amat_plane)+ 1i* imag(Amat_plane)*dot(grad_dir_pmm,fibredir);
                Amat_plane_mpp= real(Amat_plane)+ 1i* imag(Amat_plane)*dot(grad_dir_mpp,fibredir);
                Amat_plane_mpm= real(Amat_plane)+ 1i* imag(Amat_plane)*dot(grad_dir_mpm,fibredir);
                Amat_plane_mmp= real(Amat_plane)+ 1i* imag(Amat_plane)*dot(grad_dir_mmp,fibredir);
                Amat_plane_mmm= real(Amat_plane)+ 1i* imag(Amat_plane)*dot(grad_dir_mmm,fibredir);

                 RAmat_cyl_ppp = Rmat_cyl*Amat_cyl_ppp^G_ind(m);
                 RAmat_cyl_ppm = Rmat_cyl*Amat_cyl_ppm^G_ind(m);
                 RAmat_cyl_pmp = Rmat_cyl*Amat_cyl_pmp^G_ind(m);
                 RAmat_cyl_pmm = Rmat_cyl*Amat_cyl_pmm^G_ind(m); 
                 RAmat_cyl_mpp = Rmat_cyl*Amat_cyl_mpp^G_ind(m);
                 RAmat_cyl_mpm = Rmat_cyl*Amat_cyl_mpm^G_ind(m);
                 RAmat_cyl_mmp = Rmat_cyl*Amat_cyl_mmp^G_ind(m);
                 RAmat_cyl_mmm = Rmat_cyl*Amat_cyl_mmm^G_ind(m);

                 RAmat_plane_ppp = Rmat_plane*Amat_plane_ppp^G_ind(m);
                 RAmat_plane_ppm = Rmat_plane*Amat_plane_ppm^G_ind(m);
                 RAmat_plane_pmp = Rmat_plane*Amat_plane_pmp^G_ind(m);
                 RAmat_plane_pmm = Rmat_plane*Amat_plane_pmm^G_ind(m); 
                 RAmat_plane_mpp = Rmat_plane*Amat_plane_mpp^G_ind(m);
                 RAmat_plane_mpm = Rmat_plane*Amat_plane_mpm^G_ind(m);
                 RAmat_plane_mmp = Rmat_plane*Amat_plane_mmp^G_ind(m);
                 RAmat_plane_mmm = Rmat_plane*Amat_plane_mmm^G_ind(m);

                Ecal_cyl = Ecal_cyl*Rmat_cyl;  % first gradient point is 0
                Ecal_plane = Ecal_plane*Rmat_plane;  % first gradient point is 0
                for i = 2:K(m)-1
                     if sign_vec1(i) >0 && sign_vec2(i) >0 && sign_vec3(i) >0
                        Ecal_cyl = Ecal_cyl*RAmat_cyl_ppp;
                        Ecal_plane = Ecal_plane*RAmat_plane_ppp;
                     elseif sign_vec1(i) >0 && sign_vec2(i) >0 && sign_vec3(i) <0
                         Ecal_cyl = Ecal_cyl*RAmat_cyl_ppm;
                        Ecal_plane = Ecal_plane*RAmat_plane_ppm;
                     elseif sign_vec1(i) >0 && sign_vec2(i) <0 && sign_vec3(i) >0
                         Ecal_cyl = Ecal_cyl*RAmat_cyl_pmp;
                        Ecal_plane = Ecal_plane*RAmat_plane_pmp;
                     elseif sign_vec1(i) >0 && sign_vec2(i) <0 && sign_vec3(i) <0
                         Ecal_cyl = Ecal_cyl*RAmat_cyl_pmm;
                        Ecal_plane = Ecal_plane*RAmat_plane_pmm;
                     elseif sign_vec1(i) <0 && sign_vec2(i) >0 && sign_vec3(i) >0
                        Ecal_cyl = Ecal_cyl*RAmat_cyl_mpp;
                        Ecal_plane = Ecal_plane*RAmat_plane_mpp;
                     elseif sign_vec1(i) <0 && sign_vec2(i) >0 && sign_vec3(i) <0
                         Ecal_cyl = Ecal_cyl*RAmat_cyl_mpm;
                        Ecal_plane = Ecal_plane*RAmat_plane_mpm;
                     elseif sign_vec1(i) <0 && sign_vec2(i) <0 && sign_vec3(i) >0
                         Ecal_cyl = Ecal_cyl*RAmat_cyl_mmp;
                        Ecal_plane = Ecal_plane*RAmat_plane_mmp;
                     else 
                         Ecal_cyl = Ecal_cyl*RAmat_cyl_mmm; 
                        Ecal_plane = Ecal_plane*RAmat_plane_mmm; 
                     end                             
                end
                Ecal_cyl = Ecal_cyl*Rmat_cyl;            
                Ecal_plane = Ecal_plane*Rmat_plane;

              Rmid_cyl=Rmat_cyl^dstmp(m);
              ePerp(m,indj)=Ecal_cyl*Rmid_cyl*Ecal_cyl';    
              Rmid_plane=Rmat_plane^dstmp(m);
              ePar(m,indj)=Ecal_plane*Rmid_plane*Ecal_plane'; 

            end
         elseif protocol.mirror == 0
            for m=1:M

                time_vec = tau*(0:K(m)-1);
                sign_vec1 = (-1).^floor((protocol.omegax(m)*time_vec-protocol.phix(m))./pi-1E-10);
                sign_vec2 = (-1).^floor((protocol.omegay(m)*time_vec-protocol.phiy(m))./pi-1E-10);
                sign_vec3 = (-1).^floor((protocol.omegaz(m)*time_vec-protocol.phiz(m))./pi-1E-10);             

                Ecal_cyl=Smat_cyl';   
                Ecal_plane=Smat_plane'; 

                grad_dir_ppp = [qhat_X(m); qhat_Y(m); qhat_Z(m); ];
                grad_dir_ppm = [qhat_X(m); qhat_Y(m); -qhat_Z(m); ];
                grad_dir_pmp = [qhat_X(m); -qhat_Y(m); qhat_Z(m); ];
                grad_dir_pmm = [qhat_X(m); -qhat_Y(m); -qhat_Z(m); ];
                grad_dir_mpp = [-qhat_X(m); qhat_Y(m); qhat_Z(m); ];
                grad_dir_mpm = [-qhat_X(m); qhat_Y(m); -qhat_Z(m); ];
                grad_dir_mmp = [-qhat_X(m); -qhat_Y(m); qhat_Z(m); ];
                grad_dir_mmm = [-qhat_X(m); -qhat_Y(m); -qhat_Z(m); ];


                Amat_cyl_ppp=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dir_ppp,w))+1i*imag(Amat_cyl)*dot(grad_dir_ppp,v);
                Amat_cyl_ppm=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dir_ppm,w))+1i*imag(Amat_cyl)*dot(grad_dir_ppm,v); 
                Amat_cyl_pmp=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dir_pmp,w))+1i*imag(Amat_cyl)*dot(grad_dir_pmp,v); 
                Amat_cyl_pmm=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dir_pmm,w))+1i*imag(Amat_cyl)*dot(grad_dir_pmm,v);
                Amat_cyl_mpp=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dir_mpp,w))+1i*imag(Amat_cyl)*dot(grad_dir_mpp,v);
                Amat_cyl_mpm=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dir_mpm,w))+1i*imag(Amat_cyl)*dot(grad_dir_mpm,v); 
                Amat_cyl_mmp=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dir_mmp,w))+1i*imag(Amat_cyl)*dot(grad_dir_mmp,v); 
                Amat_cyl_mmm=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dir_mmm,w))+1i*imag(Amat_cyl)*dot(grad_dir_mmm,v);

                Amat_plane_ppp= real(Amat_plane)+ 1i* imag(Amat_plane)*dot(grad_dir_ppp,fibredir);
                Amat_plane_ppm= real(Amat_plane)+ 1i* imag(Amat_plane)*dot(grad_dir_ppm,fibredir);
                Amat_plane_pmp= real(Amat_plane)+ 1i* imag(Amat_plane)*dot(grad_dir_pmp,fibredir);
                Amat_plane_pmm= real(Amat_plane)+ 1i* imag(Amat_plane)*dot(grad_dir_pmm,fibredir);
                Amat_plane_mpp= real(Amat_plane)+ 1i* imag(Amat_plane)*dot(grad_dir_mpp,fibredir);
                Amat_plane_mpm= real(Amat_plane)+ 1i* imag(Amat_plane)*dot(grad_dir_mpm,fibredir);
                Amat_plane_mmp= real(Amat_plane)+ 1i* imag(Amat_plane)*dot(grad_dir_mmp,fibredir);
                Amat_plane_mmm= real(Amat_plane)+ 1i* imag(Amat_plane)*dot(grad_dir_mmm,fibredir);

                 RAmat_cyl_ppp = Rmat_cyl*Amat_cyl_ppp^G_ind(m);
                 RAmat_cyl_ppm = Rmat_cyl*Amat_cyl_ppm^G_ind(m);
                 RAmat_cyl_pmp = Rmat_cyl*Amat_cyl_pmp^G_ind(m);
                 RAmat_cyl_pmm = Rmat_cyl*Amat_cyl_pmm^G_ind(m); 
                 RAmat_cyl_mpp = Rmat_cyl*Amat_cyl_mpp^G_ind(m);
                 RAmat_cyl_mpm = Rmat_cyl*Amat_cyl_mpm^G_ind(m);
                 RAmat_cyl_mmp = Rmat_cyl*Amat_cyl_mmp^G_ind(m);
                 RAmat_cyl_mmm = Rmat_cyl*Amat_cyl_mmm^G_ind(m);

                 RAmat_plane_ppp = Rmat_plane*Amat_plane_ppp^G_ind(m);
                 RAmat_plane_ppm = Rmat_plane*Amat_plane_ppm^G_ind(m);
                 RAmat_plane_pmp = Rmat_plane*Amat_plane_pmp^G_ind(m);
                 RAmat_plane_pmm = Rmat_plane*Amat_plane_pmm^G_ind(m); 
                 RAmat_plane_mpp = Rmat_plane*Amat_plane_mpp^G_ind(m);
                 RAmat_plane_mpm = Rmat_plane*Amat_plane_mpm^G_ind(m);
                 RAmat_plane_mmp = Rmat_plane*Amat_plane_mmp^G_ind(m);
                 RAmat_plane_mmm = Rmat_plane*Amat_plane_mmm^G_ind(m);

                RAmatT_cyl_ppp =  Rmat_cyl*(Amat_cyl_ppp')^G_ind(m);
                RAmatT_cyl_ppm =  Rmat_cyl*(Amat_cyl_ppm')^G_ind(m);
                RAmatT_cyl_pmp =  Rmat_cyl*(Amat_cyl_pmp')^G_ind(m);
                RAmatT_cyl_pmm =  Rmat_cyl*(Amat_cyl_pmm')^G_ind(m);
                RAmatT_cyl_mpp =  Rmat_cyl*(Amat_cyl_mpp')^G_ind(m);
                RAmatT_cyl_mpm =  Rmat_cyl*(Amat_cyl_mpm')^G_ind(m);
                RAmatT_cyl_mmp =  Rmat_cyl*(Amat_cyl_mmp')^G_ind(m);
                RAmatT_cyl_mmm =  Rmat_cyl*(Amat_cyl_mmm')^G_ind(m);

                RAmatT_plane_ppp = Rmat_plane*(Amat_plane_ppp')^G_ind(m);
                RAmatT_plane_ppm = Rmat_plane*(Amat_plane_ppm')^G_ind(m);
                RAmatT_plane_pmp = Rmat_plane*(Amat_plane_pmp')^G_ind(m);
                RAmatT_plane_pmm = Rmat_plane*(Amat_plane_pmm')^G_ind(m);
                RAmatT_plane_mpp = Rmat_plane*(Amat_plane_mpp')^G_ind(m);
                RAmatT_plane_mpm = Rmat_plane*(Amat_plane_mpm')^G_ind(m);
                RAmatT_plane_mmp = Rmat_plane*(Amat_plane_mmp')^G_ind(m);
                RAmatT_plane_mmm = Rmat_plane*(Amat_plane_mmm')^G_ind(m);

                Ecal_cyl = Ecal_cyl*Rmat_cyl;  % first gradient point is 0
                ProdMat2_cyl = Rmat_cyl;
                Ecal_plane = Ecal_plane*Rmat_plane;  % first gradient point is 0
                ProdMat2_plane = Rmat_plane;
                for i = 2:K(m)-1
                     if sign_vec1(i) >0 && sign_vec2(i) >0 && sign_vec3(i) >0
                        Ecal_cyl = Ecal_cyl*RAmat_cyl_ppp;
                        ProdMat2_cyl = ProdMat2_cyl*RAmatT_cyl_ppp;
                        Ecal_plane = Ecal_plane*RAmat_plane_ppp;
                        ProdMat2_plane = ProdMat2_plane*RAmatT_plane_ppp;
                     elseif sign_vec1(i) >0 && sign_vec2(i) >0 && sign_vec3(i) <0
                         Ecal_cyl = Ecal_cyl*RAmat_cyl_ppm;
                          ProdMat2_cyl = ProdMat2_cyl*RAmatT_cyl_ppm; 
                         Ecal_plane = Ecal_plane*RAmat_plane_ppm;
                          ProdMat2_plane = ProdMat2_plane*RAmatT_plane_ppm; 
                     elseif sign_vec1(i) >0 && sign_vec2(i) <0 && sign_vec3(i) >0
                         Ecal_cyl = Ecal_cyl*RAmat_cyl_pmp;
                         ProdMat2_cyl = ProdMat2_cyl*RAmatT_cyl_pmp;
                         Ecal_plane = Ecal_plane*RAmat_plane_pmp;
                         ProdMat2_plane = ProdMat2_plane*RAmatT_plane_pmp;
                     elseif sign_vec1(i) >0 && sign_vec2(i) <0 && sign_vec3(i) <0
                         Ecal_cyl = Ecal_cyl*RAmat_cyl_pmm;
                         ProdMat2_cyl = ProdMat2_cyl*RAmatT_cyl_pmm;
                         Ecal_plane = Ecal_plane*RAmat_plane_pmm;
                         ProdMat2_plane = ProdMat2_plane*RAmatT_plane_pmm;
                     elseif sign_vec1(i) <0 && sign_vec2(i) >0 && sign_vec3(i) >0
                        Ecal_cyl = Ecal_cyl*RAmat_cyl_mpp;
                        ProdMat2_cyl = ProdMat2_cyl*RAmatT_cyl_mpp;
                        Ecal_plane = Ecal_plane*RAmat_plane_mpp;
                        ProdMat2_plane = ProdMat2_plane*RAmatT_plane_mpp;
                     elseif sign_vec1(i) <0 && sign_vec2(i) >0 && sign_vec3(i) <0
                         Ecal_cyl = Ecal_cyl*RAmat_cyl_mpm;
                         ProdMat2_cyl = ProdMat2_cyl*RAmatT_cyl_mpm;
                         Ecal_plane = Ecal_plane*RAmat_plane_mpm;
                         ProdMat2_plane = ProdMat2_plane*RAmatT_plane_mpm;
                     elseif sign_vec1(i) <0 && sign_vec2(i) <0 && sign_vec3(i) >0
                         Ecal_cyl = Ecal_cyl*RAmat_cyl_mmp;
                         ProdMat2_cyl = ProdMat2_cyl*RAmatT_cyl_mmp;
                         Ecal_plane = Ecal_plane*RAmat_plane_mmp;
                         ProdMat2_plane = ProdMat2_plane*RAmatT_plane_mmp;
                     else 
                         Ecal_cyl = Ecal_cyl*RAmat_cyl_mmm; 
                         ProdMat2_cyl = ProdMat2_cyl*RAmatT_cyl_mmm;
                         Ecal_plane = Ecal_plane*RAmat_plane_mmm; 
                         ProdMat2_plane = ProdMat2_plane*RAmatT_plane_mmm;
                     end                            
                end
                Ecal_cyl = Ecal_cyl*Rmat_cyl;            
                Ecal_plane = Ecal_plane*Rmat_plane;
                ProdMat2_cyl = ProdMat2_cyl*Rmat_cyl;
                ProdMat2_plane = ProdMat2_plane*Rmat_plane;
              Rmid_cyl=Rmat_cyl^dstmp(m);
              ePerp(m,indj)=Ecal_cyl*Rmid_cyl*ProdMat2_cyl*Smat_cyl;    
              Rmid_plane=Rmat_plane^dstmp(m);
              ePar(m,indj)=Ecal_plane*Rmid_plane*ProdMat2_plane*Smat_plane; 

            end

        end
    
    end
  end
E = ePar.*ePerp;
%disp(ePerp)
if strcmp(protocol.complex,'complex')
  E=[real(sum(E.*Rweight,2));imag(sum(E.*Rweight,2))];
elseif strcmp(protocol.complex,'real')
  E=real(sum(E.*Rweight,2));
elseif strcmp(protocol.complex,'abs')
  E=abs(sum(E.*Rweight,2));
end


% Compute the Jacobian matrix
if(nargout>1)
    J = zeros(length(E),length(x));
    dx = protocol.pert;
    if nargin < 3 
        
        for i = 1:length(x) 
        xpert = x;
            if i <=4 
                protocol.diff=i;
           else
                protocol.diff = 0; % fibre direction
           end
        xpert(i) = xpert(i)*(1+dx);    
        Epert = GammaFiniteCylinders_MM_SWOGSE_3D(xpert,protocol);
        dEtdx = (Epert - E)/(xpert(i)*dx);
        J(:,i) = dEtdx;
        end       
      
    else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
      
        for i = 1:length(x)
            if x_deriv(i) ~= 0  
               
                xpert = x;
                if i <=4
                protocol.diff = i;
                else 
                protocol.diff = 0;
                end
                
                xpert(i) = xpert(i)*(1+dx);    
                Epert = GammaFiniteCylinders_MM_SWOGSE_3D(xpert,protocol,x_deriv);
                dEtdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEtdx;
            end
        end
        
         
    end   
    protocol.diff=0;
   
end
    


    
    