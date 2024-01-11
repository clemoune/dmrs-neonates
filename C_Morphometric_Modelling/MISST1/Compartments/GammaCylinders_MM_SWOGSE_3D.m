function [E J]=GammaCylinders_MM_SWOGSE_3D(x,protocol,x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a GammaCylinders compartment.
% 
% [E,J]=GammaCylinders_MM_SWOGSE_3D(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate consisting of parallel
% cylinders with a Gamma distribution of radii and a diffusion 
% protocol specified in the input
% Substrate: Parallel, impermeable cylinders with a Gamma distribution of radii
% Diffusion pulse sequence: Square wave oscillating gradients with varying
% gradient orientation (SWOGSE_3D)
% Signal approximation: Matrix method (MM) - propagator expressed via 
%       eigenmode expansion
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% x - size 5 vector of model parameters in SI units for GammaCylinders:
%       x(1) - free diffusivity of the material inside the cylinders.
%       x(2) - mean radius of the distribution
%       x(3) - shape parameter of the Gamma distribution
%       x(4) - polar angle theta in spherical coordinates desbribing the fibre
% direction
%       x(5) - azimuthal angle phi in spherical coordinates describing the
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

GAMMA = 2.675987E8;
% model parameters
dRes = x(1);
theta = x(4);
phi = x(5);
fibredir = [sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];

cos_theta = cos(theta);
sin_theta = sqrt(1-cos_theta^2);
cos_phi = cos(phi);
sin_phi = sqrt(1-cos_phi^2);

 tau = protocol.tau;
 v = [cos_theta*cos_phi^2+sin_phi^2; -(1-cos_theta)*sin_phi*cos_phi; -sin_theta*cos_phi]; % vectors for the new x and y directions; see Ozarslan 2010
 w = [ -(1-cos_theta)*sin_phi*cos_phi; cos_theta*sin_phi^2+cos_phi^2; -sin_theta*sin_phi];
 K = floor((protocol.smalldel+1E-10)./tau)+1;
 dstmp =floor((protocol.delta-1E-10)./tau)-floor((protocol.smalldel+1E-10)./tau);


if protocol.angle == 1 % gradient is only along x
    Gx = protocol.Gx';
    M = size(Gx,1);
    ePerp =zeros(M,length(protocol.A_cyl));
    Bval=zeros(M,1);
    % calculate parallel signal
    G_dot_fibre = Gx.*fibredir(1);
   
  
      
    %disp('Angle method using exp(i(k-n)theta)')
    G_mag=Gx; 
    qhat_X=Gx./G_mag; % qhatx

    G_ind=round(G_mag./protocol.gstep);   
    if protocol.mirror == 1
        for m=1:M
             time_vec = tau*(0:K(m)-1);
            sign_vec1 = (-1).^floor((protocol.omegax(m)*time_vec-protocol.phix(m))./pi-1E-10);
            vec1 = zeros(size(sign_vec1));
            vec1(sign_vec1>=0) = G_dot_fibre(m); vec1(sign_vec1<0) = -G_dot_fibre(m); vec1(1)=0;vec1(end)=0;
            Gpar_vec = [vec1 zeros(1,dstmp(m)) -fliplr(vec1)];
       
            % Calculating hindered component Eh
            Fpar = cumsum(Gpar_vec)*tau;
            Bval(m)=sum((Fpar.^2)*tau);
        end
     elseif protocol.mirror == 0
        for m=1:M
             time_vec = tau*(0:K(m)-1);
            sign_vec1 = (-1).^floor((protocol.omegax(m)*time_vec-protocol.phix(m))./pi-1E-10);
            vec1 = zeros(size(sign_vec1));
            vec1(sign_vec1>=0) = G_dot_fibre(m); vec1(sign_vec1<0) = -G_dot_fibre(m); vec1(1)=0; vec1(end)=0;
            Gpar_vec = [vec1 zeros(1,dstmp(m)) -vec1];
            
            % Calculating hindered component Eh
            Fpar = cumsum(Gpar_vec)*tau;
            Bval(m)=sum((Fpar.^2)*tau);   
        end
    end
    
    for indj = 1:length(protocol.A_cyl) 
 
        % Perpendicular
        if protocol.diff==0
            Amat=protocol.A_cyl{indj};
            Smat=protocol.S_cyl{indj};
            Rmat=protocol.R_cyl{indj};
            Rweight = repmat(protocol.Rweight,[M,1]);

        elseif protocol.diff==1
            Amat=protocol.A_cyl{indj};
            Smat=protocol.S_cyl{indj};
            Rmat=protocol.RpertD_cyl{indj};
            Rweight = repmat(protocol.Rweight,[M,1]);

        elseif protocol.diff==2
            Amat=protocol.Aperta_cyl{indj};
            Smat=protocol.Sperta_cyl{indj};
            Rmat=protocol.Rperta_cyl{indj};
            Rweight = repmat(protocol.Rweight_perta,[M,1]);

        elseif protocol.diff==3
            Amat=protocol.Apertsh_cyl{indj};
            Smat=protocol.Spertsh_cyl{indj};
            Rmat=protocol.Rpertsh_cyl{indj};  
            Rweight = repmat(protocol.Rweight_pertsh,[M,1]);

        else
            error('protocol.diff is not suitable')
        end
        kDn=protocol.kDn_cyl{indj};
    
        if protocol.mirror == 1
            for m=1:M            

                Ecal=Smat';
                grad_dir_p = [qhat_X(m); 0; 0];
                Amat_p=(real(Amat)+(imag(Amat).*kDn)*dot(grad_dir_p,w))+1i*imag(Amat)*dot(grad_dir_p,v); 
                Amat_m=Amat_p'; 
                RAmat_p = Rmat*Amat_p^G_ind(m);
                RAmat_m = Rmat*Amat_m^G_ind(m);

                Ecal = Ecal*Rmat;  % first gradient point is 0
                for i = 2:K(m)-1
                     if sign_vec1(i) > 0  
                         Ecal = Ecal*RAmat_p;                  
                     else
                         Ecal = Ecal*RAmat_m;
                     end                        
                end
                 Ecal = Ecal*Rmat; % last gradient point is 0
                  Rmid=Rmat^dstmp(m);
                  ePerp(m,indj)=Ecal*Rmid*Ecal';                  
             end

        elseif protocol.mirror == 0     
             for m=1:M
                Ecal=Smat';
                grad_dir_p = [qhat_X(m); 0; 0];
                Amat_p=(real(Amat)+(imag(Amat).*kDn)*dot(grad_dir_p,w))+1i*imag(Amat)*dot(grad_dir_p,v); 
                Amat_m=Amat_p'; 
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
                 Ecal = Ecal*Rmat; % last gradient point is 0
                 ProdMat2 = ProdMat2*Rmat;
                  Rmid=Rmat^dstmp(m);
                  ePerp(m,indj)=Ecal*Rmid*ProdMat2*Smat;                 
             end     

        end
    end
    Bval=GAMMA^2*repmat(Bval,[1,length(protocol.A_cyl)]);
    ePar=exp(-Bval.*dRes);
    E_r = ePar.*ePerp;
    %Total Signal
    E=E_r;
   
elseif protocol.angle == 2 % Gx and Gy
     Gx = protocol.Gx';
     Gy = protocol.Gy';
    M = size(Gx,1);
    ePerp=zeros(M,length(protocol.A_cyl));
    Bval=zeros(M,1);
    % calculate parallel signal
    G_dot_fibre_pp = Gx.*fibredir(1) + Gy.*fibredir(2);
    G_dot_fibre_pm = Gx.*fibredir(1) - Gy.*fibredir(2);
    G_dot_fibre_mp = -Gx.*fibredir(1) + Gy.*fibredir(2);
    G_dot_fibre_mm = -G_dot_fibre_pp;
    

      %disp('Angle method using exp(i(k-n)theta)')
    G_mag=sqrt(Gx.^2+Gy.^2); 
    qhat_X=Gx./G_mag; % qhatx
    qhat_Y=Gy./G_mag; % qhaty

    G_ind=round(G_mag./protocol.gstep);   
    if protocol.mirror == 1
        for m=1:M
            time_vec = tau*(0:K(m)-1);
            sign_vec1 = (-1).^floor((protocol.omegax(m)*time_vec-protocol.phix(m))./pi-1E-10);
            sign_vec2 = (-1).^floor((protocol.omegay(m)*time_vec-protocol.phiy(m))./pi-1E-10);
            vec1 = zeros(size(sign_vec1)); 
            vec1(sign_vec1>0 & sign_vec2>0) = G_dot_fibre_pp(m);
            vec1(sign_vec1>0 & sign_vec2<0) = G_dot_fibre_pm(m);
            vec1(sign_vec1<0 & sign_vec2>0) = G_dot_fibre_mp(m);
            vec1(sign_vec1<0 & sign_vec2<0) = G_dot_fibre_mm(m);
            vec1(1)=0; vec1(end)=0;
            Gpar_vec = [vec1 zeros(1,dstmp(m)) -fliplr(vec1)];
            
             % Calculating parallel component 
            Fpar = cumsum(Gpar_vec)*tau;
            Bval(m)=sum((Fpar.^2)*tau);
        end
    elseif protocol.mirror == 0
        for m=1:M
            time_vec = tau*(0:K(m)-1);
            sign_vec1 = (-1).^floor((protocol.omegax(m)*time_vec-protocol.phix(m))./pi-1E-10);
            sign_vec2 = (-1).^floor((protocol.omegay(m)*time_vec-protocol.phiy(m))./pi-1E-10);
            vec1 = zeros(size(sign_vec1)); 
            vec1(sign_vec1>0 & sign_vec2>0) = G_dot_fibre_pp(m);
            vec1(sign_vec1>0 & sign_vec2<0) = G_dot_fibre_pm(m);
            vec1(sign_vec1<0 & sign_vec2>0) = G_dot_fibre_mp(m);
            vec1(sign_vec1<0 & sign_vec2<0) = G_dot_fibre_mm(m);
            vec1(1)=0; vec1(end)=0;
            Gpar_vec = [vec1 zeros(1,dstmp(m)) -vec1];
            % Calculating parallel component 
            Fpar = cumsum(Gpar_vec)*tau;
            Bval(m)=sum((Fpar.^2)*tau);
        end
    end
    
    for indj = 1:length(protocol.A_cyl) 
 
        % Perpendicular
        if protocol.diff==0
            Amat=protocol.A_cyl{indj};
            Smat=protocol.S_cyl{indj};
            Rmat=protocol.R_cyl{indj};
            Rweight = repmat(protocol.Rweight,[M,1]);

        elseif protocol.diff==1
            Amat=protocol.A_cyl{indj};
            Smat=protocol.S_cyl{indj};
            Rmat=protocol.RpertD_cyl{indj};
            Rweight = repmat(protocol.Rweight,[M,1]);

        elseif protocol.diff==2
            Amat=protocol.Aperta_cyl{indj};
            Smat=protocol.Sperta_cyl{indj};
            Rmat=protocol.Rperta_cyl{indj};
            Rweight = repmat(protocol.Rweight_perta,[M,1]);

        elseif protocol.diff==3
            Amat=protocol.Apertsh_cyl{indj};
            Smat=protocol.Spertsh_cyl{indj};
            Rmat=protocol.Rpertsh_cyl{indj};  
            Rweight = repmat(protocol.Rweight_pertsh,[M,1]);

        else
            error('protocol.diff is not suitable')
        end
        kDn=protocol.kDn_cyl{indj};
 
           
         if protocol.mirror == 1
            for m=1:M       
                Ecal=Smat';
                grad_dir_pp = [qhat_X(m); qhat_Y(m); 0];
                grad_dir_pm = [qhat_X(m); -qhat_Y(m); 0];
                grad_dir_mp = [-qhat_X(m); qhat_Y(m); 0];

                Amat_pp=(real(Amat)+(imag(Amat).*kDn)*dot(grad_dir_pp,w))+1i*imag(Amat)*dot(grad_dir_pp,v);
                Amat_pm=(real(Amat)+(imag(Amat).*kDn)*dot(grad_dir_pm,w))+1i*imag(Amat)*dot(grad_dir_pm,v); 
                Amat_mp=(real(Amat)+(imag(Amat).*kDn)*dot(grad_dir_mp,w))+1i*imag(Amat)*dot(grad_dir_mp,v); 
                Amat_mm=Amat_pp'; 
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
                 Ecal = Ecal*Rmat; % last gradient point is 0

                  Rmid=Rmat^dstmp(m);
                  ePerp(m,indj)=Ecal*Rmid*Ecal';                  
             end
        elseif protocol.mirror == 0        
            for m=1:M

                 Ecal=Smat';
                grad_dir_pp = [qhat_X(m); qhat_Y(m); 0];
                grad_dir_pm = [qhat_X(m); -qhat_Y(m); 0];
                grad_dir_mp = [-qhat_X(m); qhat_Y(m); 0];

                Amat_pp=(real(Amat)+(imag(Amat).*kDn)*dot(grad_dir_pp,w))+1i*imag(Amat)*dot(grad_dir_pp,v);
                Amat_pm=(real(Amat)+(imag(Amat).*kDn)*dot(grad_dir_pm,w))+1i*imag(Amat)*dot(grad_dir_pm,v); 
                Amat_mp=(real(Amat)+(imag(Amat).*kDn)*dot(grad_dir_mp,w))+1i*imag(Amat)*dot(grad_dir_mp,v); 
                Amat_mm=Amat_pp'; 
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

                 Ecal = Ecal*Rmat; % last gradient point is 0
                 ProdMat2 = ProdMat2*Rmat;
                  Rmid=Rmat^dstmp(m);
                  ePerp(m,indj)=Ecal*Rmid*ProdMat2*Smat;                 
             end     

         end
    end
    Bval=GAMMA^2*repmat(Bval,[1,length(protocol.A_cyl)]);
    ePar=exp(-Bval.*dRes);
    E_r = ePar.*ePerp;
    %Total Signal
    E=E_r;
    
elseif protocol.angle == 3 % Gx and Gz
     Gx = protocol.Gx';
     Gz = protocol.Gz';
    M = size(Gx,1);
    ePerp=zeros(M,length(protocol.A_cyl));
    Bval=zeros(M,1);
    % calculate parallel signal
    G_dot_fibre_pp = Gx.*fibredir(1) + Gz.*fibredir(3);
    G_dot_fibre_pm = Gx.*fibredir(1) - Gz.*fibredir(3);
    G_dot_fibre_mp = -Gx.*fibredir(1) + Gz.*fibredir(3);
    G_dot_fibre_mm = -G_dot_fibre_pp;
    

      %disp('Angle method using exp(i(k-n)theta)')
    G_mag=sqrt(Gx.^2+Gz.^2); 
    qhat_X=Gx./G_mag; % qhatx
    qhat_Z=Gz./G_mag; % qhaty

    G_ind=round(G_mag./protocol.gstep);   
    if protocol.mirror == 1       
        for m=1:M
            time_vec = tau*(0:K(m)-1);
            sign_vec1 = (-1).^floor((protocol.omegax(m)*time_vec-protocol.phix(m))./pi-1E-10);
            sign_vec2 = (-1).^floor((protocol.omegaz(m)*time_vec-protocol.phiz(m))./pi-1E-10);
            vec1 = zeros(size(sign_vec1)); 
            vec1(sign_vec1>0 & sign_vec2>0) = G_dot_fibre_pp(m);
            vec1(sign_vec1>0 & sign_vec2<0) = G_dot_fibre_pm(m);
            vec1(sign_vec1<0 & sign_vec2>0) = G_dot_fibre_mp(m);
            vec1(sign_vec1<0 & sign_vec2<0) = G_dot_fibre_mm(m);
            vec1(1)=0; vec1(end)=0;
            Gpar_vec = [vec1 zeros(1,dstmp(m)) -fliplr(vec1)];
 
            % Calculating parallel component 
            Fpar = cumsum(Gpar_vec)*tau;
            Bval(m)=sum((Fpar.^2)*tau);
        end
      elseif protocol.mirror == 0
        for m=1:M
            time_vec = tau*(0:K(m)-1);
            sign_vec1 = (-1).^floor((protocol.omegax(m)*time_vec-protocol.phix(m))./pi-1E-10);
            sign_vec2 = (-1).^floor((protocol.omegaz(m)*time_vec-protocol.phiz(m))./pi-1E-10);
            vec1 = zeros(size(sign_vec1)); 
            vec1(sign_vec1>0 & sign_vec2>0) = G_dot_fibre_pp(m);
            vec1(sign_vec1>0 & sign_vec2<0) = G_dot_fibre_pm(m);
            vec1(sign_vec1<0 & sign_vec2>0) = G_dot_fibre_mp(m);
            vec1(sign_vec1<0 & sign_vec2<0) = G_dot_fibre_mm(m);
            vec1(1)=0; vec1(end)=0;
            Gpar_vec = [vec1 zeros(1,dstmp(m)) -vec1];
            % Calculating parallel component 
            Fpar = cumsum(Gpar_vec)*tau;
            Bval(m)=sum((Fpar.^2)*tau);  
        end
    end
     for indj = 1:length(protocol.A_cyl) 
 
        % Perpendicular
        if protocol.diff==0
            Amat=protocol.A_cyl{indj};
            Smat=protocol.S_cyl{indj};
            Rmat=protocol.R_cyl{indj};
            Rweight = repmat(protocol.Rweight,[M,1]);

        elseif protocol.diff==1
            Amat=protocol.A_cyl{indj};
            Smat=protocol.S_cyl{indj};
            Rmat=protocol.RpertD_cyl{indj};
            Rweight = repmat(protocol.Rweight,[M,1]);

        elseif protocol.diff==2
            Amat=protocol.Aperta_cyl{indj};
            Smat=protocol.Sperta_cyl{indj};
            Rmat=protocol.Rperta_cyl{indj};
            Rweight = repmat(protocol.Rweight_perta,[M,1]);

        elseif protocol.diff==3
            Amat=protocol.Apertsh_cyl{indj};
            Smat=protocol.Spertsh_cyl{indj};
            Rmat=protocol.Rpertsh_cyl{indj};  
            Rweight = repmat(protocol.Rweight_pertsh,[M,1]);

        else
            error('protocol.diff is not suitable')
        end
        kDn=protocol.kDn_cyl{indj};
   
         if protocol.mirror == 1       
            for m=1:M

                Ecal=Smat';
                grad_dir_pp = [qhat_X(m); 0; qhat_Z(m); ];
                grad_dir_pm = [qhat_X(m); 0; -qhat_Z(m); ];
                grad_dir_mp = [-qhat_X(m); 0; qhat_Z(m); ];

                Amat_pp=(real(Amat)+(imag(Amat).*kDn)*dot(grad_dir_pp,w))+1i*imag(Amat)*dot(grad_dir_pp,v);
                Amat_pm=(real(Amat)+(imag(Amat).*kDn)*dot(grad_dir_pm,w))+1i*imag(Amat)*dot(grad_dir_pm,v); 
                Amat_mp=(real(Amat)+(imag(Amat).*kDn)*dot(grad_dir_mp,w))+1i*imag(Amat)*dot(grad_dir_mp,v); 
                Amat_mm=Amat_pp'; 
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
                 Ecal = Ecal*Rmat; % last gradient point is 0

                  Rmid=Rmat^dstmp(m);
                  ePerp(m,indj)=Ecal*Rmid*Ecal';                  
            end
         elseif protocol.mirror == 0
            for m=1:M    

                 Ecal=Smat';
                grad_dir_pp = [qhat_X(m); 0; qhat_Z(m); ];
                grad_dir_pm = [qhat_X(m); 0; -qhat_Z(m); ];
                grad_dir_mp = [-qhat_X(m); 0; qhat_Z(m); ];

                Amat_pp=(real(Amat)+(imag(Amat).*kDn)*dot(grad_dir_pp,w))+1i*imag(Amat)*dot(grad_dir_pp,v);
                Amat_pm=(real(Amat)+(imag(Amat).*kDn)*dot(grad_dir_pm,w))+1i*imag(Amat)*dot(grad_dir_pm,v); 
                Amat_mp=(real(Amat)+(imag(Amat).*kDn)*dot(grad_dir_mp,w))+1i*imag(Amat)*dot(grad_dir_mp,v); 
                Amat_mm=Amat_pp'; 
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
                   Ecal = Ecal*Rmat; % last gradient point is 0
                 ProdMat2 = ProdMat2*Rmat; 

                  Rmid=Rmat^dstmp(m);
                  ePerp(m,indj)=Ecal*Rmid*ProdMat2*Smat;                 
             end     

         end
    end
    Bval=GAMMA^2*repmat(Bval,[1,length(protocol.A_cyl)]);
    ePar=exp(-Bval.*dRes);
    E_r = ePar.*ePerp;
    %Total Signal
    E=E_r;
    if strcmp(protocol.complex,'complex')
      E=[real(E);imag(E)];
    elseif strcmp(protocol.complex,'real')
      E=real(E);
    elseif strcmp(protocol.complex,'abs')
      E=abs(E);
    end  
elseif protocol.angle == 4 % Gx, Gy and Gz
     Gx = protocol.Gx';
     Gy = protocol.Gy';
     Gz = protocol.Gz';
    M = size(Gx,1);
    ePerp=zeros(M,length(protocol.A_cyl));
    Bval=zeros(M,1);
    % calculate parallel signal
    G_dot_fibre_ppp = Gx.*fibredir(1) +Gy.*fibredir(2) + Gz.*fibredir(3);
    G_dot_fibre_ppm = Gx.*fibredir(1) +Gy.*fibredir(2) - Gz.*fibredir(3);
    G_dot_fibre_pmp = Gx.*fibredir(1) -Gy.*fibredir(2) + Gz.*fibredir(3);
    G_dot_fibre_pmm = Gx.*fibredir(1) -Gy.*fibredir(2) - Gz.*fibredir(3);
    G_dot_fibre_mpp = -Gx.*fibredir(1) +Gy.*fibredir(2) + Gz.*fibredir(3);
    G_dot_fibre_mpm = -Gx.*fibredir(1) +Gy.*fibredir(2) - Gz.*fibredir(3);
    G_dot_fibre_mmp = -Gx.*fibredir(1) -Gy.*fibredir(2) + Gz.*fibredir(3);
    G_dot_fibre_mmm = -Gx.*fibredir(1) -Gy.*fibredir(2) - Gz.*fibredir(3);
    

      %disp('Angle method using exp(i(k-n)theta)')
    G_mag=sqrt(Gx.^2+Gy.^2+Gz.^2); 
    qhat_X=Gx./G_mag; % qhatx
    qhat_Y=Gy./G_mag; % qhatx
    qhat_Z=Gz./G_mag; % qhaty


    G_ind=round(G_mag./protocol.gstep);   
    if protocol.mirror == 1
        for m=1:M
            time_vec = tau*(0:K(m)-1);
            sign_vec1 = (-1).^floor((protocol.omegax(m)*time_vec-protocol.phix(m))./pi-1E-10);
            sign_vec2 = (-1).^floor((protocol.omegay(m)*time_vec-protocol.phiy(m))./pi-1E-10);
            sign_vec3 = (-1).^floor((protocol.omegaz(m)*time_vec-protocol.phiz(m))./pi-1E-10);
            vec1 = zeros(size(sign_vec1)); 
            vec1(sign_vec1>0 & sign_vec2>0 & sign_vec3>0) = G_dot_fibre_ppp(m);
            vec1(sign_vec1>0 & sign_vec2>0 & sign_vec3<0) = G_dot_fibre_ppm(m);
            vec1(sign_vec1>0 & sign_vec2<0 & sign_vec3>0) = G_dot_fibre_pmp(m);
            vec1(sign_vec1>0 & sign_vec2<0 & sign_vec3<0) = G_dot_fibre_pmm(m);
            vec1(sign_vec1<0 & sign_vec2>0 & sign_vec3>0) = G_dot_fibre_mpp(m);
            vec1(sign_vec1<0 & sign_vec2>0 & sign_vec3<0) = G_dot_fibre_mpm(m);
            vec1(sign_vec1<0 & sign_vec2<0 & sign_vec3>0) = G_dot_fibre_mmp(m);
            vec1(sign_vec1<0 & sign_vec2<0 & sign_vec3<0) = G_dot_fibre_mmm(m);
            vec1(1) = 0; vec1(end)=0;           
            Gpar_vec = [vec1 zeros(1,dstmp(m)) -fliplr(vec1)];
           
            % Calculating parallel component 
            Fpar = cumsum(Gpar_vec*tau);
            Bval(m)=sum((Fpar.^2)*tau);
        end
      elseif protocol.mirror == 0
        for m=1:M        
            time_vec = tau*(0:K(m)-1);
            sign_vec1 = (-1).^floor((protocol.omegax(m)*time_vec-protocol.phix(m))./pi-1E-10);
            sign_vec2 = (-1).^floor((protocol.omegay(m)*time_vec-protocol.phiy(m))./pi-1E-10);
            sign_vec3 = (-1).^floor((protocol.omegaz(m)*time_vec-protocol.phiz(m))./pi-1E-10);
            vec1 = zeros(size(sign_vec1)); 
            vec1(sign_vec1>0 & sign_vec2>0 & sign_vec3>0) = G_dot_fibre_ppp(m);
            vec1(sign_vec1>0 & sign_vec2>0 & sign_vec3<0) = G_dot_fibre_ppm(m);
            vec1(sign_vec1>0 & sign_vec2<0 & sign_vec3>0) = G_dot_fibre_pmp(m);
            vec1(sign_vec1>0 & sign_vec2<0 & sign_vec3<0) = G_dot_fibre_pmm(m);
            vec1(sign_vec1<0 & sign_vec2>0 & sign_vec3>0) = G_dot_fibre_mpp(m);
            vec1(sign_vec1<0 & sign_vec2>0 & sign_vec3<0) = G_dot_fibre_mpm(m);
            vec1(sign_vec1<0 & sign_vec2<0 & sign_vec3>0) = G_dot_fibre_mmp(m);
            vec1(sign_vec1<0 & sign_vec2<0 & sign_vec3<0) = G_dot_fibre_mmm(m);
            vec1(1) = 0; vec1(end)=0;            
            Gpar_vec = [vec1 zeros(1,dstmp(m)) -(vec1)];
           
            % Calculating parallel component 
            Fpar = cumsum(Gpar_vec*tau);
            Bval(m)=sum((Fpar.^2)*tau);
        end
     end
     for indj = 1:length(protocol.A_cyl) 
 
        % Perpendicular
        if protocol.diff==0
            Amat=protocol.A_cyl{indj};
            Smat=protocol.S_cyl{indj};
            Rmat=protocol.R_cyl{indj};
            Rweight = repmat(protocol.Rweight,[M,1]);

        elseif protocol.diff==1
            Amat=protocol.A_cyl{indj};
            Smat=protocol.S_cyl{indj};
            Rmat=protocol.RpertD_cyl{indj};
            Rweight = repmat(protocol.Rweight,[M,1]);

        elseif protocol.diff==2
            Amat=protocol.Aperta_cyl{indj};
            Smat=protocol.Sperta_cyl{indj};
            Rmat=protocol.Rperta_cyl{indj};
            Rweight = repmat(protocol.Rweight_perta,[M,1]);

        elseif protocol.diff==3
            Amat=protocol.Apertsh_cyl{indj};
            Smat=protocol.Spertsh_cyl{indj};
            Rmat=protocol.Rpertsh_cyl{indj};  
            Rweight = repmat(protocol.Rweight_pertsh,[M,1]);

        else
            error('protocol.diff is not suitable')
        end
        kDn=protocol.kDn_cyl{indj};       
          if protocol.mirror == 1
            for m=1:M      
                Ecal=Smat';
                grad_dir_ppp = [qhat_X(m); qhat_Y(m); qhat_Z(m); ];
                grad_dir_ppm = [qhat_X(m); qhat_Y(m); -qhat_Z(m); ];
                grad_dir_pmp = [qhat_X(m); -qhat_Y(m); qhat_Z(m); ];
                grad_dir_pmm = [qhat_X(m); -qhat_Y(m); -qhat_Z(m); ];
                grad_dir_mpp = [-qhat_X(m); qhat_Y(m); qhat_Z(m); ];
                grad_dir_mpm = [-qhat_X(m); qhat_Y(m); -qhat_Z(m); ];
                grad_dir_mmp = [-qhat_X(m); -qhat_Y(m); qhat_Z(m); ];
                grad_dir_mmm = [-qhat_X(m); -qhat_Y(m); -qhat_Z(m); ];

                Amat_ppp=(real(Amat)+(imag(Amat).*kDn)*dot(grad_dir_ppp,w))+1i*imag(Amat)*dot(grad_dir_ppp,v);
                Amat_ppm=(real(Amat)+(imag(Amat).*kDn)*dot(grad_dir_ppm,w))+1i*imag(Amat)*dot(grad_dir_ppm,v);
                Amat_pmp=(real(Amat)+(imag(Amat).*kDn)*dot(grad_dir_pmp,w))+1i*imag(Amat)*dot(grad_dir_pmp,v);
                Amat_pmm=(real(Amat)+(imag(Amat).*kDn)*dot(grad_dir_pmm,w))+1i*imag(Amat)*dot(grad_dir_pmm,v);
                Amat_mpp=(real(Amat)+(imag(Amat).*kDn)*dot(grad_dir_mpp,w))+1i*imag(Amat)*dot(grad_dir_mpp,v);
                Amat_mpm=(real(Amat)+(imag(Amat).*kDn)*dot(grad_dir_mpm,w))+1i*imag(Amat)*dot(grad_dir_mpm,v);
                Amat_mmp=(real(Amat)+(imag(Amat).*kDn)*dot(grad_dir_mmp,w))+1i*imag(Amat)*dot(grad_dir_mmp,v);
                Amat_mmm=(real(Amat)+(imag(Amat).*kDn)*dot(grad_dir_mmm,w))+1i*imag(Amat)*dot(grad_dir_mmm,v); 

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
                 Ecal = Ecal*Rmat; % last gradient point is 0

                  Rmid=Rmat^dstmp(m);
                  ePerp(m,indj)=Ecal*Rmid*Ecal';                  
             end
        elseif protocol.mirror == 0
            for m=1:M     

                 Ecal=Smat';
                grad_dir_ppp = [qhat_X(m); qhat_Y(m); qhat_Z(m); ];
                grad_dir_ppm = [qhat_X(m); qhat_Y(m); -qhat_Z(m); ];
                grad_dir_pmp = [qhat_X(m); -qhat_Y(m); qhat_Z(m); ];
                grad_dir_pmm = [qhat_X(m); -qhat_Y(m); -qhat_Z(m); ];
                grad_dir_mpp = [-qhat_X(m); qhat_Y(m); qhat_Z(m); ];
                grad_dir_mpm = [-qhat_X(m); qhat_Y(m); -qhat_Z(m); ];
                grad_dir_mmp = [-qhat_X(m); -qhat_Y(m); qhat_Z(m); ];
                grad_dir_mmm = [-qhat_X(m); -qhat_Y(m); -qhat_Z(m); ];

                Amat_ppp=(real(Amat)+(imag(Amat).*kDn)*dot(grad_dir_ppp,w))+1i*imag(Amat)*dot(grad_dir_ppp,v);
                Amat_ppm=(real(Amat)+(imag(Amat).*kDn)*dot(grad_dir_ppm,w))+1i*imag(Amat)*dot(grad_dir_ppm,v);
                Amat_pmp=(real(Amat)+(imag(Amat).*kDn)*dot(grad_dir_pmp,w))+1i*imag(Amat)*dot(grad_dir_pmp,v);
                Amat_pmm=(real(Amat)+(imag(Amat).*kDn)*dot(grad_dir_pmm,w))+1i*imag(Amat)*dot(grad_dir_pmm,v);
                Amat_mpp=(real(Amat)+(imag(Amat).*kDn)*dot(grad_dir_mpp,w))+1i*imag(Amat)*dot(grad_dir_mpp,v);
                Amat_mpm=(real(Amat)+(imag(Amat).*kDn)*dot(grad_dir_mpm,w))+1i*imag(Amat)*dot(grad_dir_mpm,v);
                Amat_mmp=(real(Amat)+(imag(Amat).*kDn)*dot(grad_dir_mmp,w))+1i*imag(Amat)*dot(grad_dir_mmp,v);
                Amat_mmm=(real(Amat)+(imag(Amat).*kDn)*dot(grad_dir_mmm,w))+1i*imag(Amat)*dot(grad_dir_mmm,v); 

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
                 Ecal = Ecal*Rmat; % last gradient point is 0
                 ProdMat2 = ProdMat2*Rmat;
                  Rmid=Rmat^dstmp(m);
                  ePerp(m,indj)=Ecal*Rmid*ProdMat2*Smat;                 
             end     

          end 
     end
  
    Bval=GAMMA^2*repmat(Bval,[1,length(protocol.A_cyl)]);
    ePar=exp(-Bval.*dRes);
    E_r = ePar.*ePerp;
    %Total Signal
    E=E_r;  
else error('Unknown protocol.angle')
end

if strcmp(protocol.complex,'complex')
  E=[real(sum(E.*Rweight,2));imag(sum(E.*Rweight,2))];
elseif strcmp(protocol.complex,'real')
  E=real(sum(E.*Rweight,2));
elseif strcmp(protocol.complex,'abs')
  E=abs(sum(E.*Rweight,2));
end   


% Compute the Jacobian matrix
if(nargout>1)
    J = zeros(length(E),length(x)); % includes fibre direction
    dx = protocol.pert;
    if nargin < 3 
        
        for i = 1:length(x)
        xpert = x;
            if i<=3
            protocol.diff=i;
            else
            protocol.diff=0; % fibre diection
            end
        xpert(i) = xpert(i)*(1+dx);    
        Epert = GammaCylinders_MM_SWOGSE_3D(xpert,protocol);
        dEtdx = (Epert - E)/(xpert(i)*dx);
         J(:,i) = dEtdx;
        end      

        
    else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
     
        for i = 1:length(x)    
            if x_deriv(i) ~= 0  
               
                xpert = x;
                    if i<=3
                    protocol.diff=i;
                    else
                    protocol.diff = 0; % derivatives with respect to fibre direction
                    end
                xpert(i) = xpert(i)*(1+dx);    
                Epert = GammaCylinders_MM_SWOGSE_3D(xpert,protocol,x_deriv);
                dEtdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEtdx;
            end
        end
         
    end  
    protocol.diff=0;
  
   
end