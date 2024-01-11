function [E J]=GammaFiniteCylinders_MM_GEN(x,protocol,x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a GammaFiniteCylinder compartment.
% 
% [E,J]=GammaFiniteCylinder_MM_GEN(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate consisting of parallel
% finite cylinders with a Gamma distribution of radii and a diffusion 
% protocol specified in the input
% Substrate: Parallel, impermeable finite cylinders with a Gamma
% distribution of radii
% Diffusion pulse sequence: Generalized waveform (GEN)
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


if ~isfield(protocol,'complex')
  protocol.complex='real';
end
% Model parameters
theta = x(5);
phi = x(6);
fibredir = [sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];
tau=protocol.tau;
G=protocol.G;
[M totsizeG]=size(G);

% Calculating the parallel & perpendicular signal is similar  (parallel plane & cylindrical restriction )


ePerp=zeros(M,length(protocol.A_cyl)); % cylinder restriction
ePar = zeros(M,length(protocol.A_cyl)); % parallel planes 
cos_theta = cos(theta);
sin_theta = sqrt(1-cos_theta^2);
cos_phi = cos(phi);
sin_phi = sqrt(1-cos_phi^2);


 v = [cos_theta*cos_phi^2+sin_phi^2; -(1-cos_theta)*sin_phi*cos_phi; -sin_theta*cos_phi]; % vectors for the new x and y directions; see Ozarslan 2010
 w = [ -(1-cos_theta)*sin_phi*cos_phi; cos_theta*sin_phi^2+cos_phi^2; -sin_theta*sin_phi];



if isfield(protocol,'smalldel') && isfield(protocol,'delta')
     K = floor((protocol.smalldel+1E-10)./tau)+1;
     dstmp = floor((protocol.delta-1E-10)./protocol.tau)-floor((protocol.smalldel+1E-10)./protocol.tau);
    %disp('Angle method using exp(i(k-n)theta)')
      G_mag=sqrt(G(:,1:3:end).^2+G(:,2:3:end).^2+ G(:,3:3:end).^2);   
    indG=G_mag>1E-6;
    G_tmpX=G(:,1:3:end);
    G_tmpY=G(:,2:3:end);
    G_tmpZ=G(:,3:3:end);
    qhat_X = zeros(size(G_mag));
    qhat_Y = zeros(size(G_mag));
    qhat_Z = zeros(size(G_mag));
    qhat_X(indG) =  G_tmpX(indG)./G_mag(indG); % qhatx
    qhat_Y(indG) =  G_tmpY(indG)./G_mag(indG); % qhaty
    qhat_Z(indG) =  G_tmpZ(indG)./G_mag(indG); % qhatz

    Gind_mat=round(G_mag./protocol.gstep);  
    
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
        kDn_cyl=protocol.kDn_cyl{indj};
  
  
        if protocol.mirror == 1
            for m=1:M
               Ecal_cyl=Smat_cyl';
               Ecal_plane = Smat_plane';

                  for ind=1:K(m)
                    Gind=Gind_mat(m,ind);

                    if Gind~=0 % the angles are in the perpendicular plane
                      grad_dir = [qhat_X(m,ind); qhat_Y(m,ind); qhat_Z(m,ind)];
                      AmatU_cyl=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dir,w))+1i*imag(Amat_cyl)*dot(grad_dir,v); 
                      AmatUD_cyl=AmatU_cyl^Gind;
                      AmatU_plane = real(Amat_plane)+ 1i* imag(Amat_plane)*dot(grad_dir,fibredir);
                      AmatUD_plane=AmatU_plane^Gind;
                      Ecal_cyl=Ecal_cyl*Rmat_cyl*AmatUD_cyl;
                      Ecal_plane=Ecal_plane*Rmat_plane*AmatUD_plane;
                    else 
                      Ecal_cyl=Ecal_cyl*Rmat_cyl;
                      Ecal_plane=Ecal_plane*Rmat_plane;
                    end

                  end
              Rmid_cyl=Rmat_cyl^dstmp(m);
              ePerp(m,indj)=Ecal_cyl*Rmid_cyl*Ecal_cyl';   
              Rmid_plane=Rmat_plane^dstmp(m);
              ePar(m,indj)=Ecal_plane*Rmid_plane*Ecal_plane'; 

            end
        elseif protocol.mirror == 0
            for m=1:M
               Ecal_cyl=Smat_cyl';
               Ecal_plane = Smat_plane';      
               ProdMat2_cyl = eye(size(Amat_cyl)); 
               ProdMat2_plane = ProdMat2_cyl;
               for ind=1:K(m)
                    Gind=Gind_mat(m,ind);

                    if Gind~=0
                      grad_dir = [qhat_X(m,ind); qhat_Y(m,ind); qhat_Z(m,ind)];  
                      AmatU_cyl=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dir,w))+1i*imag(Amat_cyl)*dot(grad_dir,v); 
                      AmatUD_cyl=AmatU_cyl^Gind;
                      AmatUDT_cyl = AmatUD_cyl'; 
                      AmatU_plane = real(Amat_plane)+ 1i* imag(Amat_plane)*dot(grad_dir,fibredir);
                      AmatUD_plane=AmatU_plane^Gind;
                      AmatUDT_plane = AmatUD_plane';

                      Ecal_cyl=Ecal_cyl*Rmat_cyl*AmatUD_cyl;                     
                      ProdMat2_cyl = ProdMat2_cyl*Rmat_cyl*AmatUDT_cyl;  

                      Ecal_plane=Ecal_plane*Rmat_plane*AmatUD_plane;
                      ProdMat2_plane = ProdMat2_plane*Rmat_plane*AmatUDT_plane; 
                    else 
                     Ecal_cyl=Ecal_cyl*Rmat_cyl;                     
                     ProdMat2_cyl = ProdMat2_cyl*Rmat_cyl;    
                     Ecal_plane=Ecal_plane*Rmat_plane;
                     ProdMat2_plane = ProdMat2_plane*Rmat_plane; 
                    end                        

              end
              Rmid_cyl=Rmat_cyl^dstmp(m);
              Rmid_plane=Rmat_plane^dstmp(m);
              ePerp(m,indj)=Ecal_cyl*Rmid_cyl*ProdMat2_cyl*Smat_cyl;
              ePar(m,indj)=Ecal_plane*Rmid_plane*ProdMat2_plane*Smat_plane;  

            end
        end
    end
 
else
    K=totsizeG/3;
   %disp('Angle method using exp(i(k-n)theta)')
    G_mag=sqrt(G(:,1:3:end).^2+G(:,2:3:end).^2+ G(:,3:3:end).^2);   
    indG=G_mag>1E-6;
    G_tmpX=G(:,1:3:end);
    G_tmpY=G(:,2:3:end);
    G_tmpZ=G(:,3:3:end);
    qhat_X = zeros(size(G_mag));
    qhat_Y = zeros(size(G_mag));
    qhat_Z = zeros(size(G_mag));
    qhat_X(indG) =  G_tmpX(indG)./G_mag(indG); % qhatx
    qhat_Y(indG) =  G_tmpY(indG)./G_mag(indG); % qhaty
    qhat_Z(indG) =  G_tmpZ(indG)./G_mag(indG); % qhatz

    Gind_mat=round(G_mag./protocol.gstep);  

    
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
        kDn_cyl=protocol.kDn_cyl{indj};
       for m=1:M
            
           Ecal_cyl=Smat_cyl';
           Ecal_plane = Smat_plane';
           
             for ind=1:K
                Gind=Gind_mat(m,ind);
                
                if Gind~=0 % the angles are in the perpendicular plane
                  grad_dir = [qhat_X(m,ind); qhat_Y(m,ind); qhat_Z(m,ind)];
                  AmatU_cyl=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dir,w))+1i*imag(Amat_cyl)*dot(grad_dir,v); 
                  AmatUD_cyl=AmatU_cyl^Gind;
                  AmatU_plane = real(Amat_plane)+ 1i* imag(Amat_plane)*dot(grad_dir,fibredir);
                  AmatUD_plane=AmatU_plane^Gind;
                  Ecal_cyl=Ecal_cyl*Rmat_cyl*AmatUD_cyl;
                  Ecal_plane=Ecal_plane*Rmat_plane*AmatUD_plane;
                else 
                  Ecal_cyl=Ecal_cyl*Rmat_cyl;
                  Ecal_plane=Ecal_plane*Rmat_plane;
                end
               
              end
              ePerp(m,indj)=Ecal_cyl*Rmat_cyl*Smat_cyl;
              ePar(m,indj)=Ecal_plane*Rmat_plane*Smat_plane;
           
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
        Epert = GammaFiniteCylinders_MM_GEN(xpert,protocol);
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
                Epert = GammaFiniteCylinders_MM_GEN(xpert,protocol,x_deriv);
                dEtdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEtdx;
            end
        end
        
         
    end   
    protocol.diff=0;
   
end
    


    
    