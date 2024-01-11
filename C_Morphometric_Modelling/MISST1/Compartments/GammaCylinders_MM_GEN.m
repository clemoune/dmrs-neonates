function [E J]=GammaCylinders_MM_GEN(x,protocol,x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a GammaCylinder compartment.
% 
% [E,J]=GammaCylinder_MM_GEN(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate consisting of parallel
% cylinders with a Gamma distribution of radii and a diffusion 
% protocol specified in the input
% Substrate: Parallel, impermeable cylinders with a Gamma distribution of radii
% Diffusion pulse sequence: Generalized waveform (GEN)
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

dRes = x(1);
theta = x(4);
phi = x(5);



GAMMA = 2.675987E8;
% model parameters

fibredir = [sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];
tau=protocol.tau;
G=protocol.G;
[M totsizeG]=size(G);

ePerp=zeros(M,length(protocol.A_cyl) );
cos_theta = cos(theta);
sin_theta = sqrt(1-cos_theta^2);
cos_phi = cos(phi);
sin_phi = sqrt(1-cos_phi^2);


 v = [cos_theta*cos_phi^2+sin_phi^2; -(1-cos_theta)*sin_phi*cos_phi; -sin_theta*cos_phi]; % vectors for the new x and y directions; see Ozarslan 2010
 w = [ -(1-cos_theta)*sin_phi*cos_phi; cos_theta*sin_phi^2+cos_phi^2; -sin_theta*sin_phi];


K=totsizeG/3;
% Perpendicular
Bval=zeros(M,1);
for m=1:M
G_dot_fibre = G(m,1:3:end)*fibredir(1)+ G(m,2:3:end)*fibredir(2) + G(m,3:3:end)*fibredir(3); 
    for k=1:K
        Fpar=sum(G_dot_fibre(1:k)*tau);
        Bval(m)=Bval(m)+(Fpar.^2)*tau;
    end
end
Bval=GAMMA^2*Bval;
ePar=exp(-Bval.*dRes);

 
 
 
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




    if isfield(protocol,'smalldel') && isfield(protocol,'delta')
        K = floor((protocol.smalldel+1E-10)./tau)+1;
        dstmp = floor((protocol.delta-1E-10)./protocol.tau)-floor((protocol.smalldel+1E-10)./protocol.tau);
        

         %disp('Angle method using exp(i(k-n)theta)')
        G_mag=sqrt(G(:,1:3:end).^2+G(:,2:3:end).^2+ G(:,3:3:end).^2); 
        qhat_X = zeros(size(G_mag));
        qhat_Y = zeros(size(G_mag));
        qhat_Z = zeros(size(G_mag));
        indG=G_mag>1E-6;
        G_tmpX=G(:,1:3:end);
        G_tmpY=G(:,2:3:end);
        G_tmpZ=G(:,3:3:end);
        qhat_X(indG) =  G_tmpX(indG)./G_mag(indG); % qhatx
        qhat_Y(indG) =  G_tmpY(indG)./G_mag(indG); % qhaty
        qhat_Z(indG) =  G_tmpZ(indG)./G_mag(indG); % qhatz

        Gind_mat=round(G_mag./protocol.gstep);   
        if protocol.mirror == 1
            for m=1:M
                      Ecal=Smat';
                      for ind=1:K(m)
                        Gind=Gind_mat(m,ind);
                        if Gind~=0
                          grad_dir = [qhat_X(m,ind); qhat_Y(m,ind); qhat_Z(m,ind)];
                          AmatU=(real(Amat)+(imag(Amat).*kDn)*dot(grad_dir,w))+1i*imag(Amat)*dot(grad_dir,v); 
                          AmatUD=AmatU^Gind;
                          Ecal=Ecal*Rmat*AmatUD;
                        elseif Gind==0
                          Ecal=Ecal*Rmat;
                        end
                      end
                      Rmid=Rmat^dstmp(m);
                      ePerp(m,indj)=Ecal*Rmid*Ecal';                  
             end
        elseif protocol.mirror == 0
            for m=1:M
                      Ecal=Smat';
                      ProdMat2 = eye(size(Rmat));
                      for ind=1:K(m)
                        Gind=Gind_mat(m,ind);
                        if Gind~=0
                         grad_dir = [qhat_X(m,ind); qhat_Y(m,ind); qhat_Z(m,ind)];
                          AmatU=(real(Amat)+(imag(Amat).*kDn)*dot(grad_dir,w))+1i*imag(Amat)*dot(grad_dir,v); 
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
                      ePerp(m,indj)=Ecal*Rmid*ProdMat2*Smat;                 
            end     

        end

        E_r = ePerp.*repmat(ePar,[1,length(protocol.A_cyl) ]);
       
    else
        K=totsizeG/3;
        % Perpendicular
        
        %disp('Angle method using exp(i(k-n)theta)')
       G_mag=sqrt(G(:,1:3:end).^2+G(:,2:3:end).^2+ G(:,3:3:end).^2);   
        indG=G_mag>1E-6;
        qhat_X = zeros(size(G_mag));
        qhat_Y = zeros(size(G_mag));
        qhat_Z = zeros(size(G_mag));
        G_tmpX=G(:,1:3:end);
        G_tmpY=G(:,2:3:end);
        G_tmpZ=G(:,3:3:end);
        qhat_X(indG) =  G_tmpX(indG)./G_mag(indG); % qhatx
        qhat_Y(indG) =  G_tmpY(indG)./G_mag(indG); % qhaty
        qhat_Z(indG) =  G_tmpZ(indG)./G_mag(indG); % qhatz;
        Gind_mat=round(G_mag./protocol.gstep);
        for m=1:M
                  Ecal=Smat';
                  for ind=1:K
                    Gind=Gind_mat(m,ind);
                    if Gind~=0
                      grad_dir = [qhat_X(m,ind); qhat_Y(m,ind); qhat_Z(m,ind)];
                      AmatU=(real(Amat)+(imag(Amat).*kDn)*dot(grad_dir,w))+1i*imag(Amat)*dot(grad_dir,v); 
                      AmatUD=AmatU^Gind;
                      Ecal=Ecal*Rmat*AmatUD;
                    elseif Gind==0
                      Ecal=Ecal*Rmat;
                    end
                  end
                  ePerp(m,indj)=Ecal*Rmat*Smat;
         end


        E_r = ePerp.*repmat(ePar,[1,length(protocol.A_cyl)]);
  
 
    end
end
%Total Signal
E=E_r;
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
        Epert = GammaCylinders_MM_GEN(xpert,protocol);
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
                Epert = GammaCylinders_MM_GEN(xpert,protocol,x_deriv);
                dEtdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEtdx;
            end
        end
         
    end  
    protocol.diff=0;
  
   
end
    
 
