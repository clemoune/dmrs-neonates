function [E J]=AstroCylinders_GEN(x,protocol,x_deriv)
% function [E J] = AstroCylinders_GEN(model, protocol, params_deriv)
% returns the diffusion signal and signal jacobian for the given protocol
% and AstroCylinders tissue model
%
% Substrate: astrocylinders with one radius 
% Pulse sequence: Generalized gradient spin echo.
% Signal approximation: Matrix Method
%
% x is the list of model parameters in SI units:
% x(1) is the free diffusivity of the material inside and outside the cylinders.
% x(2) is radius of the cylinder
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
if ~isfield(protocol,'GENj')
    protocol.GENj=1;
end
indj=protocol.GENj; % this is the index of the radii. 1 if only one.
if ~isfield(protocol,'complex')
  protocol.complex='real';
end
% Model parameters
dRes = x(1);
tau=protocol.tau;
G=protocol.G;
[M totsizeG]=size(G);
GAMMA = 2.675987E8; % This is what is used throughout Wuzi.
% Calculating the parallel & perpendicular signal is similar  (parallel plane & cylindrical restriction )


% Perpendicular
if protocol.diff==0
    Amat_cyl=protocol.A_cyl{indj};
    Smat_cyl=protocol.S_cyl{indj};
    Rmat_cyl=protocol.R_cyl{indj};
 
elseif protocol.diff==1
    Amat_cyl=protocol.A_cyl{indj};
    Smat_cyl=protocol.S_cyl{indj};
    Rmat_cyl=protocol.RpertD_cyl{indj};
   
elseif protocol.diff==2
    Amat_cyl=protocol.Aperta_cyl{indj};
    Smat_cyl=protocol.Sperta_cyl{indj};
    Rmat_cyl=protocol.Rperta_cyl{indj};
    
else
    error('protocol.diff is not suitable')
end
kDn_cyl=protocol.kDn_cyl{indj};
 % get Ndir directions on a spehere
Ndir = 50;
fibre_dirs = ReadCaminoElecPS(sprintf('PointSets/Elec%03i.txt',Ndir))';  
[az incl] = cart2sph(fibre_dirs(1,:),fibre_dirs(2,:),fibre_dirs(3,:)); % azimuth and inclinations
theta = pi/2-incl;

sin_theta = sin(theta);
cos_theta = cos(theta);
sin_phi = sin(az);
cos_phi = cos(az);

 v_mat = [cos_theta.*cos_phi.^2+sin_phi.^2; -(1-cos_theta).*sin_phi.*cos_phi; -sin_theta.*cos_phi];
 w_mat = [ -(1-cos_theta).*sin_phi.*cos_phi; cos_theta.*sin_phi.^2+cos_phi.^2; -sin_theta.*sin_phi];
 
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

ePerp=zeros(M,Ndir); % cylinder restriction
ePar = zeros(M,Ndir); % parallel planes 
Bval=zeros(M,Ndir);
 
 
if isfield(protocol,'smalldel') && isfield(protocol,'delta')
    K = floor((protocol.smalldel+1E-10)./tau)+1;
    dstmp = floor((protocol.delta-1E-10)./protocol.tau)-floor((protocol.smalldel+1E-10)./protocol.tau);
    %disp('Angle method using exp(i(k-n)theta)')
    
    
    if protocol.mirror == 1
        for m=1:M
            for n = 1:Ndir
                   Ecal_cyl=Smat_cyl';
                  
                   fib_dir = fibre_dirs(:,n);       
                   w = w_mat(:,n);
                   v = v_mat(:,n);
                   G_dot_fibre = G(m,1:3:end)*fib_dir(1)+ G(m,2:3:end)*fib_dir(2) + G(m,3:3:end)*fib_dir(3); 
                   
                   for k=1:(2*K(m)+dstmp(m))
                        Fpar=sum(G_dot_fibre(1:k)*tau);
                        Bval(m,n)=Bval(m,n)+(Fpar.^2)*tau;
                    end
                                 
                 for ind=1:K(m)
                    Gind=Gind_mat(m,ind);
                    if Gind~=0
                      grad_dir = [qhat_X(m,ind); qhat_Y(m,ind); qhat_Z(m,ind)];
                      AmatU_cyl=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dir,w))+1i*imag(Amat_cyl)*dot(grad_dir,v); 
                      AmatUD_cyl=AmatU_cyl^Gind;                    
                      Ecal_cyl=Ecal_cyl*Rmat_cyl*AmatUD_cyl;                     
                    elseif Gind==0
                     Ecal_cyl=Ecal_cyl*Rmat_cyl;                     
                    end
                  end
                  Rmid_cyl=Rmat_cyl^dstmp(m);                  
                  ePerp(m,n)=Ecal_cyl*Rmid_cyl*Ecal_cyl';
                 
            end
        end
        Bval=GAMMA^2*Bval;
        ePar=exp(-Bval.*dRes);
    elseif protocol.mirror == 0
        for m=1:M
              for n = 1:Ndir
                   Ecal_cyl=Smat_cyl';
                
                  fib_dir = fibre_dirs(:,n);       
                   w = w_mat(:,n);
                   v = v_mat(:,n);
                   ProdMat2_cyl = eye(size(Amat_cyl)); 
                   G_dot_fibre = G(m,1:3:end)*fib_dir(1)+ G(m,2:3:end)*fib_dir(2) + G(m,3:3:end)*fib_dir(3); 
                   for k=1:(2*K(m)+dstmp(m))
                        Fpar=sum(G_dot_fibre(1:k)*tau);
                        Bval(m,n)=Bval(m,n)+(Fpar.^2)*tau;
                    end
                 
                  for ind=1:K(m)
                      
                     Gind=Gind_mat(m,ind);
                    if Gind~=0
                      grad_dir = [qhat_X(m,ind); qhat_Y(m,ind); qhat_Z(m,ind)];
                      AmatU_cyl=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dir,w))+1i*imag(Amat_cyl)*dot(grad_dir,v); 
                      AmatUD_cyl=AmatU_cyl^Gind;
                      AmatUDT_cyl = AmatUD_cyl';                     
                      Ecal_cyl=Ecal_cyl*Rmat_cyl*AmatUD_cyl;                      
                      ProdMat2_cyl = ProdMat2_cyl*Rmat_cyl*AmatUDT_cyl;                       
                    elseif Gind==0
                     Ecal_cyl=Ecal_cyl*Rmat_cyl;                     
                     ProdMat2_cyl = ProdMat2_cyl*Rmat_cyl;                    
                    end
                  end
                  Rmid_cyl=Rmat_cyl^dstmp(m);
                  ePerp(m,n)=Ecal_cyl*Rmid_cyl*ProdMat2_cyl*Smat_cyl;                                             
              end         
        end
         Bval=GAMMA^2*Bval;
        ePar=exp(-Bval.*dRes);
    end
    E = ePar.*ePerp;
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
  
        for m=1:M
            for n = 1:Ndir
                   Ecal_cyl=Smat_cyl';
                   fib_dir = fibre_dirs(:,n);       
                   w = w_mat(:,n);
                   v = v_mat(:,n);
                   G_dot_fibre = G(m,1:3:end)*fib_dir(1)+ G(m,2:3:end)*fib_dir(2) + G(m,3:3:end)*fib_dir(3); 
                 
                  for ind=1:K
                     Fpar=sum(G_dot_fibre(1:ind)*tau);
                     Bval(m,n)=Bval(m,n)+(Fpar.^2)*tau; 
                     Gind=Gind_mat(m,ind);
                    if Gind~=0
                      grad_dir = [qhat_X(m,ind); qhat_Y(m,ind); qhat_Z(m,ind)];
                      AmatU_cyl=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dir,w))+1i*imag(Amat_cyl)*dot(grad_dir,v); 
                      AmatUD_cyl=AmatU_cyl^Gind;
                    
                      Ecal_cyl=Ecal_cyl*Rmat_cyl*AmatUD_cyl;
                      
                    elseif Gind==0
                     Ecal_cyl=Ecal_cyl*Rmat_cyl;
                     
                    end
                  end
                  ePerp(m,n)=Ecal_cyl*Rmat_cyl*Smat_cyl;
                 
            end
        end
        Bval=GAMMA^2*Bval;
        ePar=exp(-Bval.*dRes);
        E = ePar.*ePerp;
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
    J = zeros(length(E),2);
    dx = protocol.pert;
    if nargin < 3 
         
        xpert = x;
        protocol.diff=1;   
        xpert(1) = xpert(1)*(1+dx);    
        Epert = AstroCylinders_GEN(xpert,protocol);
        dEtdD = (Epert - E)/(xpert(1)*dx);

        protocol.diff=2;
        xpert = x;
        xpert(2) = xpert(2)*(1+dx);    
        Epert = AstroCylinders_GEN(xpert,protocol);
        dEtda = (Epert - E)/(xpert(2)*dx);
        
        J(:,1) = dEtdD;
        J(:,2) = dEtda;
       
    else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
      
            if x_deriv(1) ~= 0  
                
                xpert = x;
                protocol.diff=1;   
                xpert(1) = xpert(1)*(1+dx);    
                Epert = AstroCylinders_GEN(xpert,protocol,x_deriv);
                dEtdx = (Epert - E)/(xpert(1)*dx);
                J(:,1) = dEtdx;
            elseif  x_deriv(2) ~= 0                 
                 xpert = x;
                protocol.diff=2;   
                xpert(2) = xpert(2)*(1+dx);    
                Epert = AstroCylinders_GEN(xpert,protocol,x_deriv);
                dEtdx = (Epert - E)/(xpert(2)*dx);
                J(:,2) = dEtdx;
                
            end             
        
    end   
    protocol.diff=0;
   
end
    