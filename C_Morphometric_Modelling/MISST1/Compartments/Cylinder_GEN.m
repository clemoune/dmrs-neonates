function [E J]=Cylinder_GEN(x,protocol,x_deriv)
% function [E J] = Cylinder_GEN(model, protocol, params_deriv)
% returns the diffusion signal and signal jacobian for the given protocol
% and Cylinder tissue model
%
% Substrate: parallel cylinders with one radius
% Pulse sequence: Generalized gradient spin echo.
% Signal approximation: Matrix Method
%
% x is the list of model parameters in SI units:
% x(1) is the free diffusivity of the material inside and outside the cylinders.
% x(2) is radius of the cylinder
% x(3) is the angle from the z direction
% x(4) is the azymuthal angle from x direction
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

GAMMA = 2.675987E8;
% model parameters
dRes = x(1);
theta = x(3);
phi = x(4);
fibredir = [sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];
tau=protocol.tau;
G=protocol.G;
[M totsizeG]=size(G);

if protocol.diff==0
    Amat=protocol.A_cyl{indj};
    Smat=protocol.S_cyl{indj};
    Rmat=protocol.R_cyl{indj};
elseif protocol.diff==1 % matrices for cumputing derivatives with respect to diffusivity
    Amat=protocol.A_cyl{indj};
    Smat=protocol.S_cyl{indj};
    Rmat=protocol.RpertD_cyl{indj};
elseif protocol.diff==2 % matrices for cumputing derivatives with respect to radius
    Amat=protocol.Aperta_cyl{indj};
    Smat=protocol.Sperta_cyl{indj};
    Rmat=protocol.Rperta_cyl{indj};
else
    error('protocol.diff is not suitable')
end
kDn=protocol.kDn_cyl{indj};

ePerp=zeros(M,1);
cos_theta = cos(theta);
sin_theta = sqrt(1-cos_theta^2);
cos_phi = cos(phi);
sin_phi = sqrt(1-cos_phi^2);


 v = [cos_theta*cos_phi^2+sin_phi^2; -(1-cos_theta)*sin_phi*cos_phi; -sin_theta*cos_phi]; % vectors for the new x and y directions; see Ozarslan 2010
 w = [ -(1-cos_theta)*sin_phi*cos_phi; cos_theta*sin_phi^2+cos_phi^2; -sin_theta*sin_phi];
                         

if isfield(protocol,'smalldel') && isfield(protocol,'delta')
    K = floor((protocol.smalldel+1E-10)./tau)+1;
    dstmp = floor((protocol.delta-1E-10)./protocol.tau)-floor((protocol.smalldel+1E-10)./protocol.tau);
    Bval=zeros(M,1);
    for m=1:M
    G_dot_fibre = G(m,1:3:end)*fibredir(1)+ G(m,2:3:end)*fibredir(2) + G(m,3:3:end)*fibredir(3); 
        for k=1:(2*K(m)+dstmp(m))
            Fpar=sum(G_dot_fibre(1:k)*tau);
            Bval(m)=Bval(m)+(Fpar.^2)*tau;
           
        end
        
        
    end    
    Bval=GAMMA^2*Bval;
    ePar=exp(-Bval.*dRes);
    
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
                  ePerp(m)=Ecal*Rmid*Ecal';                  
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
                  ePerp(m)=Ecal*Rmid*ProdMat2*Smat;                 
        end     
        
    end

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
else
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
              ePerp(m)=Ecal*Rmat*Smat;
     end
    

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
end


% Compute the Jacobian matrix
if(nargout>1)
    J = zeros(length(E),4); % includes fibre direction
    dx = protocol.pert;
    if nargin < 3 
         
        for i = 1:4
        xpert = x;
            if i<=2
            protocol.diff=i;
            else
            protocol.diff=0; % fibre diection
            end
        xpert(i) = xpert(i)*(1+dx);    
        Epert = Cylinder_GEN(xpert,protocol);
        dEtdx = (Epert - E)/(xpert(i)*dx);
         J(:,i) = dEtdx;
        end      

        
    else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
       
        for i = 1:4    
            if x_deriv(i) ~= 0                
                xpert = x;
                    if i<=2
                    protocol.diff=i;
                    else
                    protocol.diff = 0; % derivatives with respect to fibre direction
                    end
                xpert(i) = xpert(i)*(1+dx);    
                Epert = Cylinder_GEN(xpert,protocol,x_deriv);
                dEtdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEtdx;
            end
        end
         
    end  
    protocol.diff=0;
  
   
end
    
 
