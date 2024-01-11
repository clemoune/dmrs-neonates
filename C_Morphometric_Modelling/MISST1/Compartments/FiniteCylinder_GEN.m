function [E J]=FiniteCylinder_GEN(x,protocol,x_deriv)
% Substrate: cylinders with one radius and finite lenghth, oriented along z
% direction
% Pulse sequence: Generalized gradient spin echo.
% Signal approximation: Propagator expressed via eigenmode expansion,
%
% E=FiniteCylinder_GEN(x, grad_dirs, Gdiff, fibredir, roots)
% returns the measurements E according to the model. E is a matrix [N,M].
%
% x is the list of model parameters in SI units:
% x(1) is the free diffusivity of the material inside and outside the cylinders.
% x(2) is radius of the cylinder
% x(3) is the length of the cylinder
% x(4) is the polar angle discribing the fibre direction 
% x(5) is the azimuthal angle discribing the fibre direction 
%
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
theta = x(4);
phi = x(5);
fibredir = [sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];
tau=protocol.tau;
G=protocol.G;
[M totsizeG]=size(G);

% Calculating the parallel & perpendicular signal is similar  (parallel plane & cylindrical restriction )


% Perpendicular
if protocol.diff==0
    Amat_cyl=protocol.A_cyl{indj};
    Smat_cyl=protocol.S_cyl{indj};
    Rmat_cyl=protocol.R_cyl{indj};
    Amat_plane=protocol.A_plane{indj};
    Smat_plane=protocol.S_plane{indj};
    Rmat_plane=protocol.R_plane{indj};    
elseif protocol.diff==1
    Amat_cyl=protocol.A_cyl{indj};
    Smat_cyl=protocol.S_cyl{indj};
    Rmat_cyl=protocol.RpertD_cyl{indj};
    Amat_plane=protocol.A_plane{indj};
    Smat_plane=protocol.S_plane{indj};
    Rmat_plane=protocol.RpertD_plane{indj};
elseif protocol.diff==2
    Amat_cyl=protocol.Aperta_cyl{indj};
    Smat_cyl=protocol.Sperta_cyl{indj};
    Rmat_cyl=protocol.Rperta_cyl{indj};
    Amat_plane=protocol.A_plane{indj};
    Smat_plane=protocol.S_plane{indj};
    Rmat_plane=protocol.R_plane{indj};   
elseif protocol.diff==3
    Amat_cyl=protocol.A_cyl{indj};
    Smat_cyl=protocol.S_cyl{indj};
    Rmat_cyl=protocol.R_cyl{indj};
    Amat_plane=protocol.Apertl_plane{indj};
    Smat_plane=protocol.Spertl_plane{indj};
    Rmat_plane=protocol.Rpertl_plane{indj};   
else
    error('protocol.diff is not suitable')
end
kDn_cyl=protocol.kDn_cyl{indj};
ePerp=zeros(M,1); % cylinder restriction
ePar = zeros(M,1); % parallel planes 
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
          ePerp(m)=Ecal_cyl*Rmid_cyl*Ecal_cyl';   
          Rmid_plane=Rmat_plane^dstmp(m);
          ePar(m)=Ecal_plane*Rmid_plane*Ecal_plane'; 
                   
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
          ePerp(m)=Ecal_cyl*Rmid_cyl*ProdMat2_cyl*Smat_cyl;
          ePar(m)=Ecal_plane*Rmid_plane*ProdMat2_plane*Smat_plane;  
        
        end
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
   
   
    ePerp=zeros(M,1); % cylinder restriction
    ePar = zeros(M,1); % parallel planes  
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
              ePerp(m)=Ecal_cyl*Rmat_cyl*Smat_cyl;
              ePar(m)=Ecal_plane*Rmat_plane*Smat_plane;
           
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

end


% Compute the Jacobian matrix
if(nargout>1)
     J = zeros(length(E),length(x));
    dx = protocol.pert;
    if nargin < 3 
       
        for i = 1:length(x) 
        xpert = x;
            if i <=3 
                protocol.diff=i;
           else
                protocol.diff = 0; % fibre direction
           end
        xpert(i) = xpert(i)*(1+dx);    
        Epert = FiniteCylinder_GEN(xpert,protocol);
        dEtdx = (Epert - E)/(xpert(i)*dx);
        J(:,i) = dEtdx;
        end       
      
    else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
        
        for i = 1:length(x)
            if x_deriv(i) ~= 0                  
                xpert = x;
                if i <=3
                protocol.diff = i;
                else 
                protocol.diff = 0;
                end
                
                xpert(i) = xpert(i)*(1+dx);    
                Epert = FiniteCylinder_GEN(xpert,protocol,x_deriv);
                dEtdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEtdx;
            end
        end
        
         
    end   
    protocol.diff=0;
   
end
    


    
    