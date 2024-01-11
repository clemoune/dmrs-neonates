function [E J]=GammaCylinders_MM_PGSE(x,protocol,x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a GammaCylinders compartment.
% 
% [E,J]=GammaCylinders_MM_PGSE(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate consisting of parallel
% cylinders with a Gamma distribution of radii and a diffusion 
% protocol specified in the input
% Substrate: Parallel, impermeable cylinders with a Gamma distribution of radii
% Diffusion pulse sequence: Pulsed gradient spin echo (PGSE) 
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
%       protocol.grad_dirs - is the gradient direction of the gradient for
%       each measurement. It has size [N 3] where N is the number of measurements.
%       protocol.G - gradient strength, size [1 N]
%       protocol.delta - pulse separation, size [1 N]
%       protocol.smalldel - pulse duration, size [1 N]
%       protocol.tau - sampling interval of the gradient waveform, required 
%       for MM, size [1 1] 
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


dRes = x(1);
theta = x(4);
phi = x(5);

% Find E_restricted
G=protocol.G;
smalldel = protocol.smalldel; % can be extended in the future for protoccols with different timings
delta = protocol.delta;

M=length(protocol.smalldel);

grad_dirs = protocol.grad_dirs; % direction of the first gradient



if any( abs(sqrt(sum(grad_dirs.^2,2))-1) >1E-5) 
    warning('All gradient orientations must be unit vectors')
end

fibredir = GetFibreOrientation('GammaCylinders',x);


tau = protocol.tau;

G_ind = round(G./protocol.gstep);  
tmp = round((smalldel)./tau);
dstmp =round((delta-smalldel)./tau);

ePerp=zeros(M,length(protocol.A_cyl)); % cylinder restriction


G_par = G.*(grad_dirs*fibredir)';

GAMMA = 2.675987E8; % This is what is used throughout Wuzi.


bval_par = GAMMA.^2.*(G_par.^2.*smalldel.^2.*(delta-smalldel/3));


logEPar=-bval_par'.*dRes;

ePar = exp(logEPar);
ePar = repmat(ePar,[1 length(protocol.A_cyl)]);


for indj = 1:length(protocol.A_cyl)

    % Perpendicular
    if protocol.diff==0
        Amat_cyl=protocol.A_cyl{indj};
        Smat_cyl=protocol.S_cyl{indj};
        Rmat_cyl=protocol.R_cyl{indj};
        Rweight = repmat(protocol.Rweight,[M,1]);

    elseif protocol.diff==1
        Amat_cyl=protocol.A_cyl{indj};
        Smat_cyl=protocol.S_cyl{indj};
        Rmat_cyl=protocol.RpertD_cyl{indj};
        Rweight = repmat(protocol.Rweight,[M,1]);

    elseif protocol.diff==2
        Amat_cyl=protocol.Aperta_cyl{indj};
        Smat_cyl=protocol.Sperta_cyl{indj};
        Rmat_cyl=protocol.Rperta_cyl{indj};
        Rweight = repmat(protocol.Rweight_perta,[M,1]);
        
    elseif protocol.diff==3
        Amat_cyl=protocol.Apertsh_cyl{indj};
        Smat_cyl=protocol.Spertsh_cyl{indj};
        Rmat_cyl=protocol.Rpertsh_cyl{indj};  
        Rweight = repmat(protocol.Rweight_pertsh,[M,1]);

    else
        error('protocol.diff is not suitable')
    end
    kDn_cyl=protocol.kDn_cyl{indj};



  
    % angle_vec = (0:5:175)*pi/180;
    for m=1:M


       cos_theta = cos(theta);
       sin_theta = sin(theta);
       cos_phi = cos(phi);
       sin_phi = sin(phi);
       v = [cos_theta*cos_phi^2+sin_phi^2; -(1-cos_theta)*sin_phi*cos_phi; -sin_theta*cos_phi];
       w = [ -(1-cos_theta)*sin_phi*cos_phi; cos_theta*sin_phi^2+cos_phi^2; -sin_theta*sin_phi];



       Ecal_cyl=real(Smat_cyl)'+ imag(Smat_cyl)'*dot(grad_dirs(m,:)',w)+ 1i*imag(Smat_cyl)'*dot(grad_dirs(m,:)',v);

       % for the first PGSE
       AmatU1_cyl=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dirs(m,:)',w))+1i*imag(Amat_cyl)*dot(grad_dirs(m,:)',v); 
       AmatUD1_cyl=AmatU1_cyl^G_ind(m);




   
         Ecal_cyl=Ecal_cyl*(Rmat_cyl*AmatUD1_cyl)^tmp(m)*Rmat_cyl^dstmp(m)*(Rmat_cyl*AmatUD1_cyl')^tmp(m)*real(Smat_cyl); %  PGSE
         ePerp(m,indj) = Ecal_cyl;

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

% Need to test the use of multiprod for the purpose of speed improvement
%
% Amat_cyl = zeros(size(protocol.A_cyl{1},1),size(protocol.A_cyl{1},2),length(protocol.A_cyl));
% Smat_cyl = zeros(size(protocol.S_cyl{1},1),size(protocol.S_cyl{1},2),length(protocol.S_cyl));
% Rmat_cyl = zeros(size(protocol.R_cyl{1},1),size(protocol.R_cyl{1},2),length(protocol.R_cyl));
% kDn_cyl =  zeros(size(protocol.A_cyl{1},1),size(protocol.A_cyl{1},2),length(protocol.A_cyl));
% 
% 
%     % Perpendicular
%     if protocol.diff==0
%         for indj = 1:length(protocol.A_cyl)
%             Amat_cyl(:,:,indj)=protocol.A_cyl{indj};
%             Smat_cyl(:,:,indj)=protocol.S_cyl{indj};
%             Rmat_cyl(:,:,indj)=protocol.R_cyl{indj};
%             kDn_cyl(:,:,indj)=protocol.kDn_cyl{indj};
%         end
%         Rweight = repmat(protocol.Rweight,[M,1]);
% 
%     elseif protocol.diff==1
%         for indj = 1:length(protocol.A_cyl)
%             Amat_cyl(:,:,indj)=protocol.A_cyl{indj};
%             Smat_cyl(:,:,indj)=protocol.S_cyl{indj};
%             Rmat_cyl(:,:,indj)=protocol.RpertD_cyl{indj};
%             kDn_cyl(:,:,indj)=protocol.kDn_cyl{indj};
%         end
%         Rweight = repmat(protocol.Rweight,[M,1]);
% 
%     elseif protocol.diff==2
%         for indj = 1:length(protocol.A_cyl)
%             Amat_cyl(:,:,indj)=protocol.Aperta_cyl{indj};
%             Smat_cyl(:,:,indj)=protocol.Sperta_cyl{indj};
%             Rmat_cyl(:,:,indj)=protocol.Rperta_cyl{indj};
%             kDn_cyl(:,:,indj)=protocol.kDn_cyl{indj};
%         end
%         Rweight = repmat(protocol.Rweight_perta,[M,1]);
%         
%     elseif protocol.diff==3
%         for indj = 1:length(protocol.A_cyl)
%             Amat_cyl(:,:,indj)=protocol.Apertsh_cyl{indj};
%             Smat_cyl(:,:,indj)=protocol.Spertsh_cyl{indj};
%             Rmat_cyl(:,:,indj)=protocol.Rpertsh_cyl{indj};  
%             kDn_cyl(:,:,indj)=protocol.kDn_cyl{indj};
%         end
%         Rweight = repmat(protocol.Rweight_pertsh,[M,1]);
% 
%     else
%         error('protocol.diff is not suitable')
%     end
% 
% 
% 
% if min(tm - smalldel1) < 0 % overlap between the two PGSEs
% 
%      error('Not implemented yet in the case tm < smalldel')
% else
%   
%     % angle_vec = (0:5:175)*pi/180;
%     for m=1:M
% 
% 
%        cos_theta = cos(theta);
%        sin_theta = sin(theta);
%        cos_phi = cos(phi);
%        sin_phi = sin(phi);
%        v = [cos_theta*cos_phi^2+sin_phi^2; -(1-cos_theta)*sin_phi*cos_phi; -sin_theta*cos_phi];
%        w = [ -(1-cos_theta)*sin_phi*cos_phi; cos_theta*sin_phi^2+cos_phi^2; -sin_theta*sin_phi];
% 
% 
% 
%        Ecal_cyl=permute(real(Smat_cyl),[2,1,3])+ permute(imag(Smat_cyl),[2 1 3])*dot(grad_dirs1(m,:)',w)+ permute(1i*imag(Smat_cyl),[2 1 3])*dot(grad_dirs1(m,:)',v);
% 
%        % for the first PGSE
%        AmatU1_cyl=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dirs1(m,:)',w))+1i*imag(Amat_cyl)*dot(grad_dirs1(m,:)',v); 
%        % for the second PGSE
%        AmatU2_cyl=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dirs2(m,:)',w))+1i*imag(Amat_cyl)*dot(grad_dirs2(m,:)',v); 
%        dstmp12 = round((tm(m)-smalldel1(m))./tau);
%        
%        for indj = 1:size(AmatU2_cyl,3)
%            AmatUD1_cyl=AmatU1_cyl(:,:,indj)^G_ind1(m);       
%            AmatUD2_cyl=AmatU2_cyl(:,:,indj)^G_ind2(m);
%            
%            RA1_pow(:,:,indj) = (Rmat_cyl(:,:,indj)*AmatUD1_cyl)^tmp1(m);
%            RA1t_pow(:,:,indj)= (Rmat_cyl(:,:,indj)*AmatUD1_cyl')^tmp1(m);
%            RA2_pow(:,:,indj) = (Rmat_cyl(:,:,indj)*AmatUD2_cyl)^tmp2(m);
%            RA2t_pow(:,:,indj)= (Rmat_cyl(:,:,indj)*AmatUD2_cyl')^tmp2(m);
%            R1_pow(:,:,indj) = Rmat_cyl(:,:,indj)^dstmp1(m); 
%            R2_pow(:,:,indj) = Rmat_cyl(:,:,indj)^dstmp2(m); 
%            R12_pow(:,:,indj) = Rmat_cyl(:,:,indj)^dstmp12; 
%        end
% 
% 
% 
%          
%          Ecal_cyl=Ecal_cyl*(Rmat_cyl*AmatUD1_cyl)^tmp1(m)*Rmat_cyl^dstmp1(m)*(Rmat_cyl*AmatUD1_cyl')^tmp1(m)*... % first PGSE
%              Rmat_cyl^dstmp12*(Rmat_cyl*AmatUD2_cyl')^tmp2(m)*Rmat_cyl^dstmp2(m)*(Rmat_cyl*AmatUD2_cyl)^tmp2(m)*real(Smat_cyl); % second PGSE
%          ePerp(m,indj) = Ecal_cyl;
% 
%     end
% 
% end








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
        Epert = GammaCylinders_MM_PGSE(xpert,protocol);
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
                Epert = GammaCylinders_MM_PGSE(xpert,protocol,x_deriv);
                dEtdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEtdx;
            end
        end
         
    end  
    protocol.diff=0;
  
   
end
    
end
