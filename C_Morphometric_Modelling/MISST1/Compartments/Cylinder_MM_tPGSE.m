function [E J]=Cylinder_MM_tPGSE(x,protocol,x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a cylinder compartment.
% 
% [E,J]=Cylinder_GPD_tPGSE(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate consisting of parallel, 
% impermeable cylinders with a single radius and a diffusion protocol 
% specified in the input
% Substrate: parallel, impermeable cylinders with a single radius
% Diffusion pulse sequence: Triple pulsed gradient (tPGSE) 
% Signal approximation: Gaussian Phase Distribution (GPD)  
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% x - size 4 vector of model parameters in SI units for Cylinder:
%       x(1) - free diffusivity of the material inside the cylinders
%       x(2) - cylinder radius
%       x(3) - polar angle theta in spherical coordinates desbribing the fibre
% direction
%       x(4) - azimuthal angle phi in spherical coordinates describing the
% fibre direction
% protocol - structure which includes all the information related to the 
%        diffusion protocol and required to generate the signal.    
%       protocol.grad_dirs1, protocol.grad_dirs2, protocol.grad_dirs3 - 
%       are the gradient directions of the first, second and third gradient pairs
%       for each measurement. Each has size [N 3] where N is the number of measurements.
%       protocol.G1, protocol.G2, protocol.G3  - gradient strengths of the
%       first, second and third gradient pairs. Each has size [1 N]
%       protocol.delta - pulse separation of each pair, size [1 N]
%       protocol.smalldel - pulse duration of each pair, size [1 N]
%       protocol.tm - mixing time between the first and second pairs and 
%       between the second and third pairs, size [1 N]
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
%   Daniel C. Alexander (d.alexander@ucl.ac.uk)
% 

if ~isfield(protocol,'GENj')
    protocol.GENj=1;
end
indj=protocol.GENj; % this is the index of the radii. 1 if only one.
if ~isfield(protocol,'complex')
  protocol.complex='real';
end
% Model parameters

% dPar=x(2);
% dPerp=x(3);
% dRes=dPar;  % Diffusion ceoff in restricted compartment same as parallel one in hindered.

dRes = x(1);
theta = x(3);
phi = x(4);

% Find E_restricted
G1=protocol.G1;
G2=protocol.G2;
G3=protocol.G3;
smalldel1 = protocol.smalldel; % can be extended in the future for protoccols with different timings
smalldel2 = protocol.smalldel;
smalldel3 = protocol.smalldel;
delta1 = protocol.delta;
delta2 = protocol.delta;
delta3 = protocol.delta;
tm12 = protocol.tm; % mixing time (starts from the begining of the second pulse)
tm23 = protocol.tm;

M=length(protocol.smalldel);


grad_dirs1 = protocol.grad_dirs1; % direction of the first gradient
grad_dirs2 = protocol.grad_dirs2; % direction of the second gradient
grad_dirs3 = protocol.grad_dirs3; % direction of the second gradient

if any( abs(sqrt(sum(grad_dirs1.^2,2))-1) >1E-5) || any( abs(sqrt(sum(grad_dirs2.^2,2))-1)>1E-5) || any( abs(sqrt(sum(grad_dirs3.^2,2))-1)>1E-5)
    error('All gradient orientations must be unit vectors')
end



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
fibredir = GetFibreOrientation('Cylinder',x);


tau = protocol.tau;

G_ind1 = round(G1./protocol.gstep);  
G_ind2 = round(G2./protocol.gstep); 
G_ind3 = round(G3./protocol.gstep); 
tmp1 = round((smalldel1)./tau);
tmp2 = round((smalldel2)./tau);
tmp3 = round((smalldel3)./tau);
dstmp1 =round((delta1-smalldel1)./tau);
dstmp2 =round((delta2-smalldel2)./tau);
dstmp3 =round((delta3-smalldel3)./tau);

ePerp=zeros(M,1); % cylinder restriction


G1_par = G1.*(grad_dirs1*fibredir)';
G2_par = G2.*(grad_dirs2*fibredir)';
G3_par = G3.*(grad_dirs3*fibredir)';

GAMMA = 2.675987E8; % This is what is used throughout Wuzi.

if min(tm12 - smalldel1) || min(tm23 - smalldel2)   < 0 % overlap between the two PGSEs

     error('Not implemented yet in the case tm < smalldel')
else
    bval_par = GAMMA.^2.*(G1_par.^2.*smalldel1.^2.*(delta1-smalldel1/3)+G2_par.^2.*smalldel2.^2.*(delta2-smalldel2/3)+...
        G3_par.^2.*smalldel3.^2.*(delta3-smalldel3/3));


    logEPar=-bval_par'.*dRes;

    ePar = exp(logEPar);

    % angle_vec = (0:5:175)*pi/180;
    for m=1:M


       cos_theta = cos(theta);
       sin_theta = sin(theta);
       cos_phi = cos(phi);
       sin_phi = sin(phi);
       v = [cos_theta*cos_phi^2+sin_phi^2; -(1-cos_theta)*sin_phi*cos_phi; -sin_theta*cos_phi];
       w = [ -(1-cos_theta)*sin_phi*cos_phi; cos_theta*sin_phi^2+cos_phi^2; -sin_theta*sin_phi];



       Ecal_cyl=real(Smat_cyl)'+ imag(Smat_cyl)'*dot(grad_dirs1(m,:)',w)+ 1i*imag(Smat_cyl)'*dot(grad_dirs1(m,:)',v);

       % for the first PGSE
       AmatU1_cyl=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dirs1(m,:)',w))+1i*imag(Amat_cyl)*dot(grad_dirs1(m,:)',v); 
       AmatUD1_cyl=AmatU1_cyl^G_ind1(m);


       % for the second PGSE
       AmatU2_cyl=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dirs2(m,:)',w))+1i*imag(Amat_cyl)*dot(grad_dirs2(m,:)',v); 
       AmatUD2_cyl=AmatU2_cyl^G_ind2(m);
       
       % for the third PGSE
       AmatU3_cyl=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dirs3(m,:)',w))+1i*imag(Amat_cyl)*dot(grad_dirs3(m,:)',v); 
       AmatUD3_cyl=AmatU3_cyl^G_ind3(m);



         dstmp12 = round((tm12(m)-smalldel1(m))./tau);
         dstmp23 = round((tm23(m)-smalldel2(m))./tau);
         
         Ecal_cyl=Ecal_cyl*(Rmat_cyl*AmatUD1_cyl)^tmp1(m)*Rmat_cyl^dstmp1(m)*(Rmat_cyl*AmatUD1_cyl')^tmp1(m)*... % first PGSE
             Rmat_cyl^dstmp12*(Rmat_cyl*AmatUD2_cyl')^tmp2(m)*Rmat_cyl^dstmp2(m)*(Rmat_cyl*AmatUD2_cyl)^tmp2(m)*...
              Rmat_cyl^dstmp23*(Rmat_cyl*AmatUD3_cyl')^tmp3(m)*Rmat_cyl^dstmp3(m)*(Rmat_cyl*AmatUD3_cyl)^tmp3(m)*real(Smat_cyl); % second PGSE
         ePerp(m) = Ecal_cyl;


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


% Compute the Jacobian matrix
if(nargout>1)
    J = zeros(length(E),length(x)); % includes fibre direction
    dx = protocol.pert;
    if nargin < 3 
         
        for i = 1:length(x)
        xpert = x;
            if i<=2
            protocol.diff=i;
            else
            protocol.diff=0; % fibre diection
            end
        xpert(i) = xpert(i)*(1+dx);    
        Epert = Cylinder_MM_tPGSE(xpert,protocol);
        dEtdx = (Epert - E)/(xpert(i)*dx);
         J(:,i) = dEtdx;
        end      

        
    else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
       
        for i = 1:length(x)    
            if x_deriv(i) ~= 0                
                xpert = x;
                    if i<=2
                    protocol.diff=i;
                    else
                    protocol.diff = 0; % derivatives with respect to fibre direction
                    end
                xpert(i) = xpert(i)*(1+dx);    
                Epert = Cylinder_MM_tPGSE(xpert,protocol,x_deriv);
                dEtdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEtdx;
            end
        end
         
    end  
    protocol.diff=0;
  
   
end
    
end