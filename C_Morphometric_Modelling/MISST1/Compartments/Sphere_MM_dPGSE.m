function [E J]=Sphere_MM_dPGSE(x,protocol,x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a sphere compartment.
% 
% [E,J]=Sphere_MM_dPGSE(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate consisting of impermeable 
% spheres with a single radius and a diffusion protocol specified in the input
% Substrate: impermeable spheres with a single radius
% Diffusion pulse sequence: Double pulsed gradient spin echo (dPGSE)
% Signal approximation: Matrix method (MM) - propagator expressed via 
%       eigenmode expansion 
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% x - size 2 vector of model parameters in SI units for Sphere:
%       x(1) - free diffusivity of the material inside the sphere.
%       x(2) - sphere radius   
% protocol - structure which includes all the information related to the 
%        diffusion protocol and required to generate the signal.    
%       protocol.grad_dirs1 - is the gradient direction of the first 
%       gradient pair for each measurement. It has size [N 3] where 
%       N is the number of measurements.
%       protocol.grad_dirs2 - is the gradient direction of the second 
%       gradient pair for each measurement. It has size [N 3] where 
%       N is the number of measurements.
%       protocol.G1 - gradient strength of the first gradient pair, size [1 N]
%       protocol.G2 - gradient strength of the second gradient pair, size [1 N]
%       protocol.delta - pulse separation of each pair, size [1 N]
%       protocol.smalldel - pulse duration of each pair, size [1 N]
%       protocol.tm - mixing time between the two gradient pairs, must be 
%       larger than protocol.smalldel, size [1 N]
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

if ~isfield(protocol,'GENj')
    protocol.GENj=1;
end

indj=protocol.GENj; % this is the index of the radii. 1 if only one.

tau=protocol.tau;
% Find E_restricted
G1=protocol.G1;
G2=protocol.G2;
smalldel1 = protocol.smalldel; % can be extended in the future for protoccols with different timings
smalldel2 = protocol.smalldel;
delta1 = protocol.delta;
delta2 = protocol.delta;
tm = protocol.tm; % mixing time (starts from the begining of the second pulse)

M=length(protocol.smalldel);


grad_dirs1 = protocol.grad_dirs1; % direction of the first gradient
grad_dirs2 = protocol.grad_dirs2; % direction of the second gradient


if protocol.diff==0
    Amat0=protocol.A_sph{indj,1};
    Smat0=protocol.S_sph{indj,1};
    Rmat=protocol.R_sph{indj};
    Amat90=protocol.A_sph{indj,2};
    Smat90=protocol.S_sph{indj,2};
    
elseif protocol.diff==1 % matrices for cumputing derivatives with respect to diffusivity
    Amat0=protocol.A_sph{indj,1};
    Smat0=protocol.S_sph{indj,1};
    Rmat=protocol.RpertD_sph{indj};
    Amat90=protocol.A_sph{indj,2};
    Smat90=protocol.S_sph{indj,2};
    
elseif protocol.diff==2 % matrices for cumputing derivatives with respect to radius
    Amat0=protocol.Aperta_sph{indj,1};
    Smat0=protocol.Sperta_sph{indj,1};
    Rmat=protocol.Rperta_sph{indj};
    Amat90=protocol.Aperta_sph{indj,2};
    Smat90=protocol.Sperta_sph{indj,2};
    
else
    error('protocol.diff is not suitable')
end
kDn=protocol.kDn_sph;

E=zeros(M,1);

G_ind1 = round(G1./protocol.gstep);  
G_ind2 = round(G2./protocol.gstep); 
tmp1 = round((smalldel1)./tau);
tmp2 = round((smalldel2)./tau);
dstmp1 =round((delta1-smalldel1)./tau);
dstmp2 =round((delta2-smalldel2)./tau);

 %disp('Angle method using exp(i(k-n)theta)')

cos_theta1 = zeros(size(G1));
sin_phi1 = zeros(size(cos_theta1));
cos_phi1 = zeros(size(cos_theta1));

cos_theta2 = zeros(size(G1));
sin_phi2 = zeros(size(cos_theta2));
cos_phi2 = zeros(size(cos_theta2));

cos_theta1(G1 >0) = grad_dirs1(G1 >0,3);
sin_theta1 = sqrt(1-cos_theta1.^2);   

cos_theta2(G2 >0) = grad_dirs2(G2 >0,3);
sin_theta2 = sqrt(1-cos_theta2.^2); 

sin_phi1(sin_theta1 ~= 0 & G1 >0) = grad_dirs1(sin_theta1~=0 & G1 >0,2)'./sin_theta1(sin_theta1~=0 & G1 >0);
cos_phi1(sin_theta1 ~= 0 & G1 >0) = grad_dirs1(sin_theta1~=0 & G1 >0,1)'./sin_theta1(sin_theta1~=0 & G1 >0);  
 
sin_phi2(sin_theta2 ~= 0 & G2 >0) = grad_dirs2(sin_theta2~=0 & G2 >0,2)'./sin_theta2(sin_theta2~=0 & G2 >0);
cos_phi2(sin_theta2 ~= 0 & G2 >0) = grad_dirs2(sin_theta2~=0 & G2 >0,1)'./sin_theta2(sin_theta2~=0 & G2 >0); 


if ~isfield(protocol,'slew_rate')
if min(tm - smalldel1) < 0 % overlap between the two PGSEs

     error('Not implemented yet in the case tm < smalldel')
else
    for m=1:M
        Smat = real(Smat0) + imag(Smat90)*sin_phi1(m)+1i*(imag(Smat0)*cos_theta1(m)+imag(Smat90)*cos_phi1(m));  
        Ecal=Smat';
        
        AmatU1=real(Amat0)+(imag(Amat90).*kDn).*sin_phi1(m).*sin_theta1(m)+1i*(imag(Amat0)*cos_theta1(m)+imag(Amat90)*cos_phi1(m)*sin_theta1(m));
        AmatU2=real(Amat0)+(imag(Amat90).*kDn).*sin_phi2(m).*sin_theta2(m)+1i*(imag(Amat0)*cos_theta2(m)+imag(Amat90)*cos_phi2(m)*sin_theta2(m));
        
        AmatUD1=AmatU1^G_ind1(m);
        AmatUD2=AmatU2^G_ind2(m);
        
        dstmp12 = round((tm(m)-smalldel1(m))./tau);
        
         Ecal=Ecal*(Rmat*AmatUD1)^tmp1(m)*Rmat^dstmp1(m)*(Rmat*AmatUD1')^tmp1(m)*... % first PGSE
             Rmat^dstmp12*(Rmat*AmatUD2')^tmp2(m)*Rmat^dstmp2(m)*(Rmat*AmatUD2)^tmp2(m)*Smat; % second PGSE
         E(m) = Ecal;
                           
     end
end
else
       slew_rate = protocol.slew_rate;
    rt1 = G1./slew_rate;
    rt2 = G2./slew_rate;   
       tmp1_rt = floor(rt1./tau);
    tmp2_rt = floor(rt2./tau);
    tmp1 = tmp1 - 2*tmp1_rt;
    tmp2 = tmp2 - 2*tmp2_rt;
if min(tm - smalldel1) < 0 % overlap between the two PGSEs

     error('Not implemented yet in the case tm < smalldel')
else
    for m=1:M
        
        Smat = real(Smat0) + imag(Smat90)*sin_phi1(m)+1i*(imag(Smat0)*cos_theta1(m)+imag(Smat90)*cos_phi1(m));  
        Ecal=Smat';
        
        AmatU1=real(Amat0)+(imag(Amat90).*kDn).*sin_phi1(m).*sin_theta1(m)+1i*(imag(Amat0)*cos_theta1(m)+imag(Amat90)*cos_phi1(m)*sin_theta1(m));
        AmatU2=real(Amat0)+(imag(Amat90).*kDn).*sin_phi2(m).*sin_theta2(m)+1i*(imag(Amat0)*cos_theta2(m)+imag(Amat90)*cos_phi2(m)*sin_theta2(m));
        
        AmatUD1=AmatU1^G_ind1(m);
        AmatUD2=AmatU2^G_ind2(m);
        
        dstmp12 = round((tm(m)-smalldel1(m))./tau);
        
       G1up = linspace(0,G_ind1(m),tmp1_rt(m)+1);
       G1up = G1up(2:end);
       G1down = linspace(G_ind1(m),0,tmp1_rt(m)+1);
       G1down = G1down(2:end);
       
       G2up = linspace(0,G_ind2(m),tmp2_rt(m)+1);
       G2up = G2up(2:end);
       G2down = linspace(G_ind2(m),0,tmp2_rt(m)+1);
       G2down = G2down(2:end);
       
       if length(G1up) >= 1
       
       mat1up_pos = (Rmat*AmatU1^G1up(1));
       mat1up_neg = (Rmat*AmatU1'^G1up(1));
       mat1down_pos = (Rmat*AmatU1^G1down(1));
       mat1down_neg = (Rmat*AmatU1'^G1down(1));
       
       for it = 2:length(G1up)       
        mat1up_pos = mat1up_pos*(Rmat*AmatU1^G1up(it));
        mat1up_neg = mat1up_neg*(Rmat*AmatU1'^G1up(it));
        mat1down_pos = mat1down_pos*(Rmat*AmatU1^G1down(it));
        mat1down_neg = mat1down_neg*(Rmat*AmatU1'^G1down(it));
       end
       
       else
            mat1up_pos = eye(size(Rmat));
            mat1up_neg = eye(size(Rmat));
            mat1down_pos = eye(size(Rmat));
            mat1down_neg = eye(size(Rmat));
       end
       
       if length(G2up) >= 1
       
       mat2up_pos = (Rmat*AmatU2^G2up(1));
       mat2up_neg = (Rmat*AmatU2'^G2up(1));
       mat2down_pos = (Rmat*AmatU2^G2down(1));
       mat2down_neg = (Rmat*AmatU2'^G2down(1));
       
       for it = 2:length(G1up)       
        mat2up_pos = mat2up_pos*(Rmat*AmatU2^G2up(it));
        mat2up_neg = mat2up_neg*(Rmat*AmatU2'^G2up(it));
        mat2down_pos = mat2down_pos*(Rmat*AmatU2^G2down(it));
        mat2down_neg = mat2down_neg*(Rmat*AmatU2'^G2down(it));
       end
       else
            mat2up_pos = eye(size(Rmat));
            mat2up_neg = eye(size(Rmat));
            mat2down_pos = eye(size(Rmat));
            mat2down_neg = eye(size(Rmat));           
       end
        
         Ecal=Ecal*mat1up_pos*(Rmat*AmatUD1)^tmp1(m)*mat1down_pos*Rmat^dstmp1(m)*...
             mat1up_neg*(Rmat*AmatUD1')^tmp1(m)*mat1down_neg*... % first PGSE
             Rmat^dstmp12*mat2up_neg*(Rmat*AmatUD2')^tmp2(m)*mat2down_neg*...
             Rmat^dstmp2(m)*mat2up_pos*(Rmat*AmatUD2)^tmp2(m)*mat2down_pos*real(Smat); % second PGSE
         E(m) = Ecal;
                           
    end
	
end 
end

    if strcmp(protocol.complex,'complex')
      E=[real(mean(E,2));imag(mean(E,2))];
    elseif strcmp(protocol.complex,'real')
      E=real(mean(E,2));
    elseif strcmp(protocol.complex,'abs')
      E=abs(mean(E,2));
    end

% Compute the Jacobian matrix
if(nargout>1)
    J = zeros(length(E),2);
    dx = protocol.pert;
    if nargin < 3 
         
        protocol.diff=1;   
        xpert = x;
        xpert(1) = xpert(1)*(1+dx);    
        Epert = Sphere_MM_dPGSE(xpert,protocol);
        dEtdD = (Epert - E)/(xpert(1)*dx);

        protocol.diff=2;
        xpert = x;
        xpert(2) = xpert(2)*(1+dx);    
        Epert = Sphere_MM_dPGSE(xpert,protocol);
        dEtda = (Epert - E)/(xpert(2)*dx);
        
        J(:,1) = dEtdD;
        J(:,2) = dEtda;
        
    else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
       
        
            if x_deriv(1) ~= 0  
               
                 xpert = x;
                protocol.diff=1;   
                xpert(1) = xpert(1)*(1+dx);    
                Epert = Sphere_MM_dPGSE(xpert,protocol,x_deriv);
                dEtdx = (Epert - E)/(xpert(1)*dx);
                J(:,1) = dEtdx;
            elseif  x_deriv(2) ~= 0  
                
                 xpert = x;
                protocol.diff=2;   
                xpert(2) = xpert(2)*(1+dx);    
                Epert = Sphere_MM_dPGSE(xpert,protocol,x_deriv);
                dEtdx = (Epert - E)/(xpert(2)*dx);
                J(:,2) = dEtdx;
                
            end             

    end  
    protocol.diff=0;
  
   
end