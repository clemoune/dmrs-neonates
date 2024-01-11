function [E J]=Cylinder_MM_dPGSE(x,protocol,x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a cylinder compartment.
% 
% [E,J]=Cylinder_MM_dPGSE(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate consisting of parallel, 
% impermeable cylinders with a single radius and a diffusion protocol 
% specified in the input
% Substrate: parallel, impermeable cylinders with a single radius
% Diffusion pulse sequence: Double pulsed gradient spin echo (dPGSE)
% Signal approximation: Matrix method (MM) - propagator expressed via 
%       eigenmode expansion 
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
smalldel1 = protocol.smalldel; % can be extended in the future for protoccols with different timings
smalldel2 = protocol.smalldel;
delta1 = protocol.delta;
delta2 = protocol.delta;
tm = protocol.tm; % mixing time (starts from the begining of the second pulse)

M=length(protocol.smalldel);


grad_dirs1 = protocol.grad_dirs1; % direction of the first gradient
grad_dirs2 = protocol.grad_dirs2; % direction of the second gradient


if any( abs(sqrt(sum(grad_dirs1.^2,2))-1) >1E-5) || any( abs(sqrt(sum(grad_dirs2.^2,2))-1)>1E-5)
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

ePerp=zeros(M,1); % cylinder restriction


G1_par = G1.*(grad_dirs1*fibredir)';
G2_par = G2.*(grad_dirs2*fibredir)';

GAMMA = 2.675987E8; % This is what is used throughout Wuzi.

if ~isfield(protocol,'slew_rate')
    
G_ind1 = round(G1./protocol.gstep);  
G_ind2 = round(G2./protocol.gstep); 
tmp1 = round((smalldel1)./tau);
tmp2 = round((smalldel2)./tau);
dstmp1 =round((delta1-smalldel1)./tau);
dstmp2 =round((delta2-smalldel2)./tau);


if min(tm - smalldel1) < 0 % overlap between the two PGSEs

     error('Not implemented yet in the case tm < smalldel')
else
    bval_par = GAMMA.^2.*(G1_par.^2.*smalldel1.^2.*(delta1-smalldel1/3)+G2_par.^2.*smalldel2.^2.*(delta2-smalldel2/3));


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



         dstmp12 = round((tm(m)-smalldel1(m))./tau);
         Ecal_cyl=Ecal_cyl*(Rmat_cyl*AmatUD1_cyl)^tmp1(m)*Rmat_cyl^dstmp1(m)*(Rmat_cyl*AmatUD1_cyl')^tmp1(m)*... % first PGSE
             Rmat_cyl^dstmp12*(Rmat_cyl*AmatUD2_cyl')^tmp2(m)*Rmat_cyl^dstmp2(m)*(Rmat_cyl*AmatUD2_cyl)^tmp2(m)*real(Smat_cyl); % second PGSE
         ePerp(m) = Ecal_cyl;


    end

end

else
    slew_rate = protocol.slew_rate;
    rt1 = G1./slew_rate;
    rt2 = G2./slew_rate;   

    G_ind1 = floor(G1./protocol.gstep);  
    G_ind2 = floor(G2./protocol.gstep); 
    tmp1 = round((smalldel1)./tau);
    tmp2 = round((smalldel2)./tau);
    dstmp1 =round((delta1-smalldel1)./tau);
    dstmp2 =round((delta2-smalldel2)./tau);
    dstmp12 = round((tm-smalldel1)./tau);
    tmp1_rt = floor(rt1./tau);
    tmp2_rt = floor(rt2./tau);
    tmp1 = tmp1 - 2*tmp1_rt;
    tmp2 = tmp2 - 2*tmp2_rt;

       bval_par = 2*GAMMA.^2.*(G1_par.^2.*smalldel1.^3./15.*(5-15.*rt1./smalldel1./2 - ...
        5.*rt1.^2./smalldel1./4 + 4.*rt1.^3./smalldel1.^3)) +...  
        GAMMA.^2.*G1_par.^2.*(delta1 - smalldel1).*(smalldel1-rt1).^2 + ...
        2*GAMMA.^2.*(G2_par.^2.*smalldel2.^3./15.*(5-15.*rt2./smalldel2./2 - ...
        5.*rt2.^2./smalldel2./4 + 4.*rt2.^3./smalldel2.^3)) +...  
        GAMMA.^2.*G2_par.^2.*(delta2 - smalldel2).*(smalldel2-rt2).^2;

    logEPar=-bval_par'.*dRes;
    ePar = exp(logEPar);
    for m=1:M
         
     
  
       cos_theta = cos(theta);
       sin_theta = sin(theta);
       cos_phi = cos(phi);
       sin_phi = sin(phi);
       v = [cos_theta*cos_phi^2+sin_phi^2; -(1-cos_theta)*sin_phi*cos_phi; -sin_theta*cos_phi]; % vectors for the new x and y directions; see Ozarslan 2010
       w = [ -(1-cos_theta)*sin_phi*cos_phi; cos_theta*sin_phi^2+cos_phi^2; -sin_theta*sin_phi];

       
       
       Ecal_cyl=real(Smat_cyl)'+ imag(Smat_cyl)'*dot(grad_dirs1(m,:)',w)+ 1i*imag(Smat_cyl)'*dot(grad_dirs1(m,:)',v);

       % for the first PGSE
       AmatU1_cyl=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dirs1(m,:)',w))+1i*imag(Amat_cyl)*dot(grad_dirs1(m,:)',v); 
             
       AmatUD1_cyl=AmatU1_cyl^G_ind1(m);
        


       % for the second PGSE
       AmatU2_cyl=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dirs2(m,:)',w))+1i*imag(Amat_cyl)*dot(grad_dirs2(m,:)',v); 
       AmatUD2_cyl=AmatU2_cyl^G_ind2(m);

   
       G1up = linspace(0,G_ind1(m),tmp1_rt(m)+1);
       G1up = G1up(2:end);
       G1down = linspace(G_ind1(m),0,tmp1_rt(m)+1);
       G1down = G1down(2:end);
       
       G2up = linspace(0,G_ind2(m),tmp2_rt(m)+1);
       G2up = G2up(2:end);
       G2down = linspace(G_ind2(m),0,tmp2_rt(m)+1);
       G2down = G2down(2:end);
       
       if length(G1up) >= 1
       
       mat1up_pos = (Rmat_cyl*AmatU1_cyl^G1up(1));
       mat1up_neg = (Rmat_cyl*AmatU1_cyl'^G1up(1));
       mat1down_pos = (Rmat_cyl*AmatU1_cyl^G1down(1));
       mat1down_neg = (Rmat_cyl*AmatU1_cyl'^G1down(1));
       
       for it = 2:length(G1up)       
        mat1up_pos = mat1up_pos*(Rmat_cyl*AmatU1_cyl^G1up(it));
        mat1up_neg = mat1up_neg*(Rmat_cyl*AmatU1_cyl'^G1up(it));
        mat1down_pos = mat1down_pos*(Rmat_cyl*AmatU1_cyl^G1down(it));
        mat1down_neg = mat1down_neg*(Rmat_cyl*AmatU1_cyl'^G1down(it));
       end
       
       else
            mat1up_pos = eye(size(Rmat_cyl));
            mat1up_neg = eye(size(Rmat_cyl));
            mat1down_pos = eye(size(Rmat_cyl));
            mat1down_neg = eye(size(Rmat_cyl));
       end
       
       if length(G2up) >= 1
       
       mat2up_pos = (Rmat_cyl*AmatU2_cyl^G2up(1));
       mat2up_neg = (Rmat_cyl*AmatU2_cyl'^G2up(1));
       mat2down_pos = (Rmat_cyl*AmatU2_cyl^G2down(1));
       mat2down_neg = (Rmat_cyl*AmatU2_cyl'^G2down(1));
       
       for it = 2:length(G1up)       
        mat2up_pos = mat2up_pos*(Rmat_cyl*AmatU2_cyl^G2up(it));
        mat2up_neg = mat2up_neg*(Rmat_cyl*AmatU2_cyl'^G2up(it));
        mat2down_pos = mat2down_pos*(Rmat_cyl*AmatU2_cyl^G2down(it));
        mat2down_neg = mat2down_neg*(Rmat_cyl*AmatU2_cyl'^G2down(it));
       end
       else
            mat2up_pos = eye(size(Rmat_cyl));
            mat2up_neg = eye(size(Rmat_cyl));
            mat2down_pos = eye(size(Rmat_cyl));
            mat2down_neg = eye(size(Rmat_cyl));           
       end
       
       
         Ecal_cyl=Ecal_cyl*mat1up_pos*(Rmat_cyl*AmatUD1_cyl)^tmp1(m)*mat1down_pos*...
              Rmat_cyl^dstmp1(m)*mat1up_neg*(Rmat_cyl*AmatUD1_cyl')^tmp1(m)*mat1down_neg*... % first ODE
             Rmat_cyl^dstmp12(m)*mat2up_neg*(Rmat_cyl*AmatUD2_cyl')^tmp2(m)*mat2down_neg*...
             Rmat_cyl^dstmp2(m)*mat2up_pos*(Rmat_cyl*AmatUD2_cyl)^tmp2(m)*mat2down_pos*real(Smat_cyl); % second ODE


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
        Epert = Cylinder_MM_dPGSE(xpert,protocol);
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
                Epert = Cylinder_MM_dPGSE(xpert,protocol,x_deriv);
                dEtdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEtdx;
            end
        end
         
    end  
    protocol.diff=0;
  
   
end
    
end
