function [E J]=Cylinder_MM_PGSE(x,protocol,x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a cylinder compartment.
% 
% [E,J]=Cylinder_MM_PGSE(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate consisting of parallel, 
% impermeable cylinders with a single radius and a diffusion protocol 
% specified in the input
% Substrate: parallel, impermeable cylinders with a single radius
% Diffusion pulse sequence: Pulsed gradient spin echo (PGSE) 
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
G=protocol.G;
smalldel = protocol.smalldel; % can be extended in the future for protoccols with different timings
delta = protocol.delta;

M=length(protocol.smalldel);


grad_dirs = protocol.grad_dirs; % direction of the first gradient


if any( abs(sqrt(sum(grad_dirs.^2,2))-1) >1E-5)
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

G_ind = round(G./protocol.gstep);  
tmp = round((smalldel)./tau);
dstmp =round((delta-smalldel)./tau);


ePerp=zeros(M,1); % cylinder restriction


G_par = G.*(grad_dirs*fibredir)';

GAMMA = 2.675987E8; % This is what is used throughout Wuzi.

    bval_par = GAMMA.^2.*(G_par.^2.*smalldel.^2.*(delta-smalldel/3));


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



       Ecal_cyl=real(Smat_cyl)'+ imag(Smat_cyl)'*dot(grad_dirs(m,:)',w)+ 1i*imag(Smat_cyl)'*dot(grad_dirs(m,:)',v);

       % for the first PGSE
       AmatU_cyl=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dirs(m,:)',w))+1i*imag(Amat_cyl)*dot(grad_dirs(m,:)',v); 
       AmatUD_cyl=AmatU_cyl^G_ind(m);

   
         Ecal_cyl=Ecal_cyl*(Rmat_cyl*AmatUD_cyl)^tmp(m)*Rmat_cyl^dstmp(m)*(Rmat_cyl*AmatUD_cyl')^tmp(m)*real(Smat_cyl); %  PGSE
         ePerp(m) = Ecal_cyl;


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
        Epert = Cylinder_MM_PGSE(xpert,protocol);
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
                Epert = Cylinder_MM_PGSE(xpert,protocol,x_deriv);
                dEtdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEtdx;
            end
        end
         
    end  
    protocol.diff=0;
  
   
end
    
end
