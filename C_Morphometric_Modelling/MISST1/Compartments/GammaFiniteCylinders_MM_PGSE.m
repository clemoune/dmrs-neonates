function [E J]=GammaFiniteCylinders_MM_PGSE(x,protocol,x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a GammaFiniteCylinders compartment.
% 
% [E,J]=GammaFiniteCylinders_MM_PGSE(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate consisting of parallel
% finite cylinders with a Gamma distribution of radii and a diffusion 
% protocol specified in the input
% Substrate: Parallel, impermeable finite cylinders with a Gamma
% distribution of radii
% Diffusion pulse sequence: Pulsed gradient spin echo (PGSE) 
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

if ~isfield(protocol,'complex')
  protocol.complex='real';
end
% Model parameters

% dPar=x(2);
% dPerp=x(3);
% dRes=dPar;  % Diffusion ceoff in restricted compartment same as parallel one in hindered.
theta = x(5);
phi = x(6);

% Find E_restricted
G=protocol.G;
smalldel = protocol.smalldel; % can be extended to different durations
delta = protocol.delta;


M=length(protocol.G);


grad_dirs = protocol.grad_dirs; % direction of the first gradient



if any( abs(sqrt(sum(grad_dirs.^2,2))-1) >1E-5) 
    error('All gradient orientations must be unit vectors')
end

tau = protocol.tau;

G_ind = round(G./protocol.gstep);  
tmp = round((smalldel)./tau);

dstmp =round((delta-smalldel)./tau);



ePerp=zeros(M,length(protocol.A_cyl)); % cylinder restriction
ePar = zeros(M,length(protocol.A_cyl)); % parallel planes
fib_dir = GetFibreOrientation('GammaFiniteCylinders',x);

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

    % angle_vec = (0:5:175)*pi/180;
    for m=1:M


           cos_theta = cos(theta);
           sin_theta = sin(theta);
           cos_phi = cos(phi);
           sin_phi = sin(phi);
           v = [cos_theta*cos_phi^2+sin_phi^2; -(1-cos_theta)*sin_phi*cos_phi; -sin_theta*cos_phi];
           w = [ -(1-cos_theta)*sin_phi*cos_phi; cos_theta*sin_phi^2+cos_phi^2; -sin_theta*sin_phi];



           Ecal_cyl=real(Smat_cyl)'+ imag(Smat_cyl)'*dot(grad_dirs(m,:)',w)+ 1i*imag(Smat_cyl)'*dot(grad_dirs(m,:)',v);
           Ecal_plane = real(Smat_plane)'+ 1i* imag(Smat_plane)'*dot(grad_dirs(m,:)',fib_dir);  
           % for the first PGSE
           AmatU_cyl=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dirs(m,:)',w))+1i*imag(Amat_cyl)*dot(grad_dirs(m,:)',v); 
           AmatUD_cyl=AmatU_cyl^G_ind(m);

           AmatU_plane = real(Amat_plane)+ 1i* imag(Amat_plane)*dot(grad_dirs(m,:)',fib_dir);
           AmatUD_plane=AmatU_plane^G_ind(m);

        

        
          
             Ecal_cyl=Ecal_cyl*(Rmat_cyl*AmatUD_cyl)^tmp(m)*Rmat_cyl^dstmp(m)*(Rmat_cyl*AmatUD_cyl')^tmp(m)*real(Smat_cyl); %  PGSE

             Ecal_plane=Ecal_plane*(Rmat_plane*AmatUD_plane)^tmp(m)*Rmat_plane^dstmp(m)*(Rmat_plane*AmatUD_plane')^tmp(m)*real(Smat_plane); % PGSE
             ePerp(m,indj) = Ecal_cyl;
             ePar(m,indj) = Ecal_plane;
          

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
        Epert = GammaFiniteCylinders_MM_PGSE(xpert,protocol);
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
                Epert = GammaFiniteCylinders_MM_PGSE(xpert,protocol,x_deriv);
                dEtdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEtdx;
            end
        end
        
         
    end   
    protocol.diff=0;
   
end
end