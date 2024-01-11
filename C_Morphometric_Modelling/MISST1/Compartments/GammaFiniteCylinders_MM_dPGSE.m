function [E J]=GammaFiniteCylinders_MM_dPGSE(x,protocol,x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating the GammaFiniteCylinders compartment.
% 
% [E,J]=GammaFiniteCylinders_MM_dPGSE(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate consisting of parallel
% finite cylinders with a Gamma distribution of radii and a diffusion 
% protocol specified in the input
% Substrate: Parallel, impermeable finite cylinders with a Gamma
% distribution of radii
% Diffusion pulse sequence: Double pulsed gradient spin echo (dPGSE)
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
G1=protocol.G1;
G2=protocol.G2;
smalldel1 = protocol.smalldel; % can be extended to different durations
smalldel2 = protocol.smalldel;
delta1 = protocol.delta;
delta2 = protocol.delta;
tm = protocol.tm; % mixing time (starts from the begining of the second pulse)

M=length(protocol.G1);


grad_dirs1 = protocol.grad_dirs1; % direction of the first gradient
grad_dirs2 = protocol.grad_dirs2; % direction of the second gradient


if any( abs(sqrt(sum(grad_dirs1.^2,2))-1) >1E-5) || any( abs(sqrt(sum(grad_dirs2.^2,2))-1)>1E-5)
    error('All gradient orientations must be unit vectors')
end

tau = protocol.tau;

G_ind1 = round(G1./protocol.gstep);  
G_ind2 = round(G2./protocol.gstep); 
tmp1 = round((smalldel1)./tau);
tmp2 = round((smalldel2)./tau);
dstmp1 =round((delta1-smalldel1)./tau);
dstmp2 =round((delta2-smalldel2)./tau);


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



           Ecal_cyl=real(Smat_cyl)'+ imag(Smat_cyl)'*dot(grad_dirs1(m,:)',w)+ 1i*imag(Smat_cyl)'*dot(grad_dirs1(m,:)',v);
           Ecal_plane = real(Smat_plane)'+ 1i* imag(Smat_plane)'*dot(grad_dirs1(m,:)',fib_dir);  
           % for the first PGSE
           AmatU1_cyl=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dirs1(m,:)',w))+1i*imag(Amat_cyl)*dot(grad_dirs1(m,:)',v); 
           AmatUD1_cyl=AmatU1_cyl^G_ind1(m);

           AmatU1_plane = real(Amat_plane)+ 1i* imag(Amat_plane)*dot(grad_dirs1(m,:)',fib_dir);
           AmatUD1_plane=AmatU1_plane^G_ind1(m);

           % for the second PGSE
           AmatU2_cyl=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot(grad_dirs2(m,:)',w))+1i*imag(Amat_cyl)*dot(grad_dirs2(m,:)',v); 
           AmatUD2_cyl=AmatU2_cyl^G_ind2(m);

           AmatU2_plane = real(Amat_plane)+ 1i* imag(Amat_plane)*dot(grad_dirs2(m,:)',fib_dir);
           AmatUD2_plane=AmatU2_plane^G_ind2(m);

             if tm(m) >= smalldel1(m) % no overlap between the two PGSEs
                 dstmp12 = round((tm(m)-smalldel1(m))./tau);
                 Ecal_cyl=Ecal_cyl*(Rmat_cyl*AmatUD1_cyl)^tmp1(m)*Rmat_cyl^dstmp1(m)*(Rmat_cyl*AmatUD1_cyl')^tmp1(m)*... % first PGSE
                     Rmat_cyl^dstmp12*(Rmat_cyl*AmatUD2_cyl')^tmp2(m)*Rmat_cyl^dstmp2(m)*(Rmat_cyl*AmatUD2_cyl)^tmp2(m)*real(Smat_cyl); % second PGSE

                 Ecal_plane=Ecal_plane*(Rmat_plane*AmatUD1_plane)^tmp1(m)*Rmat_plane^dstmp1(m)*(Rmat_plane*AmatUD1_plane')^tmp1(m)*... % first PGSE
                     Rmat_plane^dstmp12*(Rmat_plane*AmatUD2_plane')^tmp2(m)*Rmat_plane^dstmp2(m)*(Rmat_plane*AmatUD2_plane)^tmp2(m)*real(Smat_plane); % second PGSE
                 ePerp(m,indj) = Ecal_cyl;
                 ePar(m,indj) = Ecal_plane;
             else
                 error('Not implemented yet in the case tm < smalldel')
             end

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
        Epert = GammaFiniteCylinders_MM_dPGSE(xpert,protocol);
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
                Epert = GammaFiniteCylinders_MM_dPGSE(xpert,protocol,x_deriv);
                dEtdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEtdx;
            end
        end
        
         
    end   
    protocol.diff=0;
   
end
end
