function [E J]=FiniteCylinder_MM_dPGSE(x,protocol,x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a finite cylinder compartment.
% 
% [E,J]=FiniteCylinder_MM_dPGSE(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate consisting of parallel, 
% impermeable finite cylinders with a single radius and a diffusion protocol 
% specified in the input
% Substrate: parallel, impermeable finite cylinders with a single radius
% Diffusion pulse sequence: Double pulsed gradient spin echo (dPGSE)
% Signal approximation: Matrix method (MM) - propagator expressed via 
%       eigenmode expansion
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% x - size 5 vector of model parameters in SI units for FiniteCylinder:
%       x(1) - free diffusivity of the material inside the cylinders.
%       x(2) - radius of the cylinders.
%       x(3) - eccentricity (ratio between length and diameter)
%       x(4) - polar angle theta in spherical coordinates desbribing the fibre
% direction
%       x(5) - azimuthal angle phi in spherical coordinates describing the
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
theta = x(4);
phi = x(5);

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
%     Amat_plane=protocol.Aperta_plane{indj};
%     Smat_plane=protocol.Sperta_plane{indj};
%     Rmat_plane=protocol.Rperta_plane{indj};  
elseif protocol.diff==3
    Amat_cyl=protocol.A_cyl{indj};
    Smat_cyl=protocol.S_cyl{indj};
    Rmat_cyl=protocol.R_cyl{indj};
    Amat_plane=protocol.ApertEcc_plane{indj};
    Smat_plane=protocol.SpertEcc_plane{indj};
    Rmat_plane=protocol.RpertEcc_plane{indj};   
else
    error('protocol.diff is not suitable')
end
kDn_cyl=protocol.kDn_cyl{indj};



tau = protocol.tau;

if ~isfield(protocol,'slew_rate')

G_ind1 = round(G1./protocol.gstep);  
G_ind2 = round(G2./protocol.gstep); 
tmp1 = round((smalldel1)./tau);
tmp2 = round((smalldel2)./tau);
dstmp1 =round((delta1-smalldel1)./tau);
dstmp2 =round((delta2-smalldel2)./tau);


ePerp=zeros(M,1); % cylinder restriction
ePar = zeros(M,1); % parallel planes
fib_dir = GetFibreOrientation('FiniteCylinder',x);

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
             ePerp(m) = Ecal_cyl;
             ePar(m) = Ecal_plane;
         elseif tm(m) == 0 && tmp1(m) == tmp2(m)
              % for the intersection of the two 
       AmatU3_cyl=(real(Amat_cyl)+(imag(Amat_cyl).*kDn_cyl)*dot((-grad_dirs1(m,:)-grad_dirs2(m,:))',w))+1i*imag(Amat_cyl)*dot((-grad_dirs1(m,:)-grad_dirs2(m,:))',v); 
       AmatUD3_cyl=AmatU3_cyl^G_ind2(m);

       AmatU3_plane = real(Amat_plane)+ 1i* imag(Amat_plane)*dot((-grad_dirs1(m,:)-grad_dirs2(m,:))',fib_dir);
       AmatUD3_plane=AmatU3_plane^G_ind2(m);
       
       Ecal_cyl=Ecal_cyl*(Rmat_cyl*AmatUD1_cyl)^tmp1(m)*Rmat_cyl^dstmp1(m)*... % first PGSE
                 (Rmat_cyl*AmatUD3_cyl)^tmp1(m)*Rmat_cyl^dstmp2(m)*(Rmat_cyl*AmatUD2_cyl)^tmp2(m)*real(Smat_cyl); % second PGSE
            
             Ecal_plane=Ecal_plane*(Rmat_plane*AmatUD1_plane)^tmp1(m)*Rmat_plane^dstmp1(m)*... % first PGSE
                 (Rmat_plane*AmatUD3_plane)^tmp1(m)*Rmat_plane^dstmp2(m)*(Rmat_plane*AmatUD2_plane)^tmp2(m)*real(Smat_plane); % second PGSE
             ePerp(m) = Ecal_cyl;
             ePar(m) = Ecal_plane;
             
         else
             error('Not implemented yet in the case tm < smalldel')
         end
          
end
else
    slew_rate = protocol.slew_rate;
    rt1 = G1./slew_rate;
    rt2 = G2./slew_rate; 
    G_ind1 = round(G1./protocol.gstep);  
    G_ind2 = round(G2./protocol.gstep); 
    tmp1 = round((smalldel1)./tau);
    tmp2 = round((smalldel2)./tau);
    dstmp1 =round((delta1-smalldel1)./tau);
    dstmp2 =round((delta2-smalldel2)./tau);
    dstmp12 = round((tm-smalldel1)./tau);
    tmp1_rt = floor(rt1./tau);
    tmp2_rt = floor(rt2./tau);
    tmp1 = tmp1 - 2*tmp1_rt;
    tmp2 = tmp2 - 2*tmp2_rt;

    ePerp=zeros(M,1); % cylinder restriction
    ePar = zeros(M,1); % parallel planes
    fib_dir = GetFibreOrientation('FiniteCylinder',x);
    
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

           G1up = linspace(0,G_ind1(m),tmp1_rt(m)+1);
           G1up = G1up(2:end);
           G1down = linspace(G_ind1(m),0,tmp1_rt(m)+1);
           G1down = G1down(2:end);

           G2up = linspace(0,G_ind2(m),tmp2_rt(m)+1);
           G2up = G2up(2:end);
           G2down = linspace(G_ind2(m),0,tmp2_rt(m)+1);
           G2down = G2down(2:end);
           
           
       if length(G1up) >= 1    
           
           mat1up_pos_cyl = (Rmat_cyl*AmatU1_cyl^G1up(1));
           mat1up_neg_cyl = (Rmat_cyl*AmatU1_cyl'^G1up(1));
           mat1down_pos_cyl = (Rmat_cyl*AmatU1_cyl^G1down(1));
           mat1down_neg_cyl = (Rmat_cyl*AmatU1_cyl'^G1down(1));

           for it = 2:length(G1up)       
            mat1up_pos_cyl = mat1up_pos_cyl*(Rmat_cyl*AmatU1_cyl^G1up(it));
            mat1up_neg_cyl = mat1up_neg_cyl*(Rmat_cyl*AmatU1_cyl'^G1up(it));
            mat1down_pos_cyl = mat1down_pos_cyl*(Rmat_cyl*AmatU1_cyl^G1down(it));
            mat1down_neg_cyl = mat1down_neg_cyl*(Rmat_cyl*AmatU1_cyl'^G1down(it));
           end

             mat1up_pos_plane = (Rmat_plane*AmatU1_plane^G1up(1));
           mat1up_neg_plane = (Rmat_plane*AmatU1_plane'^G1up(1));
           mat1down_pos_plane = (Rmat_plane*AmatU1_plane^G1down(1));
           mat1down_neg_plane = (Rmat_plane*AmatU1_plane'^G1down(1));

           for it = 2:length(G1up)       
            mat1up_pos_plane = mat1up_pos_plane*(Rmat_plane*AmatU1_plane^G1up(it));
            mat1up_neg_plane = mat1up_neg_plane*(Rmat_plane*AmatU1_plane'^G1up(it));
            mat1down_pos_plane = mat1down_pos_plane*(Rmat_plane*AmatU1_plane^G1down(it));
            mat1down_neg_plane = mat1down_neg_plane*(Rmat_plane*AmatU1_plane'^G1down(it));
           end
       else
            mat1up_pos_cyl = eye(size(Rmat_cyl));
            mat1up_neg_cyl = eye(size(Rmat_cyl));
            mat1down_pos_cyl = eye(size(Rmat_cyl));
            mat1down_neg_cyl = eye(size(Rmat_cyl));
            
            mat1up_pos_plane = eye(size(Rmat_plane));
            mat1up_neg_plane = eye(size(Rmat_plane));
            mat1down_pos_plane = eye(size(Rmat_plane));
            mat1down_neg_plane = eye(size(Rmat_plane));
       end
       
       if length(G2up) >= 1  
       
           mat2up_pos_cyl = (Rmat_cyl*AmatU2_cyl^G2up(1));
           mat2up_neg_cyl = (Rmat_cyl*AmatU2_cyl'^G2up(1));
           mat2down_pos_cyl = (Rmat_cyl*AmatU2_cyl^G2down(1));
           mat2down_neg_cyl = (Rmat_cyl*AmatU2_cyl'^G2down(1));

           for it = 2:length(G1up)       
            mat2up_pos_cyl = mat2up_pos_cyl*(Rmat_cyl*AmatU2_cyl^G2up(it));
            mat2up_neg_cyl = mat2up_neg_cyl*(Rmat_cyl*AmatU2_cyl'^G2up(it));
            mat2down_pos_cyl = mat2down_pos_cyl*(Rmat_cyl*AmatU2_cyl^G2down(it));
            mat2down_neg_cyl = mat2down_neg_cyl*(Rmat_cyl*AmatU2_cyl'^G2down(it));
           end    

           mat2up_pos_plane = (Rmat_plane*AmatU2_plane^G2up(1));
           mat2up_neg_plane = (Rmat_plane*AmatU2_plane'^G2up(1));
           mat2down_pos_plane = (Rmat_plane*AmatU2_plane^G2down(1));
           mat2down_neg_plane = (Rmat_plane*AmatU2_plane'^G2down(1));

           for it = 2:length(G1up)       
            mat2up_pos_plane = mat2up_pos_plane*(Rmat_plane*AmatU2_plane^G2up(it));
            mat2up_neg_plane = mat2up_neg_plane*(Rmat_plane*AmatU2_plane'^G2up(it));
            mat2down_pos_plane = mat2down_pos_plane*(Rmat_plane*AmatU2_plane^G2down(it));
            mat2down_neg_plane = mat2down_neg_plane*(Rmat_plane*AmatU2_plane'^G2down(it));
           end
           
       else
            mat2up_pos_cyl = eye(size(Rmat_cyl));
            mat2up_neg_cyl = eye(size(Rmat_cyl));
            mat2down_pos_cyl = eye(size(Rmat_cyl));
            mat2down_neg_cyl = eye(size(Rmat_cyl));
            
            mat2up_pos_plane = eye(size(Rmat_plane));
            mat2up_neg_plane = eye(size(Rmat_plane));
            mat2down_pos_plane = eye(size(Rmat_plane));
            mat2down_neg_plane = eye(size(Rmat_plane));
       end    
         if tm(m) >= smalldel1(m) % no overlap between the two PGSEs    
           Ecal_cyl=Ecal_cyl* mat1up_pos_cyl*(Rmat_cyl*AmatUD1_cyl)^tmp1(m)*mat1down_pos_cyl*...
             Rmat_cyl^dstmp1(m)*mat1up_neg_cyl*(Rmat_cyl*AmatUD1_cyl')^tmp1(m)*mat1down_neg_cyl*... % first ODE
             Rmat_cyl^dstmp12(m)*mat2up_neg_cyl*(Rmat_cyl*AmatUD2_cyl')^tmp2(m)*mat2down_neg_cyl*...
             Rmat_cyl^dstmp2(m)*mat2up_pos_cyl*(Rmat_cyl*AmatUD2_cyl)^tmp2(m)*mat2down_pos_cyl*real(Smat_cyl); % second ODE
 
         
         Ecal_plane=Ecal_plane* mat1up_pos_plane*(Rmat_plane*AmatUD1_plane)^tmp1(m)*mat1down_pos_plane*...
             Rmat_plane^dstmp1(m)*mat1up_neg_plane*(Rmat_plane*AmatUD1_plane')^tmp1(m)*mat1down_neg_plane*... % first ODE
             Rmat_plane^dstmp12(m)*mat2up_neg_plane*(Rmat_plane*AmatUD2_plane')^tmp2(m)*mat2down_neg_plane*...
             Rmat_plane^dstmp2(m)*mat2up_pos_plane*(Rmat_plane*AmatUD2_plane)^tmp2(m)*mat2down_pos_plane*real(Smat_plane); % second ODE
         else
             error('not implemented yet');
         end
             ePerp(m) = Ecal_cyl;
             ePar(m) = Ecal_plane;
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
        Epert = FiniteCylinder_MM_dPGSE(xpert,protocol);
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
                Epert = FiniteCylinder_MM_dPGSE(xpert,protocol,x_deriv);
                dEtdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEtdx;
            end
        end
        
         
    end   
    protocol.diff=0;
   
end
end
