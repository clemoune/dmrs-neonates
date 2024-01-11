function [E J]=Cylinder_MM_DODE(x,protocol,x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a finite cylinder compartment.
% 
% [E,J]=Cylinder_MM_DODE(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate consisting of parallel, 
% impermeable finite cylinders with a single radius and a diffusion protocol 
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
%       protocol.Nosc - number of periods of the waveform
%       protocol.G1 - gradient strength of the first gradient pair, size [1 N]
%       protocol.G2 - gradient strength of the second gradient pair, size [1 N]
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

dRes = x(1);
theta = x(3);
phi = x(4);

% Find E_restricted
G1=protocol.G1;
G2=protocol.G2;
smalldel1 = protocol.smalldel; % can be extended to different durations
smalldel2 = protocol.smalldel; % can be extended to different durations
tm = protocol.tm; % mixing time (time interval between the pulses)
Nosc = protocol.Nosc;
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


tau = protocol.tau;
ePerp=zeros(M,1); % cylinder restriction

fib_dir = GetFibreOrientation('Cylinder',x);

G1_par = G1.*(grad_dirs1*fib_dir)';
G2_par = G2.*(grad_dirs2*fib_dir)';

% angle_vec = (0:5:175)*pi/180;
GAMMA = 2.675987E8;

if ~isfield(protocol,'slew_rate')


if ~isfield(protocol,'phase') || protocol.phase == 0;
    
    G_ind1 = floor(G1./protocol.gstep);  
    G_ind2 = floor(G2./protocol.gstep); 
    tmp1 = floor((smalldel1./2./Nosc)./tau);
    tmp2 = floor((smalldel2./2./Nosc)./tau);
    dstmp =floor(tm./tau);

    bval_par = GAMMA.^2.*(G1_par.^2./2.*smalldel1.^3./6./Nosc.^2 + G2_par.^2./2.*smalldel2.^3./6./Nosc.^2);
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



             Ecal_cyl=Ecal_cyl*((Rmat_cyl*AmatUD1_cyl)^tmp1(m)*(Rmat_cyl*AmatUD1_cyl')^tmp1(m))^Nosc(m)*... % first ODE
                 Rmat_cyl^dstmp(m)*((Rmat_cyl*AmatUD2_cyl')^tmp2(m)*(Rmat_cyl*AmatUD2_cyl)^tmp2(m))^Nosc(m)*real(Smat_cyl); % second ODE


             ePerp(m) = Ecal_cyl;



    end


elseif protocol.phase == pi/2 || protocol.phse == -pi/2; 

    G_ind1 = floor(G1./protocol.gstep);  
    G_ind2 = floor(G2./protocol.gstep); 

    tmp11 = floor((smalldel1./4./Nosc)./tau);
    tmp12 = floor((smalldel1./2./Nosc)./tau);
    tmp21 = floor((smalldel2./4./Nosc)./tau);
    tmp22 = floor((smalldel2./2./Nosc)./tau);

    Nhp = (Nosc*2); % number of half periods
    dstmp =floor(tm./tau);  

    bval_par = GAMMA.^2.*(G1_par.^2./2.*smalldel1.^3./24./Nosc.^2 + G2_par.^2./2.*smalldel2.^3./24./Nosc.^2);
    logEPar=-bval_par'.*dRes;
    ePar = exp(logEPar);

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


   
         mat11_cyl = (Rmat_cyl*AmatUD1_cyl)^tmp11(m);
         mat11_inv_cyl = (Rmat_cyl*AmatUD1_cyl')^tmp11(m);
         mat21_cyl = (Rmat_cyl*AmatUD2_cyl)^tmp21(m);
         mat21_inv_cyl = (Rmat_cyl*AmatUD2_cyl')^tmp21(m);   
            
             if mod(Nhp(m),2) == 0 % even number of half periods
            Ecal_cyl=Ecal_cyl*(mat11_cyl*mat11_inv_cyl^2*mat11_cyl)^Nosc(m)*... % first ODE
             Rmat_cyl^dstmp(m)*(mat21_inv_cyl*mat21_cyl^2*mat21_inv_cyl)^Nosc(m)*real(Smat_cyl); % second ODE

      
             else % odd number of half periods
              
                 Ecal_cyl=Ecal_cyl*(mat11_cyl*mat11_inv_cyl^2*mat11_cyl)^floor(Nosc(m))*... 
                     mat11_cyl * mat11_inv_cyl *...% first ODE
             Rmat_cyl^dstmp(m)*(mat21_inv_cyl*mat21_cyl^2*mat21_inv_cyl)^floor(Nosc(m))*...
             mat21_inv_cyl*mat21_cyl*real(Smat_cyl); % second ODE

             end
         ePerp(m) = Ecal_cyl;
  
    end
  

end

else % finite slew rate
 slew_rate = protocol.slew_rate;
rt1 = G1./slew_rate;
rt2 = G2./slew_rate;   
    if ~isfield(protocol,'phase') || protocol.phase == 0;
    
    G_ind1 = floor(G1./protocol.gstep);  
    G_ind2 = floor(G2./protocol.gstep); 
    tmp1 = floor((smalldel1./2./Nosc)./tau);
    tmp2 = floor((smalldel2./2./Nosc)./tau);
    dstmp =floor(tm./tau);
    tmp1_rt = floor(rt1./tau);
    tmp2_rt = floor(rt2./tau);
    tmp1 = tmp1 - 2*tmp1_rt;
    tmp2 = tmp2 - 2*tmp2_rt;


    bval_par = GAMMA.^2.*(G1_par.^2.*smalldel1.^3./15./Nosc.^2./4.*(5-15.*rt1.*Nosc./smalldel1 - ...
        5.*rt1.^2.*Nosc.^2./smalldel1+4.*rt1.^3.*8.*Nosc.^3./smalldel1.^3) +...    
         G2_par.^2.*smalldel2.^3./15./Nosc.^2./4.*(5-15.*rt2.*Nosc./smalldel2 - ...
        5.*rt2.^2.*Nosc.^2./smalldel2+4.*rt2.^3.*8.*Nosc.^3./smalldel2.^3));

    logEPar=-bval_par'.*dRes;
    ePar = exp(logEPar);

    for m=1:M
         
       if tmp1(m) < 0 || tmp2(m) < 0
            error(['Frequency of the ' num2str(m) 'th measurement is higher than the one permitted by the slew rate'])
        end
  
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
       
       
       
         Ecal_cyl=Ecal_cyl*( mat1up_pos*(Rmat_cyl*AmatUD1_cyl)^tmp1(m)*mat1down_pos*...
             mat1up_neg*(Rmat_cyl*AmatUD1_cyl')^tmp1(m)*mat1down_neg)^Nosc(m)*... % first ODE
             Rmat_cyl^dstmp(m)*(mat2up_neg*(Rmat_cyl*AmatUD2_cyl')^tmp2(m)*mat2down_neg*...
             mat2up_pos*(Rmat_cyl*AmatUD2_cyl)^tmp2(m)*mat2down_pos)^Nosc(m)*real(Smat_cyl); % second ODE


         ePerp(m) = Ecal_cyl;

     
          
    end


elseif protocol.phase == pi/2 || protocol.phse == -pi/2; 

G_ind1 = floor(G1./protocol.gstep);  
G_ind2 = floor(G2./protocol.gstep); 

tmp11 = floor(((smalldel1./4./Nosc)+rt1/2)./tau);
tmp12 = floor((smalldel1./2./Nosc)./tau);
tmp21 = floor(((smalldel2./4./Nosc)+rt2/2)./tau);
tmp22 = floor((smalldel2./2./Nosc)./tau);

Nhp = (Nosc*2); % number of half periods
dstmp =floor(tm./tau);  
tmp1_rt = floor(rt1./tau);
tmp2_rt = floor(rt2./tau);
tmp11 = tmp11 - 2*tmp1_rt;
tmp12 = tmp12 - 2*tmp1_rt;
tmp21 = tmp21 - 2*tmp2_rt;
tmp22 = tmp22 - 2*tmp2_rt;

bval_par = GAMMA.^2.*(G1_par.^2.*smalldel1.^3./12./4./Nosc.^2.*((1+16.*2.*Nosc).*rt1.^3.*4.*Nosc.^2./5./smalldel1.^3-...
    4.*rt1.^2.*4.*Nosc.^2./smalldel1.^2+1) +...
     G2_par.^2.*smalldel2.^3./12./4./Nosc.^2.*((1+16.*2.*Nosc).*rt2.^3.*4.*Nosc.^2./5./smalldel2.^3-...
    4.*rt2.^2.*4.*Nosc.^2./smalldel2.^2+1));
 
logEPar=-bval_par'.*dRes;
ePar = exp(logEPar);

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
      
            
             if mod(Nhp(m),2) == 0 % even number of half periods
            Ecal_cyl=Ecal_cyl*mat1up_pos*(Rmat_cyl*AmatUD1_cyl)^tmp11(m)*mat1down_pos*...
                (mat1up_neg*(Rmat_cyl*AmatUD1_cyl')^tmp12(m)*mat1down_neg*...
                mat1up_pos*(Rmat_cyl*AmatUD1_cyl)^tmp12(m)*mat1down_pos)^(Nosc(m)-1)*...
                mat1up_neg*(Rmat_cyl*AmatUD1_cyl')^tmp12(m)*mat1down_neg*...
                mat1up_pos*(Rmat_cyl*AmatUD1_cyl)^tmp11(m)*mat1down_pos*...  % first ODE
                Rmat_cyl^dstmp(m)*mat2up_neg*(Rmat_cyl*AmatUD2_cyl')^tmp21(m)*mat2down_neg*...
                (mat2up_pos*(Rmat_cyl*AmatUD2_cyl)^tmp22(m)*mat2down_pos*...
                mat2up_neg*(Rmat_cyl*AmatUD2_cyl')^tmp22(m)*mat2down_neg)^(Nosc(m)-1)...
                *mat2up_pos*(Rmat_cyl*AmatUD2_cyl)^tmp22(m)*mat2down_pos*...
                mat2up_neg*(Rmat_cyl*AmatUD2_cyl')^tmp21(m)*mat2down_neg*real(Smat_cyl); % second ODE

      
             else % odd number of half periods
              
            Ecal_cyl=Ecal_cyl*mat1up_pos*(Rmat_cyl*AmatUD1_cyl)^tmp11(m)*mat1down_pos*...
                (mat1up_neg*(Rmat_cyl*AmatUD1_cyl')^tmp12(m)*mat1down_neg*...
                mat1up_pos*(Rmat_cyl*AmatUD1_cyl)^tmp12(m)*mat1down_pos)^floor(Nosc(m))*...
                mat1up_neg*(Rmat_cyl*AmatUD1_cyl')^tmp11(m)*mat1down_neg*... % first ODE
                Rmat_cyl^dstmp(m)*mat2up_neg*(Rmat_cyl*AmatUD2_cyl')^tmp21(m)*mat2down_neg*...
                (mat2up_pos*(Rmat_cyl*AmatUD2_cyl)^tmp22(m)*mat2down_pos*...
                mat2up_neg*(Rmat_cyl*AmatUD2_cyl')^tmp22(m)*mat2down_neg)^floor(Nosc(m))...
                *mat2up_pos*(Rmat_cyl*AmatUD2_cyl)^tmp21(m)*mat2down_pos*real(Smat_cyl); % second ODE

             end
         ePerp(m) = Ecal_cyl;
  
    end
  

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
        Epert = Cylinder_MM_DODE(xpert,protocol);
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
                Epert = Cylinder_MM_DODE(xpert,protocol,x_deriv);
                dEtdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEtdx;
            end
        end
        
         
    end   
    protocol.diff=0;
   
end
end
