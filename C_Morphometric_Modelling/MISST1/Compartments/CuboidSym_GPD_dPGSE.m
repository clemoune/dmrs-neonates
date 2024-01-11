function [E,J]=CuboidSym_GPD_dPGSE(x, protocol,x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a cuboid/ensemble of cuboids with 2 equal
% sides
% 
% [E,J]=CuboidSym_GPD_dPGSE(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate consisting of a cuboid/ 
% an ensemble of cuboids with 2 equal sides and a diffusion protocol 
% specified in the input
% Substrate: cuboid / ensemble of cuboids with 2 equal sides
% Diffusion pulse sequence: Double pulsed gradient spin echo (dPGSE)
% Signal approximation: Gaussian Phase Distribution (GPD)  
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% x - size 3 vector of model parameters in SI units for CuboidSym:
%       x(1) - free diffusivity of the material inside the cuboids.
%       x(2) - length of the cuboid (lx = ly)
%       x(3) - eccentricity (ratio between height of the cuboid lz and 
%       length lx)  
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
%       protocol.roots_plane - roots of the diffusion equation for planar
%       geometry
%       protocol.cube_rotation - cell array which stores rotation matrices
%       describing the orientation of the cuboids. If it has more than one
%       matrix, than the average is computed; The default is the identity
%       matrix.
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

dRes=x(1);
l1=x(2);
Ecc=x(3);
l2 =l1;
l3 = l1*Ecc;

roots = protocol.roots_plane;

% get the relevant pusle sequence parameters from protocol
G1 = protocol.G1';
G2 = protocol.G2';
smalldel = protocol.smalldel';
delta = protocol.delta';
tm = protocol.tm';
grad_dirs1 = protocol.grad_dirs1;
grad_dirs2 = protocol.grad_dirs2;
cube_rotation = protocol.cube_rotation;

if min(tm-smalldel) <0
    error('not yet implemented for tm < smalldel')
end

% Radial wavenumbers
GAMMA = 2.675987E8; % This is what is used throughout Wuzi.
%GAMMA = 2.6751525E8; % This is the latest best estimate of GAMMA (used in Camino)


l_q=size(grad_dirs1,1);
l_a=size(cube_rotation,1);
k_max=numel(roots);

l1_mat=repmat(l1,[l_q l_a k_max]);
l2_mat=repmat(l2,[l_q l_a k_max]);
l3_mat=repmat(l3,[l_q l_a k_max]);


G1x = G1.*grad_dirs1(:,1);
G1y = G1.*grad_dirs1(:,2);
G1z = G1.*grad_dirs1(:,3);

G2x = G2.*grad_dirs2(:,1);
G2y = G2.*grad_dirs2(:,2);
G2z = G2.*grad_dirs2(:,3);

G1xmx = zeros(l_q,l_a);
G1ymx = zeros(l_q,l_a);
G1zmx = zeros(l_q,l_a);

G2xmx = zeros(l_q,l_a);
G2ymx = zeros(l_q,l_a);
G2zmx = zeros(l_q,l_a);

for i = 1:l_q
    for j = 1:l_a
        rot_mat = cube_rotation{j};
         vect1 = rot_mat'*[G1x(i) G1y(i) G1z(i)]';
        G1xmx(i,j) = vect1(1);
        G1ymx(i,j) = vect1(2);
        G1zmx(i,j) = vect1(3);
        
        vect2 = rot_mat'*[G2x(i) G2y(i) G2z(i)]';
        G2xmx(i,j) = vect2(1);
        G2ymx(i,j) = vect2(2);
        G2zmx(i,j) = vect2(3);
        
    end
end



root_m=reshape(roots,[1 1 k_max]);
root_mat2 = repmat(root_m,[l_q l_a 1]).^2;
lambda1_mat=pi^2.*root_mat2./l1_mat.^2;
lambda2_mat=pi^2.*root_mat2./l2_mat.^2;
lambda3_mat=pi^2.*root_mat2./l3_mat.^2;

B1_mat = 8*l1_mat.^2./root_mat2.^2./pi^4;
B2_mat = 8*l2_mat.^2./root_mat2.^2./pi^4;
B3_mat = 8*l3_mat.^2./root_mat2.^2./pi^4;

deltamx=repmat(delta,[1,l_a, k_max]);
smalldelmx=repmat(smalldel,[1,l_a, k_max]);
tmmx=repmat(tm,[1,l_a, k_max]);



%Intracell space
exp_smalldel1 = exp(-lambda1_mat.*smalldelmx.*dRes);
exp_delta1 = exp(-lambda1_mat.*deltamx.*dRes);
exp_delta_smalldel1 = exp(-lambda1_mat.*(deltamx-smalldelmx).*dRes);
exp_tm1 = exp(-lambda1_mat*dRes.*tmmx);
exp_tm_smalldel1 = exp(-lambda1_mat.*dRes.*(tmmx-smalldelmx));


exp_smalldel2 = exp(-lambda2_mat.*smalldelmx.*dRes);
exp_delta2 = exp(-lambda2_mat.*deltamx.*dRes);
exp_delta_smalldel2 = exp(-lambda2_mat.*(deltamx-smalldelmx).*dRes);
exp_tm2 = exp(-lambda2_mat*dRes.*tmmx);
exp_tm_smalldel2 = exp(-lambda2_mat.*dRes.*(tmmx-smalldelmx));



exp_smalldel3 = exp(-lambda3_mat.*smalldelmx.*dRes);
exp_delta3 = exp(-lambda3_mat.*deltamx.*dRes);
exp_delta_smalldel3 = exp(-lambda3_mat.*(deltamx-smalldelmx).*dRes);
exp_tm3 = exp(-lambda3_mat*dRes.*tmmx);
exp_tm_smalldel3 = exp(-lambda3_mat.*dRes.*(tmmx-smalldelmx));

am_dRes1 = lambda1_mat.*dRes;
am2_dRes1 = lambda1_mat.^2.*dRes.^2;

am_dRes2 = lambda2_mat.*dRes;
am2_dRes2 = lambda2_mat.^2.*dRes.^2;

am_dRes3 = lambda3_mat.*dRes;
am2_dRes3 = lambda3_mat.^2.*dRes.^2;


G1xmx_mat = repmat(G1xmx,[1 1 k_max]);
G1ymx_mat = repmat(G1ymx,[1 1 k_max]);
G1zmx_mat = repmat(G1zmx,[1 1 k_max]);

G2xmx_mat = repmat(G2xmx,[1 1 k_max]);
G2ymx_mat = repmat(G2ymx,[1 1 k_max]);
G2zmx_mat = repmat(G2zmx,[1 1 k_max]);

% x
  % S11
sumterms_x = 2.*(G1xmx_mat.^2+G2xmx_mat.^2).*(2.*(am_dRes1.*smalldelmx+exp_smalldel1-1)./am2_dRes1);
% S12
sumterms_x = sumterms_x+ 2.*(G1xmx_mat.^2+G2xmx_mat.^2).*(-exp_delta_smalldel1.*(1-exp_smalldel1).^2)./am2_dRes1;
% S13
sumterms_x = sumterms_x- 4.*G1xmx_mat.*G2xmx_mat.*((1-exp_smalldel1).^2.*exp_tm_smalldel1.*exp_delta1)./am2_dRes1;
% S14
sumterms_x = sumterms_x+ 2.*G1xmx_mat.*G2xmx_mat.*((1-exp_smalldel1).^2.*exp_tm1.*exp_delta1.*exp_delta_smalldel1)./am2_dRes1;
% S23
sumterms_x = sumterms_x+ 2*G1xmx_mat.*G2xmx_mat.*((1-exp_smalldel1).^2.*exp_tm_smalldel1)./am2_dRes1;

 sumterms_x = B1_mat.*sumterms_x;

% y
  % S11
sumterms_y = 2.*(G1ymx_mat.^2+G2ymx_mat.^2).*(2.*(am_dRes2.*smalldelmx+exp_smalldel2-1)./am2_dRes2);
% S12
sumterms_y = sumterms_y+ 2.*(G1ymx_mat.^2+G2ymx_mat.^2).*(-exp_delta_smalldel2.*(1-exp_smalldel2).^2)./am2_dRes2;
% S13
sumterms_y = sumterms_y- 4.*G1ymx_mat.*G2ymx_mat.*((1-exp_smalldel2).^2.*exp_tm_smalldel2.*exp_delta2)./am2_dRes2;
% S14
sumterms_y = sumterms_y+ 2.*G1ymx_mat.*G2ymx_mat.*((1-exp_smalldel2).^2.*exp_tm2.*exp_delta2.*exp_delta_smalldel2)./am2_dRes2;
% S23
sumterms_y = sumterms_y+ 2.*G1ymx_mat.*G2ymx_mat.*((1-exp_smalldel2).^2.*exp_tm_smalldel2)./am2_dRes2;

 sumterms_y = B2_mat.*sumterms_y;
 
 
% z
  % S11
sumterms_z = 2.*(G1zmx_mat.^2+G2zmx_mat.^2).*(2.*(am_dRes3.*smalldelmx+exp_smalldel3-1)./am2_dRes3);
% S12
sumterms_z = sumterms_z+ 2.*(G1zmx_mat.^2+G2zmx_mat.^2).*(-exp_delta_smalldel3.*(1-exp_smalldel3).^2)./am2_dRes3;
% S13
sumterms_z = sumterms_z- 4.*G1zmx_mat.*G2zmx_mat.*((1-exp_smalldel3).^2.*exp_tm_smalldel3.*exp_delta3)./am2_dRes3;
% S14
sumterms_z = sumterms_z+ 2.*G1zmx_mat.*G2zmx_mat.*((1-exp_smalldel3).^2.*exp_tm3.*exp_delta3.*exp_delta_smalldel3)./am2_dRes3;
% S23
sumterms_z = sumterms_z+ 2.*G1zmx_mat.*G2zmx_mat.*((1-exp_smalldel3).^2.*exp_tm_smalldel3)./am2_dRes3;

 sumterms_z = B3_mat.*sumterms_z;
 

 
s = sum(sumterms_x,3)+sum(sumterms_y,3)+sum(sumterms_z,3);
if(min(s(:))<0)
    warning('Negative sums found in GPD sum.  Setting to zero.');
    s(find(s<0))=0;
    x;
end

logE = -0.5*GAMMA^2.*s;
E_intra_mat = exp(logE);
E = mean(E_intra_mat,2); 


% 
% % keep the notation from Ianus et.al. JMR2012
% Gamma1_mat = lambda1_mat.*dRes.*smalldelmx-1+exp_smalldel1+exp_delta1 -...
%     exp_delta_smalldel1/2-exp_smalldel1.*exp_delta1/2;  
% Gamma2_mat = lambda2_mat.*dRes.*smalldelmx-1+exp_smalldel2+exp_delta2 -...
%     exp_delta_smalldel2/2-exp_smalldel2.*exp_delta2/2; 
% Gamma3_mat = lambda3_mat.*dRes.*smalldelmx-1+exp_smalldel3+exp_delta3 -...
%     exp_delta_smalldel3/2-exp_smalldel3.*exp_delta3/2; 
% 
% sumterms = Gxmx.^2.*sum(B1_mat./lambda1_mat.^2.*Gamma1_mat,3)+...
%     Gymx.^2.*sum(B2_mat./lambda2_mat.^2.*Gamma2_mat,3)+...
%     Gzmx.^2.*sum(B3_mat./lambda3_mat.^2.*Gamma3_mat,3);
% E_intra_mat = exp(-2*GAMMA^2./dRes^2.*sumterms);
% 
% E_intra = mean(E_intra_mat,2); 
% 
% E = E_intra;
 

 % Compute the Jacobian matrix; computed numerically
if(nargout>1)
    dx = 0.00001;
     J = zeros(length(E), 3);
    if nargin < 3 
        
        for i = 1:3; % compute the derivatives for all model parameters
            xpert = x;
            xpert(i) = xpert(i)*(1+dx);
            Epert = CuboidSym_GPD_dPGSE(xpert, protocol);
            dEdx = (Epert - E)/(xpert(i)*dx);
            J(:,i) = dEdx;
        end
    else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
        
        for i = 1:length(x_deriv);
            if x_deriv(i) ~= 0  
                
                xpert = x;
                xpert(i) = xpert(i)*(1+dx);
                Epert = CuboidSym_GPD_dPGSE(xpert, protocol);
                dEdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEdx;
            end
        end
       

    end   
 end
