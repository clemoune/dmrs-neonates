function [E,J]=CuboidSym_GPD_PGSE(x, protocol,x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a cuboid/ensemble of cuboids with 2 equal
% sides
% 
% [E,J]=CuboidSym_GPD_PGSE(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate consisting of a cuboid/ 
% an ensemble of cuboids with 2 equal sides and a diffusion protocol 
% specified in the input
% Substrate: cuboid / ensemble of cuboids with 2 equal sides
% Diffusion pulse sequence: Pulsed gradient spin echo (PGSE)
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
%       protocol.grad_dirs - is the gradient direction for each measurement.
%           It has size [N 3] where N is the number of measurements.
%       protocol.G - gradient strength, size [1 N]
%       protocol.delta - pulse separation, size [1 N]
%       protocol.smalldel - pulse duration, size [1 N]
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

roots = (1:10)*2-1;

G = protocol.G';
smalldel = protocol.smalldel';
delta = protocol.delta';
grad_dirs = protocol.grad_dirs;
cube_rotation = protocol.cube_rotation;

% Radial wavenumbers
GAMMA = 2.675987E8; % This is what is used throughout Wuzi.
%GAMMA = 2.6751525E8; % This is the latest best estimate of GAMMA (used in Camino)


l_q=size(grad_dirs,1);
l_a=size(cube_rotation,1);
k_max=numel(roots);

l1_mat=repmat(l1,[l_q l_a k_max]);
l2_mat=repmat(l2,[l_q l_a k_max]);
l3_mat=repmat(l3,[l_q l_a k_max]);


Gx = G.*grad_dirs(:,1);
Gy = G.*grad_dirs(:,2);
Gz = G.*grad_dirs(:,3);

Gxmx = zeros(l_q,l_a);
Gymx = zeros(l_q,l_a);
Gzmx = zeros(l_q,l_a);

for i = 1:l_q
    for j = 1:l_a
        rot_mat = cube_rotation{j};
        vect = rot_mat'*[Gx(i) Gy(i) Gz(i)]';
        Gxmx(i,j) = vect(1);
        Gymx(i,j) = vect(2);
        Gzmx(i,j) = vect(3);
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

deltamx=repmat(delta,[1,l_a k_max]);
smalldelmx=repmat(smalldel,[1,l_a, k_max]);



%Intracell space
exp_smalldel1 = exp(-lambda1_mat.*smalldelmx.*dRes);
exp_delta1 = exp(-lambda1_mat.*deltamx.*dRes);
exp_delta_smalldel1 = exp(-lambda1_mat.*(deltamx-smalldelmx).*dRes);

exp_smalldel2 = exp(-lambda2_mat.*smalldelmx.*dRes);
exp_delta2 = exp(-lambda2_mat.*deltamx.*dRes);
exp_delta_smalldel2 = exp(-lambda2_mat.*(deltamx-smalldelmx).*dRes);

exp_smalldel3 = exp(-lambda3_mat.*smalldelmx.*dRes);
exp_delta3 = exp(-lambda3_mat.*deltamx.*dRes);
exp_delta_smalldel3 = exp(-lambda3_mat.*(deltamx-smalldelmx).*dRes);

% keep the notation from Ianus et.al. JMR2012
Gamma1_mat = lambda1_mat.*dRes.*smalldelmx-1+exp_smalldel1+exp_delta1 -...
    exp_delta_smalldel1/2-exp_smalldel1.*exp_delta1/2;  
Gamma2_mat = lambda2_mat.*dRes.*smalldelmx-1+exp_smalldel2+exp_delta2 -...
    exp_delta_smalldel2/2-exp_smalldel2.*exp_delta2/2; 
Gamma3_mat = lambda3_mat.*dRes.*smalldelmx-1+exp_smalldel3+exp_delta3 -...
    exp_delta_smalldel3/2-exp_smalldel3.*exp_delta3/2; 

sumterms = Gxmx.^2.*sum(B1_mat./lambda1_mat.^2.*Gamma1_mat,3)+...
    Gymx.^2.*sum(B2_mat./lambda2_mat.^2.*Gamma2_mat,3)+...
    Gzmx.^2.*sum(B3_mat./lambda3_mat.^2.*Gamma3_mat,3);
E_intra_mat = exp(-2*GAMMA^2./dRes^2.*sumterms);

E_intra = mean(E_intra_mat,2); 

E = E_intra;
 

 % Compute the Jacobian matrix; computed numerically
if(nargout>1)
    dx = 0.00001;
     J = zeros(length(E), 3);
    if nargin < 3 
        
        for i = 1:3; % compute the derivatives for all model parameters
            xpert = x;
            xpert(i) = xpert(i)*(1+dx);
            Epert = CuboidSym_GPD_PGSE(xpert, protocol);
            dEdx = (Epert - E)/(xpert(i)*dx);
            J(:,i) = dEdx;
        end
    else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
        
        for i = 1:length(x_deriv);
            if x_deriv(i) ~= 0  
                
                xpert = x;
                xpert(i) = xpert(i)*(1+dx);
                Epert = CuboidSym_GPD_PGSE(xpert, protocol);
                dEdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEdx;
            end
        end
       

    end   
 end
