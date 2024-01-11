function [E,J]=FiniteCylinder_GPD_PGSE(x, protocol, x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a finite cylinder compartment.
% 
% [E,J]=FiniteCylinder_GPD_PGSE(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate consisting of parallel, 
% impermeable finite cylinders with a single radius and a diffusion protocol 
% specified in the input
% Substrate: parallel, impermeable finite cylinders with a single radius
% Diffusion pulse sequence: Pulsed gradient spin echo (PGSE) 
% Signal approximation: Gaussian Phase Distribution (GPD)  
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
%       protocol.grad_dirs - is the gradient direction of the gradient for
%       each measurement. It has size [N 3] where N is the number of measurements.
%       protocol.G - gradient strength, size [1 N]
%       protocol.delta - pulse separation, size [1 N]
%       protocol.smalldel - pulse duration, size [1 N]
%       protocol.roots_cyl - solutions to the Bessel function equation from 
%       function BesselJ_RootsCyl.m
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

roots_cyl = protocol.roots_cyl;
roots_plane = protocol.roots_plane;

% Check the roots array is correct
if(abs(roots_cyl(1) - 1.8412)>0.0001)
    error('Looks like the roots array is wrong.  First value should be 1.8412, but is %f', roots_cyl(1));
end

if iscell(x)
    dRes=cell2mat(x(1));
    R=cell2mat(x(2)); 
    ecc = cell2mat(x(3));
    theta = cell2mat(x(4));
    phi = cell2mat(x(5));
else
    dRes=x(1);
    R=x(2); 
    ecc = x(3);
    theta = x(4);
    phi = x(5);
end

L = 2*R*ecc;

if length(roots_cyl) ~=length(roots_plane)
    error('the root vectors should have the same length');
end

% calculate fibre direction from the specified angles
fibredir = [sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];

% get the relevalt pusle sequence parameters from protocol
G = protocol.G';
smalldel = protocol.smalldel';
delta = protocol.delta';
grad_dirs = protocol.grad_dirs;


% Radial wavenumbers
GAMMA = 2.675987E8; % This is what is used throughout Wuzi.
%GAMMA = 2.6751525E8; % This is the latest best estimate of GAMMA (used in Camino)
modQ = GAMMA*smalldel.*G;
modQ_Sq = modQ.^2;

l_q=size(grad_dirs,1);
l_a=numel(R);
k_max=numel(roots_cyl);
k_max_plane = numel(roots_plane);

if isfield(protocol,'Rweight')
    weight_matrix=repmat(protocol.Rweight,[l_q 1]);
else
    weight_matrix = repmat(ones(1,l_a)./l_a,[l_q 1]);
end

R_mat=repmat(R,[l_q 1]); % cylinder radius
R_mat=R_mat(:);
R_mat=repmat(R_mat,[1 1 k_max]);
R_matSq=R_mat.^2;

L_mat=repmat(L,[l_q 1]); % cylinder length
L_mat=L_mat(:);
L_mat=repmat(L_mat,[1 1 k_max_plane]);

root_m_cyl=reshape(roots_cyl,[1 1 k_max]);
alpha_mat_cyl=repmat(root_m_cyl,[l_q*l_a 1 1])./R_mat;
amSq_cyl=alpha_mat_cyl.^2;
amP6_cyl=amSq_cyl.^3;

root_m_plane=reshape(roots_plane,[1 1 k_max_plane]);
root_mat2_plane = repmat(root_m_plane,[l_q*l_a 1 1]).^2;
lambda_mat_plane=pi^2.*root_mat2_plane./L_mat.^2;
B_mat_plane = 8*L_mat.^2./root_mat2_plane.^2./pi^4;

deltamx=repmat(delta,[1,l_a]);
deltamx_rep = deltamx(:);
deltamx_rep = repmat(deltamx_rep,[1 1 k_max]);

smalldelmx=repmat(smalldel,[1,l_a]);
smalldelmx_rep = smalldelmx(:);
smalldelmx_rep = repmat(smalldelmx_rep,[1 1 k_max]);


% Angles between gradient directions and fibre direction.
cosTheta = grad_dirs*fibredir;
cosThetaSq = cosTheta.^2;
cosThetaSq_matrix=repmat(cosThetaSq,[1,l_a]);
sinThetaSq_matrix=1-cosThetaSq_matrix;

Gmx=repmat(G,[1,l_a]);
GmxSq = Gmx.^2;

% Parallel component - planar restriction
exp_smalldel_plane = exp(-lambda_mat_plane.*smalldelmx_rep.*dRes);
exp_delta_plane = exp(-lambda_mat_plane.*deltamx_rep.*dRes);
exp_delta_smalldel_plane = exp(-lambda_mat_plane.*(deltamx_rep-smalldelmx_rep).*dRes);

% keep the notation from Ianus et.al. JMR2012
Gamma_mat_plane = lambda_mat_plane.*dRes.*smalldelmx_rep-1+exp_smalldel_plane+exp_delta_plane -...
    exp_delta_smalldel_plane/2-exp_smalldel_plane.*exp_delta_plane/2;  

s_plane = sum(B_mat_plane./lambda_mat_plane.^2.*Gamma_mat_plane,3);
s_plane = reshape(s_plane,[l_q,l_a]);

 logE_plane = -2.*GAMMA^2*GmxSq.*cosThetaSq_matrix./dRes^2.*s_plane;
 ePar = exp(logE_plane);


% Perpendicular component (VG model)
sda2_cyl = smalldelmx_rep.*amSq_cyl;
bda2_cyl = deltamx_rep.*amSq_cyl;
emdsda2_cyl = exp(-dRes*sda2_cyl);
emdbda2_cyl = exp(-dRes*bda2_cyl);
emdbdmsda2_cyl = exp(-dRes*(bda2_cyl - sda2_cyl));
emdbdpsda2_cyl = exp(-dRes*(bda2_cyl + sda2_cyl));

sumnum_cyl = 2*dRes*sda2_cyl - 2;
sumnum_cyl = sumnum_cyl + 2*emdsda2_cyl + 2*emdbda2_cyl;
sumnum_cyl = sumnum_cyl - emdbdmsda2_cyl - emdbdpsda2_cyl;

sumdenom_cyl = dRes^2*amP6_cyl.*(R_matSq.*amSq_cyl - 1);

% Check for zeros on top and bottom
%sumdenom(find(sumnum) == 0) = 1;
sumterms_cyl = sumnum_cyl./sumdenom_cyl;

testinds_cyl = find(sumterms_cyl(:,:,end)>0);
test_cyl = sumterms_cyl(testinds_cyl,1)./sumterms_cyl(testinds_cyl,end);
if(min(test_cyl)<1E4)
    warning('Ratio of largest to smallest terms in GPD model sum is <1E4.  May need more terms.');
    x;
end

s_cyl = sum(sumterms_cyl,3);
s_cyl = reshape(s_cyl,[l_q,l_a]);
if(min(s_cyl)<0)
    warning('Negative sums found in GPD sum.  Setting to zero.');
    s_cyl(find(s_cyl<0))=0;
    x;
end
%disp(s.*GmxSq)

logE = -2*GAMMA^2*GmxSq.*sinThetaSq_matrix.*s_cyl;
ePerp = exp(logE);

E_r_matrix = ePar.*ePerp;

E_r=sum(E_r_matrix.*weight_matrix,2);

E=E_r;


% Compute the Jacobian matrix; computed numerically
if(nargout>1)
    dx = 0.00001;
      J = zeros(length(E), length(x));
    if nargin< 3 
       
        for i = 1:length(x); % compute the derivatives for all model parameters
            xpert = x;
            xpert(i) = xpert(i)*(1+dx);
            Epert = FiniteCylinder_GPD_PGSE(xpert, protocol);
            dEdx = (Epert - E)/(xpert(i)*dx);
            J(:,i) = dEdx;
        end
    else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
      
        
        for i = 1:length(x_deriv);
            if x_deriv(i) ~= 0  
               
                xpert = x;
                xpert(i) = xpert(i)*(1+dx);
                Epert = FiniteCylinder_GPD_PGSE(xpert, protocol);
                dEdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEdx;
            end
        end
        
    end   
end