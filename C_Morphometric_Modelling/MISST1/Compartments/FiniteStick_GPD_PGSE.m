function [E,J]=FiniteStick_GPD_PGSE(x, protocol, x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a finite cylinder compartment.
% 
% [E,J]=FiniteStick_GPD_PGSE(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate consisting of parallel, 
% impermeable finite sticks with a single length and a diffusion protocol 
% specified in the input
% Substrate: parallel, impermeable finite sticks with a single length
% Diffusion pulse sequence: Pulsed gradient spin echo (PGSE) 
% Signal approximation: Gaussian Phase Distribution (GPD)  
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% x - size 5 vector of model parameters in SI units for FiniteCylinder:
%       x(1) - free diffusivity of the material inside the cylinders.
%       x(2) - length of stick
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

roots_plane = protocol.roots_plane;


if iscell(x)
    dRes=cell2mat(x(1));
    L=cell2mat(x(2)); 
    theta = cell2mat(x(3));
    phi = cell2mat(x(4));
else
    dRes=x(1);
    L = x(2);
    theta = x(3);
    phi = x(4);
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
l_a=numel(L);
k_max_plane = numel(roots_plane);

if isfield(protocol,'Rweight')
    weight_matrix=repmat(protocol.Rweight,[l_q 1]);
else
    weight_matrix = repmat(ones(1,l_a)./l_a,[l_q 1]);
end


L_mat=repmat(L,[l_q 1]); % cylinder length
L_mat=L_mat(:);
L_mat=repmat(L_mat,[1 1 k_max_plane]);



root_m_plane=reshape(roots_plane,[1 1 k_max_plane]);
root_mat2_plane = repmat(root_m_plane,[l_q*l_a 1 1]).^2;
lambda_mat_plane=pi^2.*root_mat2_plane./L_mat.^2;
B_mat_plane = 8*L_mat.^2./root_mat2_plane.^2./pi^4;

deltamx=repmat(delta,[1,l_a]);
deltamx_rep = deltamx(:);
deltamx_rep = repmat(deltamx_rep,[1 1 k_max_plane]);

smalldelmx=repmat(smalldel,[1,l_a]);
smalldelmx_rep = smalldelmx(:);
smalldelmx_rep = repmat(smalldelmx_rep,[1 1 k_max_plane]);


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




E_r_matrix = ePar;

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
            Epert = FiniteStick_GPD_PGSE(xpert, protocol);
            dEdx = (Epert - E)/(xpert(i)*dx);
            J(:,i) = dEdx;
        end
    else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
      
        
        for i = 1:length(x_deriv);
            if x_deriv(i) ~= 0  
               
                xpert = x;
                xpert(i) = xpert(i)*(1+dx);
                Epert = FiniteStick_GPD_PGSE(xpert, protocol);
                dEdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEdx;
            end
        end
        
    end   
end