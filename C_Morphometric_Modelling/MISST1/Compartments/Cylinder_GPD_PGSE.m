function [E,J]=Cylinder_GPD_PGSE(x, protocol,x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a cylinder compartment.
% 
% [E,J]=Cylinder_GPD_PGSE(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate consisting of parallel, 
% impermeable cylinders with a single radius and a diffusion protocol 
% specified in the input
% Substrate: parallel, impermeable cylinders with a single radius
% Diffusion pulse sequence: Pulsed gradient spin echo (PGSE) 
% Signal approximation: Gaussian Phase Distribution (GPD)  
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

roots = protocol.roots_cyl;

% Check the roots array is correct
if(abs(roots(1) - 1.8412)>0.0001)
    error('Looks like the roots array is wrong.  First value should be 1.8412, but is %f', roots(1));
end

if iscell(x)
    dRes=cell2mat(x(1));
    R=cell2mat(x(2)); 
    theta = cell2mat(x(3));
    phi = cell2mat(x(4));
else
    dRes=x(1);
    R=x(2); 
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
l_a=numel(R);
k_max=numel(roots);

if isfield(protocol,'Rweight')
    weight_matrix=repmat(protocol.Rweight,[l_q 1]);
else
    weight_matrix = repmat(ones(1,l_a)./l_a,[l_q 1]);
end

R_mat=repmat(R,[l_q 1]);
R_mat=R_mat(:);
R_mat=repmat(R_mat,[1 1 k_max]);
R_matSq=R_mat.^2;

root_m=reshape(roots,[1 1 k_max]);
alpha_mat=repmat(root_m,[l_q*l_a 1 1])./R_mat;
amSq=alpha_mat.^2;
amP6=amSq.^3;

% Angles between gradient directions and fibre direction.
cosTheta = grad_dirs*fibredir;
cosThetaSq = cosTheta.^2;
cosThetaSq_matrix=repmat(cosThetaSq,[1,l_a]);
sinThetaSq_matrix=1-cosThetaSq_matrix;

deltamx=repmat(delta,[1,l_a]);
deltamx_rep = deltamx(:);
deltamx_rep = repmat(deltamx_rep,[1 1 k_max]);

smalldelmx=repmat(smalldel,[1,l_a]);
smalldelmx_rep = smalldelmx(:);
smalldelmx_rep = repmat(smalldelmx_rep,[1 1 k_max]);

Gmx=repmat(G,[1,l_a]);
GmxSq = Gmx.^2;

% Need the diffusion times in a matrix for the computation below.
difftime = (delta-smalldel/3);
difftime_matrix = repmat(difftime,[1,l_a]);


% Find restricted signal decay
modQSq_matrix=repmat(modQ_Sq,[1,l_a]);

% Parallel component
logEPar=-modQSq_matrix.*cosThetaSq_matrix.*difftime_matrix*dRes;
ePar=exp(logEPar);

% Perpendicular component (VG model)
sda2 = smalldelmx_rep.*amSq;
bda2 = deltamx_rep.*amSq;
emdsda2 = exp(-dRes*sda2);
emdbda2 = exp(-dRes*bda2);
emdbdmsda2 = exp(-dRes*(bda2 - sda2));
emdbdpsda2 = exp(-dRes*(bda2 + sda2));

sumnum = 2*dRes*sda2 - 2;
sumnum = sumnum + 2*emdsda2 + 2*emdbda2;
sumnum = sumnum - emdbdmsda2 - emdbdpsda2;

sumdenom = dRes^2*amP6.*(R_matSq.*amSq - 1);

% Check for zeros on top and bottom
%sumdenom(find(sumnum) == 0) = 1;
sumterms = sumnum./sumdenom;

testinds = find(sumterms(:,:,end)>0);
test = sumterms(testinds,1)./sumterms(testinds,end);
if(min(test)<1E4)
    warning('Ratio of largest to smallest terms in GPD model sum is <1E4.  May need more terms.');
    x;
end

s = sum(sumterms,3);
s = reshape(s,[l_q,l_a]);
if(min(s)<0)
    warning('Negative sums found in GPD sum.  Setting to zero.');
    s(find(s<0))=0;
    x;
end
%disp(s.*GmxSq)

logE = -2*GAMMA^2*GmxSq.*sinThetaSq_matrix.*s;
ePerp = exp(logE);

E_r_matrix = ePar.*ePerp;

E_r=sum(E_r_matrix.*weight_matrix,2);

E=E_r;


% Compute the Jacobian matrix; computed numerically
if(nargout>1)
    dx = 0.00001;
     J = zeros(length(E), length(x));
    if nargin < 3 
        
        for i = 1:4; % compute the derivatives for all model parameters
            xpert = x;
            xpert(i) = xpert(i)*(1+dx);
            Epert = Cylinder_GPD_PGSE(xpert, protocol);
            dEdx = (Epert - E)/(xpert(i)*dx);
            J(:,i) = dEdx;
        end
        
       
    else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
    
        for i = 1:length(x_deriv);
            if x_deriv(i) ~= 0                  
                xpert = x;
                xpert(i) = xpert(i)*(1+dx);
                Epert = Cylinder_GPD_PGSE(xpert, protocol);
                dEdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEdx;
            end
        end
        
       
    end   
end