function [E,J]=FiniteCylinder_GPD_DSE(x, protocol, x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a finite cylinder compartment.
% 
% [E,J]=FiniteCylinder_GPD_DSE(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate consisting of parallel, 
% impermeable finite cylinders with a single radius and a diffusion protocol 
% specified in the input
% Substrate: parallel, impermeable finite cylinders with a single radius
% Diffusion pulse sequence: dual spin echo (DSE)
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
%       protocol.grad_dirs - is the gradient direction for each measurement. 
%       It has size [N 3] where 
%       protocol.G, protocol.TE, protocol.delta1, protocol.delta2, 
%       protocol.delta3, protocol.t1, protocol.t2, and protocol.t3 are the 
%       gradient strength, echo time, pulse lengths and starting points for
%       each measurement in the protocol.  Each has size [N 1].
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
delta1 = protocol.delta1';
delta2 = protocol.delta2';
delta3 = protocol.delta3';
t1 = protocol.t1';
t2 = protocol.t2';
t3 = protocol.t3';
grad_dirs = protocol.grad_dirs;



% Radial wavenumbers
GAMMA = 2.675987E8; % This is what is used throughout Wuzi.
%GAMMA = 2.6751525E8; % This is the latest best estimate of GAMMA (used in Camino)


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
L_matSq=L_mat.^2;

root_m_cyl=reshape(roots_cyl,[1 1 k_max]);
alpha_mat_cyl=repmat(root_m_cyl,[l_q*l_a 1 1])./R_mat;
amSq_cyl=alpha_mat_cyl.^2;
amP6_cyl=amSq_cyl.^3;

root_m_plane=reshape(roots_plane,[1 1 k_max]);
alpha_mat_plane=pi*repmat(root_m_plane,[l_q*l_a 1 1])./L_mat;
amSq_plane=alpha_mat_plane.^2;
amP6_plane=amSq_plane.^3;


% Angles between gradient directions and fibre direction.
cosTheta = grad_dirs*fibredir;
cosThetaSq = cosTheta.^2;
cosThetaSq_matrix=repmat(cosThetaSq,[1,l_a]);
sinThetaSq_matrix=1-cosThetaSq_matrix;

delta1mx_part=repmat(delta1,[1,l_a]);
delta1mx = delta1mx_part(:);
delta1mx = repmat(delta1mx,[1 1 k_max]);

delta2mx_part=repmat(delta2,[1,l_a]);
delta2mx = delta2mx_part(:);
delta2mx = repmat(delta2mx,[1 1 k_max]);

delta3mx_part=repmat(delta3,[1,l_a]);
delta3mx = delta3mx_part(:);
delta3mx = repmat(delta3mx,[1 1 k_max]);

t1mx_part=repmat(t1,[1,l_a]);
t1mx = t1mx_part(:);
t1mx = repmat(t1mx,[1 1 k_max]);

t2mx_part=repmat(t2,[1,l_a]);
t2mx = t2mx_part(:);
t2mx = repmat(t2mx,[1 1 k_max]);

t3mx_part=repmat(t3,[1,l_a]);
t3mx = t3mx_part(:);
t3mx = repmat(t3mx,[1 1 k_max]);

Gmx=repmat(G,[1,l_a]);
GmxSq = Gmx.^2;


% Perpendicular component (VG model)
vgterm=@(arg) exp(-dRes.*amSq_cyl.*arg);

sumnum_cyl = 2*dRes.*amSq_cyl.*(delta1mx+delta2mx) - 5;
sumnum_cyl = sumnum_cyl - (vgterm(t2mx-t1mx) - vgterm(t3mx-t1mx) - vgterm(t3mx-t2mx) - vgterm(delta1mx)...
    - vgterm(t2mx-t1mx-delta1mx) + vgterm(t3mx-t1mx-delta1mx) - 2*vgterm(delta2mx)...
    - 2*vgterm(t2mx-t1mx+delta2mx) + 2*vgterm(t2mx-t1mx+delta2mx-delta1mx) + ...
    2*vgterm(t3mx-t2mx-delta2mx) - 2*vgterm(delta3mx) + vgterm(delta2mx+delta3mx) + ...
    vgterm(t2mx-t1mx+delta2mx+delta3mx) - vgterm(t2mx-t1mx+delta2mx+delta3mx-delta1mx) - ...
    2*vgterm(t3mx-t2mx+delta1mx-delta3mx) - vgterm(t3mx-t1mx+delta2mx-delta3mx) -...
    vgterm(delta1mx+delta2mx-delta3mx) + vgterm(t3mx-t1mx+delta1mx+delta2mx-delta3mx) +...
    vgterm(t3mx-t2mx+delta1mx+delta2mx-delta3mx) - vgterm(t3mx-t2mx-delta2mx-delta3mx) +...
    vgterm(t3mx-t2mx+delta1mx-2*delta3mx));

sumdenom_cyl = dRes^2*amP6_cyl.*(R_matSq.*amSq_cyl - 1);

% Check for zeros on top and bottom
%sumdenom(find(sumnum) == 0) = 1;
sumterms_cyl = sumnum_cyl./sumdenom_cyl;

testinds = find(sumterms_cyl(:,:,end)>0);
test = sumterms_cyl(testinds,1)./sumterms_cyl(testinds,end);
if(min(test)<1E4)
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


% Parallel component (VG model)
vgterm=@(arg) exp(-dRes.*amSq_plane.*arg);

sumnum_plane = 2*dRes.*amSq_plane.*(delta1mx+delta2mx) - 5;
sumnum_plane = sumnum_plane - (vgterm(t2mx-t1mx) - vgterm(t3mx-t1mx) - vgterm(t3mx-t2mx) - vgterm(delta1mx)...
    - vgterm(t2mx-t1mx-delta1mx) + vgterm(t3mx-t1mx-delta1mx) - 2*vgterm(delta2mx)...
    - 2*vgterm(t2mx-t1mx+delta2mx) + 2*vgterm(t2mx-t1mx+delta2mx-delta1mx) + ...
    2*vgterm(t3mx-t2mx-delta2mx) - 2*vgterm(delta3mx) + vgterm(delta2mx+delta3mx) + ...
    vgterm(t2mx-t1mx+delta2mx+delta3mx) - vgterm(t2mx-t1mx+delta2mx+delta3mx-delta1mx) - ...
    2*vgterm(t3mx-t2mx+delta1mx-delta3mx) - vgterm(t3mx-t1mx+delta2mx-delta3mx) -...
    vgterm(delta1mx+delta2mx-delta3mx) + vgterm(t3mx-t1mx+delta1mx+delta2mx-delta3mx) +...
    vgterm(t3mx-t2mx+delta1mx+delta2mx-delta3mx) - vgterm(t3mx-t2mx-delta2mx-delta3mx) +...
    vgterm(t3mx-t2mx+delta1mx-2*delta3mx));

sumdenom_plane = dRes^2*amP6_plane.*(L_matSq.*amSq_plane)./4;

% Check for zeros on top and bottom
%sumdenom(find(sumnum) == 0) = 1;
sumterms_plane = sumnum_plane./sumdenom_plane;

testinds = find(sumterms_plane(:,:,end)>0);
test = sumterms_plane(testinds,1)./sumterms_plane(testinds,end);
if(min(test)<1E4)
    warning('Ratio of largest to smallest terms in GPD model sum is <1E4.  May need more terms.');
    x;
end

s_plane = sum(sumterms_plane,3);
s_plane = reshape(s_plane,[l_q,l_a]);
if(min(s_plane)<0)
    warning('Negative sums found in GPD sum.  Setting to zero.');
    s_plane(s_plane<0)=0;
    x;
end
%disp(s.*GmxSq)

logE = -2*GAMMA^2*GmxSq.*cosThetaSq_matrix.*s_plane;
ePar = exp(logE);

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
            Epert = FiniteCylinder_GPD_DSE(xpert, protocol);
            dEdx = (Epert - E)/(xpert(i)*dx);
            J(:,i) = dEdx;
        end
    else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
      
        
        for i = 1:length(x_deriv);
            if x_deriv(i) ~= 0  
               
                xpert = x;
                xpert(i) = xpert(i)*(1+dx);
                Epert = FiniteCylinder_GPD_DSE(xpert, protocol);
                dEdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEdx;
            end
        end
        
    end   
end