function [E,J]=FiniteCylinder_GPD_TWOGSE(x, protocol, x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a finite cylinder compartment.
% 
% [E,J]=FiniteCylinder_GPD_TPGSE(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate consisting of parallel, 
% impermeable finite cylinders with a single radius and a diffusion protocol 
% specified in the input
% Substrate: parallel, impermeable finite cylinders with a single radius
% Diffusion pulse sequence: Trapezoidal wave oscillating gradient spin echo (TWOGSE)
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
%       protocol.omega - gradient angular frequency, size [1 N]
%       protocol.slew_rate - slew rate of the gradient, size [1 N]
%       protocol.roots_cyl - solutions to the Bessel function equation from 
%       function BesselJ_RootsCyl.m
%       ! TWOGSE sequences require gradient waveforms with integer number of lobes 
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

G = protocol.G';
smalldel = protocol.smalldel';
delta = protocol.delta';
omega = protocol.omega';
grad_dirs = protocol.grad_dirs;
slew_rate = protocol.slew_rate';

% calculate fibre direction from the specified angles
fibredir = [sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];

GAMMA = 2.675987E8;

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

L_mat=repmat(L,[l_q 1]); % cylinder length
L_mat=L_mat(:);
L_mat=repmat(L_mat,[1 1 k_max_plane]);

root_m_cyl=reshape(roots_cyl,[1 1 k_max]);
root_m_plane=reshape(roots_plane,[1 1 k_max_plane]);
alpha_mat_cyl=repmat(root_m_cyl,[l_q*l_a 1 1])./R_mat;
alpha_mat_plane=pi*repmat(root_m_plane,[l_q*l_a 1 1])./L_mat;
amSq_cyl=alpha_mat_cyl.^2;
amSq_plane = alpha_mat_plane.^2;


%Geometric factor B
Bmat_rep_cyl= 2./amSq_cyl./(repmat(root_m_cyl,[l_q*l_a 1 1]).^2 -1);
Bmat_rep_plane = 8./pi^2./amSq_plane./repmat(root_m_plane,[l_q*l_a 1 1]).^2; 

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

for i = 1:length(omega) % safe check 
    if omega(i)<pi/smalldel(i);
        omega(i) = pi/smalldel(i); % smallest omega should correspond for a PGSE sequence with the given smalldel
        if isfield(protocol,'phase')
            protocol.phase(i) = 0;
        end
    end
end
niu = omega/(2*pi());

freqmx=repmat(niu,[1,l_a]);
freqmx_rep = freqmx(:);
freqmx_rep = repmat(freqmx_rep,[1 1 k_max]);
%rise time of the gradient
rt = G./slew_rate;
rtmx = repmat(rt,[1,l_a]);
rtmx_rep = rtmx(:);
rtmx_rep = repmat(rtmx_rep,[1 1 k_max]);

% the number of integer half periods:
NT = floor(2.*smalldel.*niu+0.00000000001); % needed for numerical precision (if NT is int, the floor rounds it down)
NTmx_rep = floor(2.*smalldelmx_rep.*freqmx_rep+0.00000000001);

if max(smalldel - NT./2./niu) > 1E-5
    warning('TWOGSE expressions work only for integer number of lobes - results might be incorrect')
end
if isfield(protocol,'phase')
    if (protocol.phase ~= 0) 
     warning('TWOGSE expressions work only for integer number of lobes - results might be incorrect')
    end
end


              
% Find restricted signal decay
% Parallel component
exp_NT_plane = exp(-amSq_plane.*dRes.*NTmx_rep./(2.*freqmx_rep));
exp_1_plane = exp(-amSq_plane.*dRes./(2.*freqmx_rep));
exp_smalldel_NT_plane = exp(-amSq_plane.*dRes.*(smalldelmx_rep-NTmx_rep./(2.*freqmx_rep)));
exp_delta_NT_plane = exp(-amSq_plane.*dRes.*(deltamx_rep-NTmx_rep./(2.*freqmx_rep)));
exp_delta_plane = exp(-amSq_plane.*dRes.*deltamx_rep);
exp_delta_smalldel_plane = exp(-amSq_plane.*dRes.*(deltamx_rep-smalldelmx_rep));
exp_smalldel_plane = exp(-amSq_plane.*dRes.*smalldelmx_rep);
exp_rt_plane = exp(-amSq_plane.*dRes.*rtmx_rep);
exp_1_rt_plane = exp(-amSq_plane.*dRes.*(1./(2.*freqmx_rep)-rtmx_rep));
exp_1_2rt_plane = exp(-amSq_plane.*dRes.*(1./(2.*freqmx_rep)-2.*rtmx_rep));
am_dRes_plane = amSq_plane.*dRes;
am2_dRes2_plane = amSq_plane.^2.*dRes.^2;
sgn = (-1).^NTmx_rep;

% corresponds to S1 in mathematica (trap_formulae_IntN)

sumterms_plane = -4.*smalldelmx_rep./(3.*am2_dRes2_plane.^2.*rtmx_rep.^2).*...
    (12.*am_dRes_plane.*rtmx_rep.*freqmx_rep - 6.*(1-exp_rt_plane).*freqmx_rep.*...
    (-exp_1_rt_plane+exp_1_2rt_plane+2)+am_dRes_plane.^3.*rtmx_rep.^2.*(8.*freqmx_rep.*rtmx_rep-3));

% corresponds to 2*S2 in Mathematica
sumterms_plane = sumterms_plane - 4.*sgn./am2_dRes2_plane.^2./rtmx_rep.^2./(1+exp_1_plane).^2.*(1-exp_rt_plane).^2.*...
    (exp_1_rt_plane-1).^2.*(exp_smalldel_plane+NTmx_rep.*sgn.*exp_1_plane+(NTmx_rep-1).*sgn);
% corresponds to S3 in Mathematica  
sumterms_plane = sumterms_plane + 2.*exp_delta_smalldel_plane./am2_dRes2_plane.^2./rtmx_rep.^2./(1+exp_1_plane).^2.*...
     (exp_1_rt_plane-1).^2.*(1-exp_rt_plane).^2.*(1-sgn.*exp_smalldel_plane).*(sgn-exp_smalldel_plane);


sumterms_plane = Bmat_rep_plane.*sumterms_plane;

testinds_plane = find(sumterms_plane(:,:,end)>0);
test_cyl = sumterms_plane(testinds_plane,1)./sumterms_plane(testinds_plane,end);
if(min(abs(test_cyl))<1E2)
    warning('Ratio of largest to smallest terms in GPD model sum is <1E2.  May need more terms.');
    x;
end

s_plane = sum(sumterms_plane,3);

s_plane = reshape(s_plane,[l_q,l_a]);
if(min(s_plane)+1E-12<0)
    warning('Negative sums found in GPD sum.  Setting to zero.');
    s_plane(find(s_plane<0))=0;
    x;
end


 logE_plane = -0.5.*GAMMA^2*GmxSq.*cosThetaSq_matrix.*s_plane;
 ePar = exp(logE_plane);



% Perpendicular component

% Find restricted signal decay
% Parallel component
exp_NT_cyl = exp(-amSq_cyl.*dRes.*NTmx_rep./(2.*freqmx_rep));
exp_1_cyl = exp(-amSq_cyl.*dRes./(2.*freqmx_rep));
exp_smalldel_NT_cyl = exp(-amSq_cyl.*dRes.*(smalldelmx_rep-NTmx_rep./(2.*freqmx_rep)));
exp_delta_NT_cyl = exp(-amSq_cyl.*dRes.*(deltamx_rep-NTmx_rep./(2.*freqmx_rep)));
exp_delta_cyl = exp(-amSq_cyl.*dRes.*deltamx_rep);
exp_delta_smalldel_cyl = exp(-amSq_cyl.*dRes.*(deltamx_rep-smalldelmx_rep));
exp_smalldel_cyl = exp(-amSq_cyl.*dRes.*smalldelmx_rep);
exp_rt_cyl = exp(-amSq_cyl.*dRes.*rtmx_rep);
exp_1_rt_cyl = exp(-amSq_cyl.*dRes.*(1./(2.*freqmx_rep)-rtmx_rep));
exp_1_2rt_cyl = exp(-amSq_cyl.*dRes.*(1./(2.*freqmx_rep)-2.*rtmx_rep));
am_dRes_cyl = amSq_cyl.*dRes;
am2_dRes2_cyl = amSq_cyl.^2.*dRes.^2;
sgn = (-1).^NTmx_rep;

% corresponds to S1 in mathematica (trap_formulae_IntN)

sumterms_cyl = -4.*smalldelmx_rep./(3.*am2_dRes2_cyl.^2.*rtmx_rep.^2).*...
    (12.*am_dRes_cyl.*rtmx_rep.*freqmx_rep - 6.*(1-exp_rt_cyl).*freqmx_rep.*...
    (-exp_1_rt_cyl+exp_1_2rt_cyl+2)+am_dRes_cyl.^3.*rtmx_rep.^2.*(8.*freqmx_rep.*rtmx_rep-3));

% corresponds to 2*S2 in Mathematica
sumterms_cyl = sumterms_cyl - 4.*sgn./am2_dRes2_cyl.^2./rtmx_rep.^2./(1+exp_1_cyl).^2.*(1-exp_rt_cyl).^2.*...
    (exp_1_rt_cyl-1).^2.*(exp_smalldel_cyl+NTmx_rep.*sgn.*exp_1_cyl+(NTmx_rep-1).*sgn);
% corresponds to S3 in Mathematica  
sumterms_cyl = sumterms_cyl + 2.*exp_delta_smalldel_cyl./am2_dRes2_cyl.^2./rtmx_rep.^2./(1+exp_1_cyl).^2.*...
     (exp_1_rt_cyl-1).^2.*(1-exp_rt_cyl).^2.*(1-sgn.*exp_smalldel_cyl).*(sgn-exp_smalldel_cyl);

sumterms_cyl = Bmat_rep_cyl.*sumterms_cyl;

testinds_cyl = find(sumterms_cyl(:,:,end)>0);
test_cyl = sumterms_cyl(testinds_cyl,1)./sumterms_cyl(testinds_cyl,end);
if(min(abs(test_cyl))<1E2)
    warning('Ratio of largest to smallest terms in GPD model sum is <1E2.  May need more terms.');
    x;
end

s_cyl = sum(sumterms_cyl,3);
s_cyl = reshape(s_cyl,[l_q,l_a]);
if(min(s_cyl)+1E-12<0)
    warning('Negative sums found in GPD sum.  Setting to zero.');
    s_cyl(find(s_cyl<0))=0;
    x;
end
logE = -0.5.*GAMMA^2*GmxSq.*sinThetaSq_matrix.*s_cyl;
ePerp = exp(logE);
E_r_matrix = ePar.*ePerp;
E_r=sum(E_r_matrix.*weight_matrix,2);
E=E_r;

  

 % Compute the Jacobian matrix; computed numerically
if(nargout>1)
    dx = 0.00001;
    J = zeros(length(E), 5);
    if nargin < 3
         
        for i = 1:5; % compute the derivatives for all model parameters
            xpert = x;
            xpert(i) = xpert(i)*(1+dx);
            Epert = FiniteCylinder_GPD_TWOGSE(xpert, protocol);
            dEdx = (Epert - E)/(xpert(i)*dx);
            J(:,i) = dEdx;
        end

       
    else  % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
      
        for i = 1:length(x_deriv);
            if x_deriv(i) ~= 0  
                
                xpert = x;
                xpert(i) = xpert(i)*(1+dx);
                Epert = FiniteCylinder_GPD_TWOGSE(xpert, protocol);
                dEdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEdx;
            end
        end
    end   
 end
