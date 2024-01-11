function [E,J]=Cylinder_GPD_TWOGSE(x, protocol,x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a cylinder compartment.
% 
% [E,J]=Cylinder_GPD_TPGSE(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate consisting of parallel, 
% impermeable cylinders with a single radius and a diffusion protocol 
% specified in the input
% Substrate: parallel, impermeable cylinders with a single radius
% Diffusion pulse sequence: Trapezoidal wave oscillating gradient spin echo (TWOGSE)
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
k_max=numel(roots);
if isfield(protocol,'Rweight')
    weight_matrix=repmat(protocol.Rweight,[l_q 1]);
else
    weight_matrix = repmat(ones(1,l_a)./l_a,[l_q 1]);
end

R_mat=repmat(R,[l_q 1]);
R_mat=R_mat(:);
R_mat=repmat(R_mat,[1 1 k_max]);


root_m=reshape(roots,[1 1 k_max]);
alpha_mat=repmat(root_m,[l_q*l_a 1 1])./R_mat;
amSq=alpha_mat.^2;


%Geometric factor B
Bmat_rep = 2./amSq./(repmat(root_m,[l_q*l_a 1 1]).^2 -1);

% Angles between gradient directions and fibre direction.
cosTheta = grad_dirs*fibredir;
cosThetaSq = cosTheta.^2;
sinThetaSq = 1-cosThetaSq;
cosThetaSq_matrix=repmat(cosThetaSq,[1,l_a]);
sinThetaSq_matrix=1-cosThetaSq_matrix;

deltamx=repmat(delta,[1,l_a]);
deltamx_rep = deltamx(:);
deltamx_rep = repmat(deltamx_rep,[1 1 k_max]);

smalldelmx=repmat(smalldel,[1,l_a]);
smalldelmx_rep = smalldelmx(:);
smalldelmx_rep = repmat(smalldelmx_rep,[1 1 k_max]);

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

if length(dRes) == 1
    dResmx=repmat(dRes,[l_q,l_a]);
    dResmx_rep = dResmx(:);
    dResmx_rep = repmat(dResmx_rep,[1 1 k_max]);
elseif length(dRes) == length(smalldel)
    dResmx=repmat(dRes,[1,l_a]);
    dResmx_rep = dResmx(:);
    dResmx_rep = repmat(dResmx_rep,[1 1 k_max]);
else 
    error('D must be specified for each measurements')
end

if max(smalldel - NT./2./niu) > 1E-5
    warning('TWOGSE expressions work only for integer number of lobes - results might be incorrect')
end
if isfield(protocol,'phase')
    if (protocol.phase ~= 0) 
     warning('TWOGSE expressions work only for integer number of lobes - results might be incorrect')
    end
end


Gmx=repmat(G,[1,l_a]);
GmxSq = Gmx.^2;

% Find hindered signals

bval = GAMMA.^2.*G.^2.*((delta-smalldel)./4.*(1-(-1).^NT).^2.*(1./2./niu-rt).^2+...
    smalldel./(240.*niu.^2).*(40-120.*rt.*niu-40.*rt.^2.*niu.^2+256.*rt.^3.*niu.^3));

% Find restricted signal decay
bval_matrix=repmat(bval,[1,l_a]);

% Parallel component
logEPar=-bval_matrix.*cosThetaSq_matrix.*dResmx;
ePar=exp(logEPar);

% Perpendicular component

%precomputing values that appear often in the expresions:
% rearrange all formulas to have only negative exponentials
exp_NT = exp(-amSq.*dResmx_rep.*NTmx_rep./(2.*freqmx_rep));
exp_1 = exp(-amSq.*dResmx_rep./(2.*freqmx_rep));
exp_smalldel_NT = exp(-amSq.*dResmx_rep.*(smalldelmx_rep-NTmx_rep./(2.*freqmx_rep)));
exp_delta_NT = exp(-amSq.*dResmx_rep.*(deltamx_rep-NTmx_rep./(2.*freqmx_rep)));
exp_delta = exp(-amSq.*dResmx_rep.*deltamx_rep);
exp_delta_smalldel = exp(-amSq.*dResmx_rep.*(deltamx_rep-smalldelmx_rep));
exp_smalldel = exp(-amSq.*dResmx_rep.*smalldelmx_rep);
exp_rt = exp(-amSq.*dResmx_rep.*rtmx_rep);
exp_1_rt = exp(-amSq.*dResmx_rep.*(1./(2.*freqmx_rep)-rtmx_rep));
exp_1_2rt = exp(-amSq.*dResmx_rep.*(1./(2.*freqmx_rep)-2.*rtmx_rep));
am_dRes = amSq.*dResmx_rep;
am2_dRes2 = amSq.^2.*dResmx_rep.^2;
sgn = (-1).^NTmx_rep;

% corresponds to S1 in mathematica (trap_formulae_IntN)

sumterms = -4.*smalldelmx_rep./(3.*am2_dRes2.^2.*rtmx_rep.^2).*...
    (12.*am_dRes.*rtmx_rep.*freqmx_rep - 6.*(1-exp_rt).*freqmx_rep.*...
    (-exp_1_rt+exp_1_2rt+2)+am_dRes.^3.*rtmx_rep.^2.*(8.*freqmx_rep.*rtmx_rep-3));

% corresponds to 2*S2 in Mathematica
sumterms = sumterms - 4.*sgn./am2_dRes2.^2./rtmx_rep.^2./(1+exp_1).^2.*(1-exp_rt).^2.*...
    (exp_1_rt-1).^2.*(exp_smalldel+NTmx_rep.*sgn.*exp_1+(NTmx_rep-1).*sgn);
% corresponds to S3 in Mathematica  
sumterms = sumterms + 2.*exp_delta_smalldel./am2_dRes2.^2./rtmx_rep.^2./(1+exp_1).^2.*...
     (exp_1_rt-1).^2.*(1-exp_rt).^2.*(1-sgn.*exp_smalldel).*(sgn-exp_smalldel);


sumterms = Bmat_rep.*sumterms;

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

% The original version had (1-cosThetaSq_matrix) squared, but I don't know
% why...  Removed for now.

logE = -0.5.*GAMMA^2*GmxSq.*sinThetaSq_matrix.*s;
ePerp = exp(logE);
ePerp(G==0)=1;

E_r_matrix = ePar.*ePerp; 

E=sum(E_r_matrix.*weight_matrix,2);


 % Compute the Jacobian matrix; computed numerically
if(nargout>1)
    dx = 0.0001;
     J = zeros(length(E), 4);
    if nargin < 3 
        
        for i = 1:4; % compute the derivatives for all model parameters
            xpert = x;
            xpert(i) = xpert(i)*(1+dx);
            Epert = Cylinder_GPD_TWOGSE(xpert, protocol);
            dEdx = (Epert - E)/(xpert(i)*dx);
            J(:,i) = dEdx;
        end
    else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
       
        for i = 1:length(x_deriv);
            if x_deriv(i) ~= 0  
                xpert = x;
                xpert(i) = xpert(i)*(1+dx);
                Epert = Cylinder_GPD_TWOGSE(xpert, protocol);
                dEdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEdx;
            end
        end
       

    end   
 end