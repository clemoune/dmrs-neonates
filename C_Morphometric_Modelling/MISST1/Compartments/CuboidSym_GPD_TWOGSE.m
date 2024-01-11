function [E,J]=CuboidSym_GPD_TWOGSE(x, protocol, x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a cuboid/ensemble of cuboids with 2 equal
% sides
% 
% [E,J]=CuboidSym_GPD_TWOGSE(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate consisting of a cuboid/ 
% an ensemble of cuboids with 2 equal sides and a diffusion protocol 
% specified in the input
% Substrate: cuboid / ensemble of cuboids with 2 equal sides
% Diffusion pulse sequence: Trapezoidal wave oscillating gradient spin echo (TWOGSE)
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
%       protocol.omega - gradient angular frequency, size [1 N]
%       protocol.slew_rate - slew rate of the gradient, size [1 N]
%       protocol.roots_plane - roots of the diffusion equation for planar
%       geometry
%       protocol.cube_rotation - cell array which stores rotation matrices
%       describing the orientation of the cuboids. If it has more than one
%       matrix, than the average is computed; The default is the identity
%       matrix.
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

dRes=x(1);
l1=x(2);
Ecc=x(3);

l2 = l1;
l3 = Ecc*l1;

roots = (1:10)*2-1;

G = protocol.G';
smalldel = protocol.smalldel';
delta = protocol.delta';
omega = protocol.omega';
grad_dirs = protocol.grad_dirs;
cube_rotation = protocol.cube_rotation;
slew_rate = protocol.slew_rate';

GAMMA = 2.675987E8;

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

deltamx=repmat(delta,[1,l_a]);
deltamx_rep = deltamx(:);
deltamx_rep = repmat(deltamx_rep,[1 1 k_max]);

smalldelmx=repmat(smalldel,[1,l_a]);
smalldelmx_rep = smalldelmx(:);
smalldelmx_rep = repmat(smalldelmx_rep,[1 1 k_max]);

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

exp_NT_x = exp(-lambda1_mat.*dRes.*NTmx_rep./(2.*freqmx_rep));
exp_1_x = exp(-lambda1_mat.*dRes./(2.*freqmx_rep));
exp_smalldel_NT_x = exp(-lambda1_mat.*dRes.*(smalldelmx_rep-NTmx_rep./(2.*freqmx_rep)));
exp_delta_NT_x = exp(-lambda1_mat.*dRes.*(deltamx_rep-NTmx_rep./(2.*freqmx_rep)));
exp_delta_x = exp(-lambda1_mat.*dRes.*deltamx_rep);
exp_delta_smalldel_x = exp(-lambda1_mat.*dRes.*(deltamx_rep-smalldelmx_rep));
exp_smalldel_x = exp(-lambda1_mat.*dRes.*smalldelmx_rep);
exp_rt_x = exp(-lambda1_mat.*dRes.*rtmx_rep);
exp_1_rt_x = exp(-lambda1_mat.*dRes.*(1./(2.*freqmx_rep)-rtmx_rep));
exp_1_2rt_x = exp(-lambda1_mat.*dRes.*(1./(2.*freqmx_rep)-2.*rtmx_rep));
am_dRes_x = lambda1_mat.*dRes;
am2_dRes2_x = lambda1_mat.^2.*dRes.^2;
sgn = (-1).^NTmx_rep;


exp_NT_y = exp(-lambda2_mat.*dRes.*NTmx_rep./(2.*freqmx_rep));
exp_1_y = exp(-lambda2_mat.*dRes./(2.*freqmx_rep));
exp_smalldel_NT_y = exp(-lambda2_mat.*dRes.*(smalldelmx_rep-NTmx_rep./(2.*freqmx_rep)));
exp_delta_NT_y = exp(-lambda2_mat.*dRes.*(deltamx_rep-NTmx_rep./(2.*freqmx_rep)));
exp_delta_y = exp(-lambda2_mat.*dRes.*deltamx_rep);
exp_delta_smalldel_y = exp(-lambda2_mat.*dRes.*(deltamx_rep-smalldelmx_rep));
exp_smalldel_y = exp(-lambda2_mat.*dRes.*smalldelmx_rep);
exp_rt_y = exp(-lambda2_mat.*dRes.*rtmx_rep);
exp_1_rt_y = exp(-lambda2_mat.*dRes.*(1./(2.*freqmx_rep)-rtmx_rep));
exp_1_2rt_y = exp(-lambda2_mat.*dRes.*(1./(2.*freqmx_rep)-2.*rtmx_rep));
am_dRes_y = lambda2_mat.*dRes;
am2_dRes2_y = lambda2_mat.^2.*dRes.^2;


exp_NT_z = exp(-lambda3_mat.*dRes.*NTmx_rep./(2.*freqmx_rep));
exp_1_z = exp(-lambda3_mat.*dRes./(2.*freqmx_rep));
exp_smalldel_NT_z = exp(-lambda3_mat.*dRes.*(smalldelmx_rep-NTmx_rep./(2.*freqmx_rep)));
exp_delta_NT_z = exp(-lambda3_mat.*dRes.*(deltamx_rep-NTmx_rep./(2.*freqmx_rep)));
exp_delta_z = exp(-lambda3_mat.*dRes.*deltamx_rep);
exp_delta_smalldel_z = exp(-lambda3_mat.*dRes.*(deltamx_rep-smalldelmx_rep));
exp_smalldel_z = exp(-lambda3_mat.*dRes.*smalldelmx_rep);
exp_rt_z = exp(-lambda3_mat.*dRes.*rtmx_rep);
exp_1_rt_z = exp(-lambda3_mat.*dRes.*(1./(2.*freqmx_rep)-rtmx_rep));
exp_1_2rt_z = exp(-lambda3_mat.*dRes.*(1./(2.*freqmx_rep)-2.*rtmx_rep));
am_dRes_z = lambda3_mat.*dRes;
am2_dRes2_z = lambda3_mat.^2.*dRes.^2;

 
% x
% corresponds to S1 in mathematica (trap_formulae_IntN)
sumterms_x = -4.*smalldelmx_rep./(3.*am2_dRes2_x.^2.*rtmx_rep.^2).*...
    (12.*am_dRes_x.*rtmx_rep.*freqmx_rep - 6.*(1-exp_rt_x).*freqmx_rep.*...
    (-exp_1_rt_x+exp_1_2rt_x+2)+am_dRes_x.^3.*rtmx_rep.^2.*(8.*freqmx_rep.*rtmx_rep-3));
% corresponds to 2*S2 in Mathematica
sumterms_x = sumterms_x - 4.*sgn./am2_dRes2_x.^2./rtmx_rep.^2./(1+exp_1_x).^2.*(1-exp_rt_x).^2.*...
    (exp_1_rt_x-1).^2.*(exp_smalldel_x+NTmx_rep.*sgn.*exp_1_x+(NTmx_rep-1).*sgn);
% corresponds to S3 in Mathematica  
sumterms_x = sumterms_x + 2.*exp_delta_smalldel_x./am2_dRes2_x.^2./rtmx_rep.^2./(1+exp_1_x).^2.*...
     (exp_1_rt_x-1).^2.*(1-exp_rt_x).^2.*(1-sgn.*exp_smalldel_x).*(sgn-exp_smalldel_x);

% y 
% corresponds to S1 in mathematica (trap_formulae_IntN)
sumterms_y = -4.*smalldelmx_rep./(3.*am2_dRes2_y.^2.*rtmx_rep.^2).*...
    (12.*am_dRes_y.*rtmx_rep.*freqmx_rep - 6.*(1-exp_rt_y).*freqmx_rep.*...
    (-exp_1_rt_y+exp_1_2rt_y+2)+am_dRes_y.^3.*rtmx_rep.^2.*(8.*freqmx_rep.*rtmx_rep-3));
% corresponds to 2*S2 in Mathematica
sumterms_y = sumterms_y - 4.*sgn./am2_dRes2_y.^2./rtmx_rep.^2./(1+exp_1_y).^2.*(1-exp_rt_y).^2.*...
    (exp_1_rt_y-1).^2.*(exp_smalldel_y+NTmx_rep.*sgn.*exp_1_y+(NTmx_rep-1).*sgn);
% corresponds to S3 in Mathematica  
sumterms_y = sumterms_y + 2.*exp_delta_smalldel_y./am2_dRes2_y.^2./rtmx_rep.^2./(1+exp_1_y).^2.*...
     (exp_1_rt_y-1).^2.*(1-exp_rt_y).^2.*(1-sgn.*exp_smalldel_y).*(sgn-exp_smalldel_y); 
 
 % z
 % corresponds to S1 in mathematica (trap_formulae_IntN)
sumterms_z = -4.*smalldelmx_rep./(3.*am2_dRes2_z.^2.*rtmx_rep.^2).*...
    (12.*am_dRes_z.*rtmx_rep.*freqmx_rep - 6.*(1-exp_rt_z).*freqmx_rep.*...
    (-exp_1_rt_z+exp_1_2rt_z+2)+am_dRes_z.^3.*rtmx_rep.^2.*(8.*freqmx_rep.*rtmx_rep-3));
% corresponds to 2*S2 in Mathematica
sumterms_z = sumterms_z - 4.*sgn./am2_dRes2_z.^2./rtmx_rep.^2./(1+exp_1_z).^2.*(1-exp_rt_z).^2.*...
    (exp_1_rt_z-1).^2.*(exp_smalldel_z+NTmx_rep.*sgn.*exp_1_z+(NTmx_rep-1).*sgn);
% corresponds to S3 in Mathematica  
sumterms_z = sumterms_z + 2.*exp_delta_smalldel_z./am2_dRes2_z.^2./rtmx_rep.^2./(1+exp_1_z).^2.*...
     (exp_1_rt_z-1).^2.*(1-exp_rt_z).^2.*(1-sgn.*exp_smalldel_z).*(sgn-exp_smalldel_z);

sumterms_x = B1_mat.*sumterms_x;
sumterms_y = B2_mat.*sumterms_y;
sumterms_z = B3_mat.*sumterms_z;

testinds = find(sumterms_x(:,:,end)>0);
test = sumterms_x(testinds,1)./sumterms_x(testinds,end);
if(min(abs(test))<1E2)
    warning('Ratio of largest to smallest terms in GPD model sum is <1E2.  May need more terms.');
    x;
end

s_x = sum(sumterms_x,3);
s_x = reshape(s_x,[l_q,l_a]);
s_y = sum(sumterms_y,3);
s_y = reshape(s_y,[l_q,l_a]);
s_z = sum(sumterms_z,3);
s_z = reshape(s_z,[l_q,l_a]);

if(min(s_x)+1E-12<0)
    warning('Negative sums found in GPD sum.  Setting to zero.');
    s_x(find(s_x<0))=0;
    x;
end
logEx = - 0.5.*GAMMA.^2.*Gxmx.^2.*s_x;
logEy = - 0.5.*GAMMA.^2.*Gymx.^2.*s_y;
logEz = - 0.5.*GAMMA.^2.*Gzmx.^2.*s_z;

logE = logEx+logEy+logEz;

E_r_matrix = exp(logE);
E_r=mean(E_r_matrix,2);
E=E_r;

       
 

 % Compute the Jacobian matrix; computed numerically
if(nargout>1)
    dx = 0.00001;
    J = zeros(length(E), 3);
    if nargin < 3 
         
        for i = 1:3; % compute the derivatives for all model parameters
            xpert = x;
            xpert(i) = xpert(i)*(1+dx);
            Epert = CuboidSym_GPD_TWOGSE(xpert, protocol);
            dEdx = (Epert - E)/(xpert(i)*dx);
            J(:,i) = dEdx;
        end
    else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
        
        for i = 1:length(x_deriv);
            if x_deriv(i) ~= 0  
                
                xpert = x;
                xpert(i) = xpert(i)*(1+dx);
                Epert = CuboidSym_GPD_TWOGSE(xpert, protocol);
                dEdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEdx;
            end
        end
       

    end   
 end