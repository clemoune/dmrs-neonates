function [Lperp,Lpar]=LogComponentsCylinder_PGSE(x, protocol)
% Substrate: Parallel, impermeable cylinders with one radius
% Pulse sequence: Pulsed gradient spin echo
% Signal approximation: Gaussian phase distribution.
%
% [Lperp,Lpar]=LogComponentsCylinder_PGSE(x, protocol)
% returns the Log of the signal when the gradient is parallel or perpendicular to the fibre direction 
% needed for watson and bingham distributions
%
% x is the list of model parameters in SI units:
% x(1) is the free diffusivity of the material inside the cylinders.
% x(2) is the radius of the cylinders.


% protocol includes all the information required to generate the signal
% protocol.grad_dirs is the gradient direction for each measurement.  It has size [N
% 3] where N is the number of measurements.
%
% protocol.G, protocol.delta and protocol.smalldel are the gradient strength, pulse separation and
% pulse length of each measurement in the protocol.  Each has
% size [N 1].
%
%
% protocol.roots_cyl contains solutions to the Bessel function equation from function
% BesselJ_RootsCyl.
%
% $Id$
roots = protocol.roots_cyl;

% Check the roots array is correct
if(abs(roots(1) - 1.8412)>0.0001)
    error('Looks like the roots array is wrong.  First value should be 1.8412, but is %f', roots(1));
end

dRes=x(1);
R=[x(2)]; 


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

R_mat=repmat(R,[l_q 1]);
R_mat=R_mat(:);
R_mat=repmat(R_mat,[1 1 k_max]);
R_matSq=R_mat.^2;

root_m=reshape(roots,[1 1 k_max]);
alpha_mat=repmat(root_m,[l_q*l_a 1 1])./R_mat;
amSq=alpha_mat.^2;
amP6=amSq.^3;


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
logEPar=-modQSq_matrix.*difftime_matrix*dRes;


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

Lperp = -2*GAMMA^2*GmxSq.*s;
Lperp(G==0) = 0;
Lpar = logEPar;


