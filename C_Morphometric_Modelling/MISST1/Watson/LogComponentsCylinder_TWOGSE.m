function [Lperp,Lpar]=LogComponentsCylinder_TWOGSE(x, protocol)
% Substrate: Parallel, impermeable cylinders with one radius 
% Pulse sequence: trapezoidal wave oscillating gradient spin echo with angular
% frequency omega - must have an integer number of lobes
% Signal approximation: Gaussian phase distribution.
%
% [Lperp,Lpar]=LogComponentsCylinder_TWOGSE(x, protocol)
% returns the Log of the signal when the gradient is parallel or perpendicular to the fibre direction 
% needed for watson and bingham distributions 
%
% x is the list of model parameters in SI units:
% x(1) is the free diffusivity of the material inside and outside the cylinders.
% x(2) is the radius of the cylinders.

%
% protocol.grad_dirs is the gradient direction for each measurement.  It has size [N
% 3] where N is the number of measurements.
%
% protocol.G, protocol.delta, protocol.smalldel and protocol.omega are the gradient strength, pulse separation,
% pulse length and gradient angular frequency of each measurement in the protocol.  Each has
% size [N 1].
%
%
% protocol.roots_cyl contains solutions to the Bessel function equation from function
% BesselJ_RootsCyl.
% $Id$

roots = protocol.roots_cyl;

% Check the roots array is correct
if(abs(roots(1) - 1.8412)>0.0001)
    error('Looks like the roots array is wrong.  First value should be 1.8412, but is %f', roots(1));
end

dRes=x(1);
R=x(2); 


G = protocol.G';
smalldel = protocol.smalldel';
delta = protocol.delta';
omega = protocol.omega';
grad_dirs = protocol.grad_dirs;
slew_rate = protocol.slew_rate';



GAMMA = 2.675987E8;

l_q=size(grad_dirs,1);
l_a=numel(R);
k_max=numel(roots);

R_mat=repmat(R,[l_q 1]);
R_mat=R_mat(:);
R_mat=repmat(R_mat,[1 1 k_max]);


root_m=reshape(roots,[1 1 k_max]);
alpha_mat=repmat(root_m,[l_q*l_a 1 1])./R_mat;
amSq=alpha_mat.^2;


%Geometric factor B
Bmat_rep = 2./amSq./(repmat(root_m,[l_q*l_a 1 1]).^2 -1);


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
Lpar=-bval_matrix.*dRes;

% Perpendicular component

%precomputing values that appear often in the expresions:
% rearrange all formulas to have only negative exponentials
exp_NT = exp(-amSq.*dRes.*NTmx_rep./(2.*freqmx_rep));
exp_1 = exp(-amSq.*dRes./(2.*freqmx_rep));
exp_smalldel_NT = exp(-amSq.*dRes.*(smalldelmx_rep-NTmx_rep./(2.*freqmx_rep)));
exp_delta_NT = exp(-amSq.*dRes.*(deltamx_rep-NTmx_rep./(2.*freqmx_rep)));
exp_delta = exp(-amSq.*dRes.*deltamx_rep);
exp_delta_smalldel = exp(-amSq.*dRes.*(deltamx_rep-smalldelmx_rep));
exp_smalldel = exp(-amSq.*dRes.*smalldelmx_rep);
exp_rt = exp(-amSq.*dRes.*rtmx_rep);
exp_1_rt = exp(-amSq.*dRes.*(1./(2.*freqmx_rep)-rtmx_rep));
exp_1_2rt = exp(-amSq.*dRes.*(1./(2.*freqmx_rep)-2.*rtmx_rep));
am_dRes = amSq.*dRes;
am2_dRes2 = amSq.^2.*dRes.^2;
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

Lperp = -0.5.*GAMMA^2*GmxSq.*s;
Lperp(G==0) = 0;

