function [E,J]=FiniteCylinder_GPD_FullSTEAM(x,protocol)
% Substrate: Parallel, impermeable finite cylinders with one radius and length in a homogeneous
% background.
% Pulse sequence: STEAM with additional imaging gradients
% Signal approximation: Gaussian phase distribution.
%
% [E,J]=FiniteCylinder_GPD_FullSTEAM(x,protocol)
% returns the measurements E according to the model and the Jacobian J of the
% measurements with respect to the parameters.  The Jacobian does not
% include derisdelrvates with respect to the fibre direction.
%
% x is the list of model parameters in SI units:
% x(1) is the free diffusivity of the material inside
% x(2) is the radius of the cylinders.
% x(3) is the eccentricity of the cylinder.
% x(4) is the polar angle
% x(5) is the azimuth angle
%
% grad_dirs is the gradient direction for each measurement.  It has size [N
% 3] where N is the number of measurements.
%
% G, TM, smalldel, cG, rG, sdelc, sdelr, gap1 and gap2 are the
% gradient strength,
% mixing time, gradient vectors for the crusher and slice-select pulses,
% lengths of the crusher and slice select pulses and gaps between the
% first and second diffusion pulse and crusher pairs, for each
% measurement in the
% protocol.
%
% fibredir is a unit vector along the cylinder axis.  It must be in
% Cartesian coordinates [x y z]' with size [3 1].
%
% roots contains solutions to the Bessel function equation from function
% BesselJ_RootsCyl.
%
% $Id$



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
% delta = protocol.delta';
grad_dirs = protocol.grad_dirs;
cG = protocol.cG;
rG = protocol.rG;
TM = protocol.TM';
sdelc = protocol.sdelc';
sdelr = protocol.sdelr';
gap1 = protocol.gap1';
gap2 = protocol.gap2';

% Check the roots array is correct
if(abs(roots_cyl(1) - 1.8412)>0.0001)
    error('Looks like the roots array is wrong.  First value should be 1.8412, but is %f', roots_cyl(1));
end

% Radial wavenumbers
GAMMA = 2.675987E8; % This is what is used throughout Wuzi.

l_q=size(grad_dirs,1);
l_a=numel(R);
k_max=numel(roots_cyl);
k_max_plane = numel(roots_plane);

if isfield(protocol,'Rweight')
    weight_matrix=repmat(protocol.Rweight,[l_q 1]);
else
    weight_matrix = repmat(ones(1,l_a)./l_a,[l_q 1]);
end




R_mat=repmat(R,[l_q 1]);
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

% Angles between gradient directions and fibre direction.
cosTheta = grad_dirs*fibredir;
% Numerical problems occasionally take this out of range
tp = find(cosTheta>1);
tm = find(cosTheta<-1);
if(max(cosTheta(tp))>(1.01) | min(cosTheta(tm)<(-1.01)))
   error('Cos theta well out of range.  Need to check this.');
else
   cosTheta(tp) = 1;
   cosTheta(tm) = -1;
end
cosThetaSq = cosTheta.^2;
sinThetaSq = 1-cosThetaSq;
sinTheta = sqrt(sinThetaSq);
cosThetaSq_matrix=repmat(cosThetaSq,[1,l_a]);
sinThetaSq_matrix=1-cosThetaSq_matrix;


sdeldmx_part=repmat(smalldel,[1,l_a]);
sdeldmx = sdeldmx_part(:);
sdeldmx = repmat(sdeldmx,[1 1 k_max]);

sdelcmx_part=repmat(sdelc,[1,l_a]);
sdelcmx = sdelcmx_part(:);
sdelcmx = repmat(sdelcmx,[1 1 k_max]);

sdelrmx_part=repmat(sdelr,[1,l_a]);
sdelrmx = sdelrmx_part(:);
sdelrmx = repmat(sdelrmx,[1 1 k_max]);

TMmx_part=repmat(TM,[1,l_a]);
TMmx = TMmx_part(:);
TMmx = repmat(TMmx,[1 1 k_max]);

gap1mx_part=repmat(gap1,[1,l_a]);
gap1mx = gap1mx_part(:);
gap1mx = repmat(gap1mx,[1 1 k_max]);

gap2mx_part=repmat(gap2,[1,l_a]);
gap2mx = gap2mx_part(:);
gap2mx = repmat(gap2mx,[1 1 k_max]);

% Components of gradients along the perpendicular component of the
% diffusion gradient vector.
fibdirrep = repmat(fibredir',[size(grad_dirs,1) 1]);
comp1_dirs = grad_dirs - repmat(sum(grad_dirs.*fibdirrep,2),[1 3]).*fibdirrep;
% Find any measurements where the gradient and fibre directions
% are parallel
fgpar = find(abs(cosTheta) == 1);
comp1_dirs(fgpar, :) = repmat([fibredir(3) fibredir(3) fibredir(1)+fibredir(2)], [length(fgpar) 1]);
comp1_dirs = comp1_dirs./repmat(sqrt(sum(comp1_dirs.^2, 2)), [1 3]);
GD_Comp1 = repmat(G.*sinTheta,[1,l_a]);
GC_Comp1 = repmat(sum(comp1_dirs.*cG,2),[1,l_a]);
GR_Comp1 = repmat(sum(comp1_dirs.*rG,2),[1,l_a]);

% Components along a perpendicular direction also perpendicular
% to the diffusion gradient.
comp2_dirs = cross(grad_dirs, fibdirrep);
comp2_dirs(fgpar, :) = cross(grad_dirs(fgpar,:), comp1_dirs(fgpar,:));
comp2_dirs = comp2_dirs./repmat(sqrt(sum(comp2_dirs.^2, 2)), [1 3]);
GC_Comp2 = repmat(sum(comp2_dirs.*cG,2),[1,l_a]);
GR_Comp2 = repmat(sum(comp2_dirs.*rG,2),[1,l_a]);

% Components parallel to fibre
GD_CompPar = G.*cosTheta;
GC_CompPar = cG*fibredir;
GR_CompPar = rG*fibredir;

% Synthesize measurements from model



tdd = gap1 + gap2 + TM + 2*sdelc + 2*smalldel/3 + 2*sdelr;
tcc = TM + 2*sdelc/3 + 2*sdelr;
trr = TM + 2*sdelr/3;
tdc = TM + sdelc + 2*sdelr;
tdr = TM + sdelr;
tcr = TM + sdelr;

% Intra-axonal signals
% Parallel component
% ePar=exp(-GAMMA^2*dRes*(GD_CompPar.^2.*smalldel.^2.*tdd + GC_CompPar.^2.*sdelc.^2.*tcc + ...
%     GR_CompPar.^2.*sdelr.^2.*trr + 2*GD_CompPar.*GC_CompPar.*smalldel.*sdelc.*tdc + ...
%     2*GD_CompPar.*GR_CompPar.*smalldel.*sdelr.*tdr + 2*GC_CompPar.*GR_CompPar.*sdelc.*sdelr.*tcr));

% Perpendicular component
% Terms in GPD sum.
a2d_cyl = amSq_cyl*dRes;
a2d_plane = lambda_mat_plane*dRes;

sumnumGG_cyl = -2 + 2*exp(-a2d_cyl.*sdeldmx) - ...
     exp(-(a2d_cyl.*(gap1mx + gap2mx + 2*(sdelcmx + sdelrmx) + TMmx))) + ...
     2*exp(-a2d_cyl.*(gap1mx + gap2mx + 2*sdelcmx + 2*sdelrmx + sdeldmx + TMmx)) - ...
     exp(-(a2d_cyl.*(gap1mx + gap2mx + 2*(sdelcmx + sdelrmx + sdeldmx) + TMmx))) + ...
     2*a2d_cyl.*sdeldmx;
 
sumnumGG_plane = -2 + 2*exp(-a2d_plane.*sdeldmx) - ...
     exp(-(a2d_plane.*(gap1mx + gap2mx + 2*(sdelcmx + sdelrmx) + TMmx))) + ...
     2*exp(-a2d_plane.*(gap1mx + gap2mx + 2*sdelcmx + 2*sdelrmx + sdeldmx + TMmx)) - ...
     exp(-(a2d_plane.*(gap1mx + gap2mx + 2*(sdelcmx + sdelrmx + sdeldmx) + TMmx))) + ...
     2*a2d_plane.*sdeldmx; 

sumnumGGc_cyl = exp(a2d_cyl.*(-gap1mx)) + ...
     exp(a2d_cyl.*(-gap2mx)) - ...
     exp(a2d_cyl.*(-gap1mx - sdelcmx)) - ...
     exp(a2d_cyl.*(-gap2mx - sdelcmx)) - ...
     exp(a2d_cyl.*(-gap1mx - sdeldmx)) - ...
     exp(a2d_cyl.*(-gap2mx - sdeldmx)) + ...
     exp(a2d_cyl.*(-gap1mx - sdelcmx - sdeldmx)) + ...
     exp(a2d_cyl.*(-gap2mx - sdelcmx - sdeldmx)) + ...
     exp(a2d_cyl.*(-gap1mx - 2*sdelcmx - 2*sdelrmx - TMmx)) + ...
     exp(a2d_cyl.*(-gap2mx - 2*sdelcmx - 2*sdelrmx - TMmx)) - ...
     exp(a2d_cyl.*(-gap1mx - sdelcmx - 2*sdelrmx - TMmx)) - ...
     exp(a2d_cyl.*(-gap2mx - sdelcmx - 2*sdelrmx - TMmx)) - ...
     exp(a2d_cyl.*(-gap1mx - 2*sdelcmx - 2*sdelrmx - sdeldmx - TMmx)) - ...
     exp(a2d_cyl.*(-gap2mx - 2*sdelcmx - 2*sdelrmx - sdeldmx - TMmx)) + ...
     exp(a2d_cyl.*(-gap1mx - sdelcmx - 2*sdelrmx - sdeldmx - TMmx)) + ...
     exp(a2d_cyl.*(-gap2mx - sdelcmx - 2*sdelrmx - sdeldmx - TMmx));

sumnumGGc_plane = exp(a2d_plane.*(-gap1mx)) + ...
     exp(a2d_plane.*(-gap2mx)) - ...
     exp(a2d_plane.*(-gap1mx - sdelcmx)) - ...
     exp(a2d_plane.*(-gap2mx - sdelcmx)) - ...
     exp(a2d_plane.*(-gap1mx - sdeldmx)) - ...
     exp(a2d_plane.*(-gap2mx - sdeldmx)) + ...
     exp(a2d_plane.*(-gap1mx - sdelcmx - sdeldmx)) + ...
     exp(a2d_plane.*(-gap2mx - sdelcmx - sdeldmx)) + ...
     exp(a2d_plane.*(-gap1mx - 2*sdelcmx - 2*sdelrmx - TMmx)) + ...
     exp(a2d_plane.*(-gap2mx - 2*sdelcmx - 2*sdelrmx - TMmx)) - ...
     exp(a2d_plane.*(-gap1mx - sdelcmx - 2*sdelrmx - TMmx)) - ...
     exp(a2d_plane.*(-gap2mx - sdelcmx - 2*sdelrmx - TMmx)) - ...
     exp(a2d_plane.*(-gap1mx - 2*sdelcmx - 2*sdelrmx - sdeldmx - TMmx)) - ...
     exp(a2d_plane.*(-gap2mx - 2*sdelcmx - 2*sdelrmx - sdeldmx - TMmx)) + ...
     exp(a2d_plane.*(-gap1mx - sdelcmx - 2*sdelrmx - sdeldmx - TMmx)) + ...
     exp(a2d_plane.*(-gap2mx - sdelcmx - 2*sdelrmx - sdeldmx - TMmx)); 

sumnumGcGc_cyl = -2 + 2*exp(-a2d_cyl.*sdelcmx) - ...
     exp(-(a2d_cyl.*(2*sdelrmx + TMmx))) + ...
     2*exp(-a2d_cyl.*(sdelcmx + 2*sdelrmx + TMmx)) - ...
     exp(-(a2d_cyl.*(2*(sdelcmx + sdelrmx) + TMmx))) + ...
     2*a2d_cyl.*sdelcmx;

 sumnumGcGc_plane = -2 + 2*exp(-a2d_plane.*sdelcmx) - ...
     exp(-(a2d_plane.*(2*sdelrmx + TMmx))) + ...
     2*exp(-a2d_plane.*(sdelcmx + 2*sdelrmx + TMmx)) - ...
     exp(-(a2d_plane.*(2*(sdelcmx + sdelrmx) + TMmx))) + ...
     2*a2d_plane.*sdelcmx;

sumnumGGr_cyl = exp(a2d_cyl.*(-gap1mx - sdelcmx)) + ...
     exp(a2d_cyl.*(-gap2mx - sdelcmx)) - ...
     exp(a2d_cyl.*(-gap1mx - sdelcmx - sdelrmx)) - ...
     exp(a2d_cyl.*(-gap2mx - sdelcmx - sdelrmx)) - ...
     exp(a2d_cyl.*(-gap1mx - sdelcmx - sdeldmx)) - ...
     exp(a2d_cyl.*(-gap2mx - sdelcmx - sdeldmx)) + ...
     exp(a2d_cyl.*(-gap1mx - sdelcmx - sdelrmx - sdeldmx)) + ...
     exp(a2d_cyl.*(-gap2mx - sdelcmx - sdelrmx - sdeldmx)) + ...
     exp(a2d_cyl.*(-gap1mx - sdelcmx - 2*sdelrmx - TMmx)) + ...
     exp(a2d_cyl.*(-gap2mx - sdelcmx - 2*sdelrmx - TMmx)) - ...
     exp(a2d_cyl.*(-gap1mx - sdelcmx - sdelrmx - TMmx)) - ...
     exp(a2d_cyl.*(-gap2mx - sdelcmx - sdelrmx - TMmx)) - ...
     exp(a2d_cyl.*(-gap1mx - sdelcmx - 2*sdelrmx - sdeldmx - TMmx)) - ...
     exp(a2d_cyl.*(-gap2mx - sdelcmx - 2*sdelrmx - sdeldmx - TMmx)) + ...
     exp(a2d_cyl.*(-gap1mx - sdelcmx - sdelrmx - sdeldmx - TMmx)) + ...
     exp(a2d_cyl.*(-gap2mx - sdelcmx - sdelrmx - sdeldmx - TMmx));
 
 sumnumGGr_plane = exp(a2d_plane.*(-gap1mx - sdelcmx)) + ...
     exp(a2d_plane.*(-gap2mx - sdelcmx)) - ...
     exp(a2d_plane.*(-gap1mx - sdelcmx - sdelrmx)) - ...
     exp(a2d_plane.*(-gap2mx - sdelcmx - sdelrmx)) - ...
     exp(a2d_plane.*(-gap1mx - sdelcmx - sdeldmx)) - ...
     exp(a2d_plane.*(-gap2mx - sdelcmx - sdeldmx)) + ...
     exp(a2d_plane.*(-gap1mx - sdelcmx - sdelrmx - sdeldmx)) + ...
     exp(a2d_plane.*(-gap2mx - sdelcmx - sdelrmx - sdeldmx)) + ...
     exp(a2d_plane.*(-gap1mx - sdelcmx - 2*sdelrmx - TMmx)) + ...
     exp(a2d_plane.*(-gap2mx - sdelcmx - 2*sdelrmx - TMmx)) - ...
     exp(a2d_plane.*(-gap1mx - sdelcmx - sdelrmx - TMmx)) - ...
     exp(a2d_plane.*(-gap2mx - sdelcmx - sdelrmx - TMmx)) - ...
     exp(a2d_plane.*(-gap1mx - sdelcmx - 2*sdelrmx - sdeldmx - TMmx)) - ...
     exp(a2d_plane.*(-gap2mx - sdelcmx - 2*sdelrmx - sdeldmx - TMmx)) + ...
     exp(a2d_plane.*(-gap1mx - sdelcmx - sdelrmx - sdeldmx - TMmx)) + ...
     exp(a2d_plane.*(-gap2mx - sdelcmx - sdelrmx - sdeldmx - TMmx));

sumnumGcGr_cyl = 2 - 2*exp(-a2d_cyl.*sdelcmx) - ...
     2*exp(-a2d_cyl.*sdelrmx) + ...
     2*exp(a2d_cyl.*(-sdelcmx - sdelrmx)) + ...
     2*exp(a2d_cyl.*(-2*sdelrmx - TMmx)) - ...
     2*exp(a2d_cyl.*(-sdelcmx - 2*sdelrmx - TMmx)) - ...
     2*exp(a2d_cyl.*(-sdelrmx - TMmx)) + ...
     2*exp(a2d_cyl.*(-sdelcmx - sdelrmx - TMmx));...
     
 sumnumGcGr_plane = 2 - 2*exp(-a2d_plane.*sdelcmx) - ...
     2*exp(-a2d_plane.*sdelrmx) + ...
     2*exp(a2d_plane.*(-sdelcmx - sdelrmx)) + ...
     2*exp(a2d_plane.*(-2*sdelrmx - TMmx)) - ...
     2*exp(a2d_plane.*(-sdelcmx - 2*sdelrmx - TMmx)) - ...
     2*exp(a2d_plane.*(-sdelrmx - TMmx)) + ...
     2*exp(a2d_plane.*(-sdelcmx - sdelrmx - TMmx));...

sumnumGrGr_cyl = -2 + 2*exp(-a2d_cyl.*sdelrmx) - ...
     exp(-(a2d_cyl.*TMmx)) + ...
     2*exp(-a2d_cyl.*(sdelrmx + TMmx)) - ...
     exp(-(a2d_cyl.*(2*sdelrmx + TMmx))) + ...
     2*a2d_cyl.*sdelrmx;
 
 sumnumGrGr_plane = -2 + 2*exp(-a2d_plane.*sdelrmx) - ...
     exp(-(a2d_plane.*TMmx)) + ...
     2*exp(-a2d_plane.*(sdelrmx + TMmx)) - ...
     exp(-(a2d_plane.*(2*sdelrmx + TMmx))) + ...
     2*a2d_plane.*sdelrmx;

sumdenom_cyl = dRes^2*amP6_cyl.*(R_matSq.*amSq_cyl - 1);

sumdenom_plane = dRes^2.*lambda_mat_plane.^2./B_mat_plane.*2;

% Check for zeros on top and bottom
%sumdenom(find(sumnum) == 0) = 1;

sumtermsGG = sumnumGG_cyl./sumdenom_cyl;
sumtermsGGc = sumnumGGc_cyl./sumdenom_cyl;
sumtermsGcGc = sumnumGcGc_cyl./sumdenom_cyl;
sumtermsGGr = sumnumGGr_cyl./sumdenom_cyl;
sumtermsGcGr = sumnumGcGr_cyl./sumdenom_cyl;
sumtermsGrGr = sumnumGrGr_cyl./sumdenom_cyl;

sumtermsGG_plane = sumnumGG_plane./sumdenom_plane;
sumtermsGGc_plane = sumnumGGc_plane./sumdenom_plane;
sumtermsGcGc_plane = sumnumGcGc_plane./sumdenom_plane;
sumtermsGGr_plane = sumnumGGr_plane./sumdenom_plane;
sumtermsGcGr_plane = sumnumGcGr_plane./sumdenom_plane;
sumtermsGrGr_plane = sumnumGrGr_plane./sumdenom_plane;

testinds = find(sumtermsGG(:,:,end)>0);
test = sumtermsGG(testinds,1)./sumtermsGG(testinds,end);
if(min(test)<1E4)
    warning('Ratio of largest to smallest terms in GPD sum is <1E4.  May need more terms.');
    x
end

sGG = reshape(sum(sumtermsGG,3),[l_q,l_a]);
sGGc = reshape(sum(sumtermsGGc,3),[l_q,l_a]);
sGcGc = reshape(sum(sumtermsGcGc,3),[l_q,l_a]);
sGGr = reshape(sum(sumtermsGGr,3),[l_q,l_a]);
sGcGr = reshape(sum(sumtermsGcGr,3),[l_q,l_a]);
sGrGr = reshape(sum(sumtermsGrGr,3),[l_q,l_a]);

sGG_plane = reshape(sum(sumtermsGG_plane,3),[l_q,l_a]);
sGGc_plane = reshape(sum(sumtermsGGc_plane,3),[l_q,l_a]);
sGcGc_plane = reshape(sum(sumtermsGcGc_plane,3),[l_q,l_a]);
sGGr_plane = reshape(sum(sumtermsGGr_plane,3),[l_q,l_a]);
sGcGr_plane = reshape(sum(sumtermsGcGr_plane,3),[l_q,l_a]);
sGrGr_plane = reshape(sum(sumtermsGrGr_plane,3),[l_q,l_a]);

logE_CompPar = -2*GAMMA^2.*(sGG_plane.*(GD_CompPar.^2) + sGGc_plane.*GD_CompPar.*GC_CompPar + ...
                                              sGcGc_plane.*(GC_CompPar.^2) + sGGr_plane.*GD_CompPar.*GR_CompPar + ...
                                              sGrGr_plane.*(GR_CompPar.^2) + sGcGr_plane.*GC_CompPar.*GR_CompPar);



% Log of attenuation from gradients along perpendicular component of
% diffusion gradient.
logE_Comp1 = -2*GAMMA^2.*(sGG.*(GD_Comp1.^2) + sGGc.*GD_Comp1.*GC_Comp1 + ...
                                              sGcGc.*(GC_Comp1.^2) + sGGr.*GD_Comp1.*GR_Comp1 + ...
                                              sGrGr.*(GR_Comp1.^2) + sGcGr.*GC_Comp1.*GR_Comp1);

% Log of attenuation from gradients along the other perpendicular
% direction.
logE_Comp2 = -2*GAMMA^2.*(sGcGc.*(GC_Comp2.^2) + sGrGr.*(GR_Comp2.^2) + sGcGr.*GC_Comp2.*GR_Comp2);

ePerp = exp(logE_Comp1 + logE_Comp2);
ePar = exp(logE_CompPar);

E_r_matrix = ePar.*ePerp;

E_r=sum(E_r_matrix.*weight_matrix,2);
E = E_r;



% Compute the Jacobian matrix
if(nargout>1)
    
    error('Jacobian not implemented');
    
    % dE_tot/df = E_r - E_h
    dEtdf = E_r - E_h;
    
    % dE_tot/ddPar
    % dE_h/ddPar
    dEhddPar = -bValue.*cosThetaSq.*E_h;
    % dE_r/ddPar
    % dEperp/ddPar
    vgdterm=@(arg) -amSq_cyl.*arg.*exp(-dRes.*amSq_cyl.*arg);
    
    sumnumD = 2*amSq_cyl.*(delta1mx+delta2mx);
    sumnumD = sumnumD - (vgdterm(t2mx-t1mx) - vgdterm(t3mx-t1mx) - vgdterm(t3mx-t2mx) - vgdterm(delta1mx) - vgdterm(t2mx-t1mx-delta1mx) + vgdterm(t3mx-t1mx-delta1mx) - 2*vgdterm(delta2mx) - 2*vgdterm(t2mx-t1mx+delta2mx) + 2*vgdterm(t2mx-t1mx+delta2mx-delta1mx) + 2*vgdterm(t3mx-t2mx-delta2mx) - 2*vgdterm(delta3mx) + vgdterm(delta2mx+delta3mx) + vgdterm(t2mx-t1mx+delta2mx+delta3mx) - vgdterm(t2mx-t1mx+delta2mx+delta3mx-delta1mx) - 2*vgdterm(t3mx-t2mx+delta1mx-delta3mx) - vgdterm(t3mx-t1mx+delta2mx-delta3mx) - vgdterm(delta1mx+delta2mx-delta3mx) + vgdterm(t3mx-t1mx+delta1mx+delta2mx-delta3mx) + vgdterm(t3mx-t2mx+delta1mx+delta2mx-delta3mx) - vgdterm(t3mx-t2mx-delta2mx-delta3mx) + vgdterm(t3mx-t2mx+delta1mx-2*delta3mx));
    sumtermsD = sumnumD./sumdenom_cyl;

    sD = sum(sumtermsD,3);
    sD = reshape(sD,[l_q,l_a]);

    dEperpddParMatrix = -2*GAMMA^2*GmxSq.*sinThetaSq_matrix.*(sD - 2*s/dRes).*ePerp;
    dEparddParMatrix = -bValue_matrix.*cosThetaSq_matrix.*ePar;
    dErddPar = sum((dEperpddParMatrix.*ePar + dEparddParMatrix.*ePerp),2);
    dEtddPar = (1-f)*dEhddPar + f*dErddPar;
    
    % dE_tot/ddPerp
    % dE_h/ddPerp
    dEhddPerp = -bValue.*sinThetaSq.*E_h;
    dEtddPerp = (1-f)*dEhddPerp;
    
    % Perturb R and recompute the signals.  Fiddly to do analytically,
    % because the alpha_m depend on R.  Not impossible though.
    dx = 0.00001;
    xpert = x;
    xpert(4) = xpert(4)*(1+dx);
    Epert = SynthMeasCylSingleRadGPD_DSE(xpert, grad_dirs, G, TE, delta1, delta2, delta3, t1, t2, t3, fibredir, roots_cyl);
    dEtdr = (Epert - E)/(xpert(4)*dx);
        
    % Construct the jacobian matrix. 
    J = zeros(length(E), 4);
    J(:,1) = dEtdf;
    J(:,2) = dEtddPar;
    J(:,3) = dEtddPerp;
    J(:,4) = dEtdr;
    
end


    
    
