function [E,J]=Cylinder_GPD_FullSTEAM(x,protocol)
% Substrate: Parallel, impermeable cylinders with one radius in a homogeneous
% background.
% Pulse sequence: STEAM with additional imaging gradients
% Signal approximation: Gaussian phase distribution.
%
% [E,J]=SynthMeasCylSingleRadGPD_FullSTEAM(x, grad_dirs, G, TM, smalldel, cG, rG, sdelc, sdelr, dcgap, fibredir, roots)
% returns the measurements E according to the model and the Jacobian J of the
% measurements with respect to the parameters.  The Jacobian does not
% include derisdelrvates with respect to the fibre direction.
%
% x is the list of model parameters in SI units:
% x(1) is the volume fraction of the intracellular space.
% x(2) is the free diffusivity of the material inside and outside the cylinders.
% x(3) is the hindered diffusivity outside the cylinders in perpendicular directions.
% x(4) is the radius of the cylinders.
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

%%(x, grad_dirs, G, cG, rG, TM, gap1, gap2, smalldel, sdelc, sdelr, fibredir, roots)


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
if(abs(roots(1) - 1.8412)>0.0001)
    error('Looks like the roots array is wrong.  First value should be 1.8412, but is %f', roots(1));
end

% Radial wavenumbers
GAMMA = 2.675987E8; % This is what is used throughout Wuzi.

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
ePar=exp(-GAMMA^2*dRes*(GD_CompPar.^2.*smalldel.^2.*tdd + GC_CompPar.^2.*sdelc.^2.*tcc + ...
    GR_CompPar.^2.*sdelr.^2.*trr + 2*GD_CompPar.*GC_CompPar.*smalldel.*sdelc.*tdc + ...
    2*GD_CompPar.*GR_CompPar.*smalldel.*sdelr.*tdr + 2*GC_CompPar.*GR_CompPar.*sdelc.*sdelr.*tcr));

% Perpendicular component
% Terms in GPD sum.
a2d = amSq*dRes;

sumnumGG = -2 + 2*exp(-a2d.*sdeldmx) - ...
     exp(-(a2d.*(gap1mx + gap2mx + 2*(sdelcmx + sdelrmx) + TMmx))) + ...
     2*exp(-a2d.*(gap1mx + gap2mx + 2*sdelcmx + 2*sdelrmx + sdeldmx + TMmx)) - ...
     exp(-(a2d.*(gap1mx + gap2mx + 2*(sdelcmx + sdelrmx + sdeldmx) + TMmx))) + ...
     2*a2d.*sdeldmx;

sumnumGGc = exp(a2d.*(-gap1mx)) + ...
     exp(a2d.*(-gap2mx)) - ...
     exp(a2d.*(-gap1mx - sdelcmx)) - ...
     exp(a2d.*(-gap2mx - sdelcmx)) - ...
     exp(a2d.*(-gap1mx - sdeldmx)) - ...
     exp(a2d.*(-gap2mx - sdeldmx)) + ...
     exp(a2d.*(-gap1mx - sdelcmx - sdeldmx)) + ...
     exp(a2d.*(-gap2mx - sdelcmx - sdeldmx)) + ...
     exp(a2d.*(-gap1mx - 2*sdelcmx - 2*sdelrmx - TMmx)) + ...
     exp(a2d.*(-gap2mx - 2*sdelcmx - 2*sdelrmx - TMmx)) - ...
     exp(a2d.*(-gap1mx - sdelcmx - 2*sdelrmx - TMmx)) - ...
     exp(a2d.*(-gap2mx - sdelcmx - 2*sdelrmx - TMmx)) - ...
     exp(a2d.*(-gap1mx - 2*sdelcmx - 2*sdelrmx - sdeldmx - TMmx)) - ...
     exp(a2d.*(-gap2mx - 2*sdelcmx - 2*sdelrmx - sdeldmx - TMmx)) + ...
     exp(a2d.*(-gap1mx - sdelcmx - 2*sdelrmx - sdeldmx - TMmx)) + ...
     exp(a2d.*(-gap2mx - sdelcmx - 2*sdelrmx - sdeldmx - TMmx));

sumnumGcGc = -2 + 2*exp(-a2d.*sdelcmx) - ...
     exp(-(a2d.*(2*sdelrmx + TMmx))) + ...
     2*exp(-a2d.*(sdelcmx + 2*sdelrmx + TMmx)) - ...
     exp(-(a2d.*(2*(sdelcmx + sdelrmx) + TMmx))) + ...
     2*a2d.*sdelcmx;

sumnumGGr = exp(a2d.*(-gap1mx - sdelcmx)) + ...
     exp(a2d.*(-gap2mx - sdelcmx)) - ...
     exp(a2d.*(-gap1mx - sdelcmx - sdelrmx)) - ...
     exp(a2d.*(-gap2mx - sdelcmx - sdelrmx)) - ...
     exp(a2d.*(-gap1mx - sdelcmx - sdeldmx)) - ...
     exp(a2d.*(-gap2mx - sdelcmx - sdeldmx)) + ...
     exp(a2d.*(-gap1mx - sdelcmx - sdelrmx - sdeldmx)) + ...
     exp(a2d.*(-gap2mx - sdelcmx - sdelrmx - sdeldmx)) + ...
     exp(a2d.*(-gap1mx - sdelcmx - 2*sdelrmx - TMmx)) + ...
     exp(a2d.*(-gap2mx - sdelcmx - 2*sdelrmx - TMmx)) - ...
     exp(a2d.*(-gap1mx - sdelcmx - sdelrmx - TMmx)) - ...
     exp(a2d.*(-gap2mx - sdelcmx - sdelrmx - TMmx)) - ...
     exp(a2d.*(-gap1mx - sdelcmx - 2*sdelrmx - sdeldmx - TMmx)) - ...
     exp(a2d.*(-gap2mx - sdelcmx - 2*sdelrmx - sdeldmx - TMmx)) + ...
     exp(a2d.*(-gap1mx - sdelcmx - sdelrmx - sdeldmx - TMmx)) + ...
     exp(a2d.*(-gap2mx - sdelcmx - sdelrmx - sdeldmx - TMmx));

sumnumGcGr = 2 - 2*exp(-a2d.*sdelcmx) - ...
     2*exp(-a2d.*sdelrmx) + ...
     2*exp(a2d.*(-sdelcmx - sdelrmx)) + ...
     2*exp(a2d.*(-2*sdelrmx - TMmx)) - ...
     2*exp(a2d.*(-sdelcmx - 2*sdelrmx - TMmx)) - ...
     2*exp(a2d.*(-sdelrmx - TMmx)) + ...
     2*exp(a2d.*(-sdelcmx - sdelrmx - TMmx));...

sumnumGrGr = -2 + 2*exp(-a2d.*sdelrmx) - ...
     exp(-(a2d.*TMmx)) + ...
     2*exp(-a2d.*(sdelrmx + TMmx)) - ...
     exp(-(a2d.*(2*sdelrmx + TMmx))) + ...
     2*a2d.*sdelrmx;

sumdenom = dRes^2*amP6.*(R_matSq.*amSq - 1);

% Check for zeros on top and bottom
%sumdenom(find(sumnum) == 0) = 1;

sumtermsGG = sumnumGG./sumdenom;
sumtermsGGc = sumnumGGc./sumdenom;
sumtermsGcGc = sumnumGcGc./sumdenom;
sumtermsGGr = sumnumGGr./sumdenom;
sumtermsGcGr = sumnumGcGr./sumdenom;
sumtermsGrGr = sumnumGrGr./sumdenom;

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

% Log of attenuation from gradients along perpendicular component of
% diffusion gradient.
logE_Comp1 = -2*GAMMA^2.*(sGG.*(GD_Comp1.^2) + sGGc.*GD_Comp1.*GC_Comp1 + ...
                                              sGcGc.*(GC_Comp1.^2) + sGGr.*GD_Comp1.*GR_Comp1 + ...
                                              sGrGr.*(GR_Comp1.^2) + sGcGr.*GC_Comp1.*GR_Comp1);

% Log of attenuation from gradients along the other perpendicular
% direction.
logE_Comp2 = -2*GAMMA^2.*(sGcGc.*(GC_Comp2.^2) + sGrGr.*(GR_Comp2.^2) + sGcGr.*GC_Comp2.*GR_Comp2);

ePerp = exp(logE_Comp1 + logE_Comp2);

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
    vgdterm=@(arg) -amSq.*arg.*exp(-dRes.*amSq.*arg);
    
    sumnumD = 2*amSq.*(delta1mx+delta2mx);
    sumnumD = sumnumD - (vgdterm(t2mx-t1mx) - vgdterm(t3mx-t1mx) - vgdterm(t3mx-t2mx) - vgdterm(delta1mx) - vgdterm(t2mx-t1mx-delta1mx) + vgdterm(t3mx-t1mx-delta1mx) - 2*vgdterm(delta2mx) - 2*vgdterm(t2mx-t1mx+delta2mx) + 2*vgdterm(t2mx-t1mx+delta2mx-delta1mx) + 2*vgdterm(t3mx-t2mx-delta2mx) - 2*vgdterm(delta3mx) + vgdterm(delta2mx+delta3mx) + vgdterm(t2mx-t1mx+delta2mx+delta3mx) - vgdterm(t2mx-t1mx+delta2mx+delta3mx-delta1mx) - 2*vgdterm(t3mx-t2mx+delta1mx-delta3mx) - vgdterm(t3mx-t1mx+delta2mx-delta3mx) - vgdterm(delta1mx+delta2mx-delta3mx) + vgdterm(t3mx-t1mx+delta1mx+delta2mx-delta3mx) + vgdterm(t3mx-t2mx+delta1mx+delta2mx-delta3mx) - vgdterm(t3mx-t2mx-delta2mx-delta3mx) + vgdterm(t3mx-t2mx+delta1mx-2*delta3mx));
    sumtermsD = sumnumD./sumdenom;

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
    Epert = SynthMeasCylSingleRadGPD_DSE(xpert, grad_dirs, G, TE, delta1, delta2, delta3, t1, t2, t3, fibredir, roots);
    dEtdr = (Epert - E)/(xpert(4)*dx);
        
    % Construct the jacobian matrix. 
    J = zeros(length(E), 4);
    J(:,1) = dEtdf;
    J(:,2) = dEtddPar;
    J(:,3) = dEtddPerp;
    J(:,4) = dEtdr;
    
end


    
    
