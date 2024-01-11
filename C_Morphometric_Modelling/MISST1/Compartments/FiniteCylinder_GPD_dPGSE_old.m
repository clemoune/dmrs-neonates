function [E,J]=FiniteCylinder_GPD_dPGSE_old(x, protocol,x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a finite cylinder compartment.
% 
% [E,J]=FiniteCylinder_GPD_dPGSE(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate consisting of parallel, 
% impermeable finite cylinders with a single radius and a diffusion protocol 
% specified in the input
% Substrate: parallel, impermeable finite cylinders with a single radius
% Diffusion pulse sequence: Double pulsed gradient spin echo (dPGSE)
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
%       protocol.grad_dirs1 - is the gradient direction of the first 
%       gradient pair for each measurement. It has size [N 3] where 
%       N is the number of measurements.
%       protocol.grad_dirs2 - is the gradient direction of the second 
%       gradient pair for each measurement. It has size [N 3] where 
%       N is the number of measurements.
%       protocol.G1 - gradient strength of the first gradient pair, size [1 N]
%       protocol.G2 - gradient strength of the second gradient pair, size [1 N]
%       protocol.delta - pulse separation of each pair, size [1 N]
%       protocol.smalldel - pulse duration of each pair, size [1 N]
%       protocol.tm - mixing time between the two gradient pairs, must be 
%       larger than protocol.smalldel, size [1 N]
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
roots_plane = protocol.roots_plane;
% Check the roots array is correct
if(abs(roots(1) - 1.8412)>0.0001)
    error('Looks like the roots array is wrong.  First value should be 1.8412, but is %f', roots(1));
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


L = 2.*R.*ecc;

% get the relevalt pusle sequence parameters from protocol
G1 = protocol.G1';
G2 = protocol.G2';
smalldel = protocol.smalldel';
delta = protocol.delta';
tm = protocol.tm';
grad_dirs1 = protocol.grad_dirs1;
grad_dirs2 = protocol.grad_dirs2;

if min(tm-smalldel) <0
    error('not yet implemented for tm < smalldel')
end

% calculate fibre direction from the specified angles
fibredir = [sin(theta).*cos(phi); sin(theta).*sin(phi); cos(theta)];

% Radial wavenumbers
GAMMA = 2.675987E8; % This is what is used throughout Wuzi.
%GAMMA = 2.6751525E8; % This is the latest best estimate of GAMMA (used in Camino)

if numel(theta) == 1
    l_q=size(grad_dirs1,1);
    l_a=numel(R);
    k_max=numel(roots);
    k_max_plane = numel(roots_plane);
    
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


    L_mat=repmat(L,[l_q 1]); % cylinder length
    L_mat=L_mat(:);
    L_mat=repmat(L_mat,[1 1 k_max_plane]);

    root_m_plane=reshape(roots_plane,[1 1 k_max_plane]);
    root_mat2_plane = repmat(root_m_plane,[l_q*l_a 1 1]).^2;
    lambda_mat_plane=pi^2.*root_mat2_plane./L_mat.^2;
    Bmat_plane = 8*L_mat.^2./root_mat2_plane.^2./pi^4;

    % Angles between gradient directions and fibre direction.
    cosTheta1 = grad_dirs1*fibredir;
    cosTheta1mx = repmat(cosTheta1,[1,l_a]);
    
    cosTheta2 = grad_dirs2*fibredir;
    cosTheta2mx = repmat(cosTheta2,[1,l_a]);
    
else
      l_q=size(grad_dirs1,1);
    l_R=numel(R);
    l_N = numel(theta);
    l_a = l_R*l_N;
    k_max=numel(roots); 
    k_max_plane = numel(roots_plane);

    if isfield(protocol,'Rweight') % all directions are repeaded for eac radius
        Rweight = protocol.Rweight;
        Rweightmx = repmat(Rweight,[l_N,1]);
        Rweightmx = Rweightmx(:)';
    else
        Rweightmx = ones(1,l_a)./l_R;
    end
    if isfield(protocol,'Nweight')
        Nweight = protocol.Nweight;
        Nweightmx = repmat(Nweight,[1,l_R]);
    else
        Nweightmx = ones(1,l_a)./l_N;
    end
    weight_matrix=repmat(Rweightmx.*Nweightmx,[l_q 1]);
    
     Rm = repmat(R,[l_N,1]);
    Rm =  Rm(:)';
    R_mat=repmat(Rm,[l_q 1]);
    R_mat=R_mat(:);
    R_mat=repmat(R_mat,[1 1 k_max]);


    root_m=reshape(roots,[1 1 k_max]);
    alpha_mat=repmat(root_m,[l_q*l_a 1 1])./R_mat;
    amSq=alpha_mat.^2;
    
    Lm = repmat(L,[l_N,1]);
    Lm =  Lm(:)';
    L_mat=repmat(Lm,[l_q 1]); % cylinder length
    L_mat=L_mat(:);
    L_mat=repmat(L_mat,[1 1 k_max_plane]);

    root_m_plane=reshape(roots_plane,[1 1 k_max_plane]);
    root_mat2_plane = repmat(root_m_plane,[l_q*l_a 1 1]).^2;
    lambda_mat_plane=pi^2.*root_mat2_plane./L_mat.^2;
    Bmat_plane = 8*L_mat.^2./root_mat2_plane.^2./pi^4;
    
        % Angles between gradient directions and fibre direction.
    cosTheta1 = grad_dirs1*fibredir;
    cosTheta1mx = repmat(cosTheta1,[1,l_R]);

    cosTheta2 = grad_dirs2*fibredir;
    cosTheta2mx = repmat(cosTheta2,[1,l_R]);   
    
    
end

deltamx=repmat(delta,[1,l_a]);
deltamx_rep = deltamx(:);
deltamx_rep = repmat(deltamx_rep,[1 1 k_max]);

smalldelmx=repmat(smalldel,[1,l_a]);
smalldelmx_rep = smalldelmx(:);
smalldelmx_rep = repmat(smalldelmx_rep,[1 1 k_max]);

tmmx=repmat(tm,[1,l_a]);
tmmx_rep = tmmx(:);
tmmx_rep = repmat(tmmx_rep,[1 1 k_max]);

G1mx=repmat(G1,[1,l_a]);
G1mx_par = G1mx.*cosTheta1mx;
G1mx_rep = G1mx(:); G1mx_rep = repmat(G1mx_rep ,[1 1 k_max]);
G1mx_par_rep = G1mx_par(:); G1mx_par_rep = repmat(G1mx_par_rep ,[1 1 k_max]);

G2mx=repmat(G2,[1,l_a]);
G2mx_par = G2mx.*cosTheta2mx;
G2mx_rep = G2mx(:); G2mx_rep = repmat(G2mx_rep ,[1 1 k_max]);
G2mx_par_rep = G2mx_par(:); G2mx_par_rep = repmat(G2mx_par_rep ,[1 1 k_max]);

CosPsi = sum(grad_dirs1.*grad_dirs2,2);
CosPsimx = repmat(CosPsi,[1,l_a]);
CosPsimx_rep = CosPsimx(:); CosPsimx_rep = repmat(CosPsimx_rep ,[1 1 k_max]);

% Parallel component - planar restriction
exp_smalldel_plane = exp(-lambda_mat_plane.*smalldelmx_rep.*dRes);
exp_delta_plane = exp(-lambda_mat_plane.*deltamx_rep.*dRes);
exp_delta_smalldel_plane = exp(-lambda_mat_plane.*(deltamx_rep-smalldelmx_rep).*dRes);
exp_tm_plane = exp(-lambda_mat_plane.*dRes.*tmmx_rep);
exp_tm_smalldel_plane = exp(-lambda_mat_plane.*dRes.*(tmmx_rep-smalldelmx_rep));

%Geometric factor B
Bmat_rep = 2./amSq./(repmat(root_m,[l_q*l_a 1 1]).^2 -1);

exp_delta = exp(-amSq.*dRes.*deltamx_rep);
exp_delta_smalldel = exp(-amSq.*dRes.*(deltamx_rep-smalldelmx_rep));
exp_smalldel = exp(-amSq.*dRes.*smalldelmx_rep);
exp_tm = exp(-amSq.*dRes.*tmmx_rep);
exp_tm_smalldel = exp(-amSq.*dRes.*(tmmx_rep-smalldelmx_rep));
am_dRes = amSq.*dRes;
am2_dRes2 = amSq.^2.*dRes.^2;

am_dRes_plane = lambda_mat_plane.*dRes;
am2_dRes2_plane = lambda_mat_plane.^2.*dRes.^2;

% % S11
% sumterms = 2.*(G1mx_perp_rep.^2+G2mx_perp_rep.^2).*(2.*(am_dRes.*smalldelmx_rep+exp_smalldel-1)./am2_dRes2);
% % S12
% sumterms = sumterms+ 2.*(G1mx_perp_rep.^2+G2mx_perp_rep.^2).*(-exp_delta_smalldel.*(1-exp_smalldel).^2)./am2_dRes2;
% % S13
% sumterms = sumterms- 4.*G1mx_perp_rep.*G2mx_perp_rep.*((1-exp_smalldel).^2.*exp_tm_smalldel.*exp_delta)./am2_dRes2;
% % S14
% sumterms = sumterms+ 2.*G1mx_perp_rep.*G2mx_perp_rep.*((1-exp_smalldel).^2.*exp_tm.*exp_delta.*exp_delta_smalldel)./am2_dRes2;
% % S23
% sumterms = sumterms+ 2*G1mx_perp_rep.*G2mx_perp_rep.*((1-exp_smalldel).^2.*exp_tm_smalldel)./am2_dRes2;
% 
%  sumterms = Bmat_rep.*sumterms;

% S11
sumterms1 = 2.*(G1mx_rep.^2+G2mx_rep.^2).*(2.*(am_dRes.*smalldelmx_rep+exp_smalldel-1)./am2_dRes2);
% S12
sumterms1 = sumterms1+ 2.*(G1mx_rep.^2+G2mx_rep.^2).*(-exp_delta_smalldel.*(1-exp_smalldel).^2)./am2_dRes2;
% S13
sumterms1 = sumterms1- 4.*G1mx_rep.*G2mx_rep.*CosPsimx_rep.*((1-exp_smalldel).^2.*exp_tm_smalldel.*exp_delta)./am2_dRes2;
% S14
sumterms1 = sumterms1+ 2.*G1mx_rep.*G2mx_rep.*CosPsimx_rep.*((1-exp_smalldel).^2.*exp_tm.*exp_delta.*exp_delta_smalldel)./am2_dRes2;
% S23
sumterms1 = sumterms1+ 2*G1mx_rep.*G2mx_rep.*CosPsimx_rep.*((1-exp_smalldel).^2.*exp_tm_smalldel)./am2_dRes2;

 sumterms1 = Bmat_rep.*sumterms1;
 
 % S11
sumterms2 = 2.*(G1mx_par_rep.^2+G2mx_par_rep.^2).*(2.*(am_dRes.*smalldelmx_rep+exp_smalldel-1)./am2_dRes2);
% S12
sumterms2 = sumterms2+ 2.*(G1mx_par_rep.^2+G2mx_par_rep.^2).*(-exp_delta_smalldel.*(1-exp_smalldel).^2)./am2_dRes2;
% S13
sumterms2 = sumterms2- 4.*G1mx_par_rep.*G2mx_par_rep.*((1-exp_smalldel).^2.*exp_tm_smalldel.*exp_delta)./am2_dRes2;
% S14
sumterms2 = sumterms2+ 2.*G1mx_par_rep.*G2mx_par_rep.*((1-exp_smalldel).^2.*exp_tm.*exp_delta.*exp_delta_smalldel)./am2_dRes2;
% S23
sumterms2 = sumterms2+ 2*G1mx_par_rep.*G2mx_par_rep.*((1-exp_smalldel).^2.*exp_tm_smalldel)./am2_dRes2;

 sumterms2 = Bmat_rep.*sumterms2;
 
 
 % planar restriction
 
  % S11
sumterms3 = 2.*(G1mx_par_rep.^2+G2mx_par_rep.^2).*(2.*(am_dRes_plane.*smalldelmx_rep+exp_smalldel_plane-1)./am2_dRes2_plane);
% S12
sumterms3 = sumterms3+ 2.*(G1mx_par_rep.^2+G2mx_par_rep.^2).*(-exp_delta_smalldel_plane.*(1-exp_smalldel_plane).^2)./am2_dRes2_plane;
% S13
sumterms3 = sumterms3- 4.*G1mx_par_rep.*G2mx_par_rep.*((1-exp_smalldel_plane).^2.*exp_tm_smalldel_plane.*exp_delta_plane)./am2_dRes2_plane;
% S14
sumterms3 = sumterms3+ 2.*G1mx_par_rep.*G2mx_par_rep.*((1-exp_smalldel_plane).^2.*exp_tm_plane.*exp_delta_plane.*exp_delta_smalldel_plane)./am2_dRes2_plane;
% S23
sumterms3 = sumterms3+ 2*G1mx_par_rep.*G2mx_par_rep.*((1-exp_smalldel_plane).^2.*exp_tm_smalldel_plane)./am2_dRes2_plane;

 sumterms3 = Bmat_plane.*sumterms3;


% testinds = find(sumterms(:,:,end)>0);
% test = sumterms(testinds,1)./sumterms(testinds,end);
% if(min(test)<1E4)
%     warning('Ratio of largest to smallest terms in GPD model sum is <1E4.  May need more terms.');
%     x;
% end
% 
% s = sum(sumterms,3);
% s = reshape(s,[l_q,l_a]);
% if(min(s)<0)
%     warning('Negative sums found in GPD sum.  Setting to zero.');
%     s(find(s<0))=0;
%     x;
% end
%disp(s.*GmxSq)

s = sum(sumterms1,3)-sum(sumterms2,3);
s = reshape(s,[l_q,l_a]);
if(min(s)<0)
    warning('Negative sums found in GPD sum.  Setting to zero.');
    s(find(s<0))=0;
    x;
end

s_plane = sum(sumterms3,3);
s_plane = reshape(s_plane,[l_q,l_a]);
if(min(s_plane)<0)
    warning('Negative sums found in GPD sum.  Setting to zero.');
    s_plane(find(s_plane<0))=0;
    x;
end

logE = -0.5*GAMMA^2.*s;
ePerp = exp(logE);

logEpar = -0.5*GAMMA^2.*s_plane;
ePar = exp(logEpar);

E_r_matrix = ePar.*ePerp;

E_r=sum(E_r_matrix.*weight_matrix,2);

E=E_r;


% Compute the Jacobian matrix; computed numerically
if(nargout>1)
    dx = 0.00001;
     J = zeros(length(E),length(x));
    if nargin < 3 
        
        for i = 1:length(x); % compute the derivatives for all model parameters
            xpert = x;
            xpert(i) = xpert(i)*(1+dx);
            Epert = FiniteCylinder_GPD_dPGSE(xpert, protocol);
            dEdx = (Epert - E)/(xpert(i)*dx);
            J(:,i) = dEdx;
        end
        
       
    else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
    
        for i = 1:length(x_deriv);
            if x_deriv(i) ~= 0                  
                xpert = x;
                xpert(i) = xpert(i)*(1+dx);
                Epert = FiniteCylinder_GPD_dPGSE(xpert, protocol);
                dEdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEdx;
            end
        end
        
       
    end   
end