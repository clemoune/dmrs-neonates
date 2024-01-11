function [E,J]=FiniteCylinder_GPD_DODE(x, protocol, x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a finite cylinder compartment.
% 
% [E,J]=FiniteCylinder_GPD_DODE(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate consisting of parallel, 
% impermeable finite cylinders with a single radius and a diffusion protocol 
% specified in the input
% Substrate: parallel, finite impermeable cylinders with a single radius
% Diffusion pulse sequence: Square wave oscillating gradient spin echo (SWOGSE)
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
%       protocol.Nosc - number of periods, size [1 N]
%       protocol.roots_cyl - solutions to the Bessel function equation from 
%       function BesselJ_RootsCyl.m
%       optional: protocol.phase - phase of the gradient waveform, size [1 N].
%       0 if not specified.
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

G1 = protocol.G1';
G2 = protocol.G2';
smalldel = protocol.smalldel';
tm = protocol.tm';
Nosc = protocol.Nosc';
grad_dirs1 = protocol.grad_dirs1;
grad_dirs2 = protocol.grad_dirs2;

% calculate fibre direction from the specified angles
fibredir = [sin(theta).*cos(phi); sin(theta).*sin(phi); cos(theta)];

GAMMA = 2.675987E8;

if numel(theta) == 1
    l_q=size(grad_dirs1,1);
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

    root_m=reshape(roots_cyl,[1 1 k_max]);
    alpha_mat=repmat(root_m,[l_q*l_a 1 1])./R_mat;
    amSq_cyl=alpha_mat.^2;  
    Bmat_rep = 2./amSq_cyl./(repmat(root_m,[l_q*l_a 1 1]).^2 -1);


    L_mat=repmat(L,[l_q 1]); % cylinder length
    L_mat=L_mat(:);
    L_mat=repmat(L_mat,[1 1 k_max_plane]);

    root_m_plane=reshape(roots_plane,[1 1 k_max_plane]);
    root_mat2_plane = repmat(root_m_plane,[l_q*l_a 1 1]).^2;
    amSq_plane=pi^2.*root_mat2_plane./L_mat.^2;
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
    k_max=numel(roots_cyl); 
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


    root_m=reshape(roots_cyl,[1 1 k_max]);
    alpha_mat=repmat(root_m,[l_q*l_a 1 1])./R_mat;
    amSq_cyl=alpha_mat.^2;
    Bmat_rep = 2./amSq_cyl./(repmat(root_m,[l_q*l_a 1 1]).^2 -1);
    
    Lm = repmat(L,[l_N,1]);
    Lm =  Lm(:)';
    L_mat=repmat(Lm,[l_q 1]); % cylinder length
    L_mat=L_mat(:);
    L_mat=repmat(L_mat,[1 1 k_max_plane]);

    root_m_plane=reshape(roots_plane,[1 1 k_max_plane]);
    root_mat2_plane = repmat(root_m_plane,[l_q*l_a 1 1]).^2;
    amSq_plane=pi^2.*root_mat2_plane./L_mat.^2;
    Bmat_plane = 8*L_mat.^2./root_mat2_plane.^2./pi^4;
    
        % Angles between gradient directions and fibre direction.
    cosTheta1 = grad_dirs1*fibredir;
    cosTheta1mx = repmat(cosTheta1,[1,l_R]);

    cosTheta2 = grad_dirs2*fibredir;
    cosTheta2mx = repmat(cosTheta2,[1,l_R]);   
    
    
end


% Angles between gradient directions and fibre direction.
G1mx=repmat(G1,[1,l_a]);
G1mx_par = G1mx.*cosTheta1mx;
G1mx_rep = G1mx(:); G1mx_rep = repmat(G1mx_rep ,[1 1 k_max]);
G1mx_par_rep = G1mx_par(:); G1mx_par_rep = repmat(G1mx_par_rep ,[1 1 k_max]);

G2mx=repmat(-G2,[1,l_a]);
G2mx_par = G2mx.*cosTheta2mx;
G2mx_rep = G2mx(:); G2mx_rep = repmat(G2mx_rep ,[1 1 k_max]);
G2mx_par_rep = G2mx_par(:); G2mx_par_rep = repmat(G2mx_par_rep ,[1 1 k_max]);

CosPsi = sum(grad_dirs1.*grad_dirs2,2);
CosPsimx = repmat(CosPsi,[1,l_a]);
CosPsimx_rep = CosPsimx(:); CosPsimx_rep = repmat(CosPsimx_rep ,[1 1 k_max]);

tmmx=repmat(tm,[1,l_a]);
tmmx_rep = tmmx(:);
tmmx_rep = repmat(tmmx_rep,[1 1 k_max]);

smalldelmx=repmat(smalldel,[1,l_a]);
smalldelmx_rep = smalldelmx(:);
smalldelmx_rep = repmat(smalldelmx_rep,[1 1 k_max]);

NTmx=repmat(2*Nosc,[1,l_a]);
NTmx_rep = NTmx(:);
NTmx_rep = repmat(NTmx_rep,[1 1 k_max]);


niu = Nosc./smalldel;

freqmx=repmat(niu,[1,l_a]);
freqmx_rep = freqmx(:);
freqmx_rep = repmat(freqmx_rep,[1 1 k_max]);

    
    if ~isfield(protocol,'mirror') || protocol.mirror == 0 % the 2nd gradient is the negative of the first one
              
        % Find restricted signal decay
        % Parallel component
        exp_NT_plane = exp(-amSq_plane.*dRes.*NTmx_rep./(2.*freqmx_rep));
        exp_1_plane = exp(-amSq_plane.*dRes./(2.*freqmx_rep));
        exp_smalldel_NT_plane = exp(-amSq_plane.*dRes.*(smalldelmx_rep-NTmx_rep./(2.*freqmx_rep)));
        exp_tm_smalldel_NT_plane = exp(-amSq_plane.*dRes.*(tmmx_rep+smalldelmx_rep-NTmx_rep./(2.*freqmx_rep)));
        exp_tm_smalldel_plane = exp(-amSq_plane.*dRes.*(tmmx_rep+smalldelmx_rep));
        exp_tm_plane = exp(-amSq_plane.*dRes.*(tmmx_rep));
        exp_smalldel_plane = exp(-amSq_plane.*dRes.*smalldelmx_rep);
        am_dRes_plane = amSq_plane.*dRes;
        am2_dRes2_plane = amSq_plane.^2.*dRes.^2;
        sgn = (-1).^NTmx_rep;

        % corresponds to the multiplication of the full oscillations with themselves
        sumterms_plane_11 = 2.*(-1 + exp_smalldel_NT_plane + NTmx_rep.*(-1+exp_1_plane) + am_dRes_plane.*smalldelmx_rep)./am2_dRes2_plane;

        % corresponds to the multiplication of the full oscillations with
        % other full oscillations from the same gradient
        sumterms_plane_11 = sumterms_plane_11 - 2.* sgn.*(1-exp_1_plane).^2./am2_dRes2_plane./(1+exp_1_plane).^2.*...
            (exp_NT_plane + sgn.*(NTmx_rep-1)+sgn.*NTmx_rep.*exp_1_plane);

        % corresponds to the multiplication of the full oscillations with
        % other full oscillations from different gradients
        sumterms_plane_12 = - 2.*exp_tm_smalldel_NT_plane.*(1-exp_1_plane).^2.*(1-sgn.*exp_NT_plane).*(sgn-exp_NT_plane)./...
            (1+exp_1_plane).^2./am2_dRes2_plane;

    
        sumterms_plane = (G1mx_par_rep.^2+G2mx_par_rep.^2).*sumterms_plane_11;
        sumterms_plane = sumterms_plane +(G1mx_par_rep.*G2mx_par_rep).*sumterms_plane_12;
        sumterms_plane = Bmat_plane.*sumterms_plane;

        testinds_plane = find(sumterms_plane(:,:,end)>0);
        test_plane = sumterms_plane(testinds_plane,1)./sumterms_plane(testinds_plane,end);
        if(min(abs(test_plane))<1E2)
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

        logEpar = -0.5.*GAMMA^2.*s_plane;
        ePar = exp(logEpar);



        
        % Perpendicular component

        %precomputing values that appear often in the expresions; rearrange all formulas to have only negative exponentials
        exp_NT_cyl = exp(-amSq_cyl.*dRes.*NTmx_rep./(2.*freqmx_rep));
        exp_1_cyl = exp(-amSq_cyl.*dRes./(2.*freqmx_rep));
        exp_smalldel_NT_cyl = exp(-amSq_cyl.*dRes.*(smalldelmx_rep-NTmx_rep./(2.*freqmx_rep)));
        exp_tm_smalldel_NT_cyl = exp(-amSq_cyl.*dRes.*(tmmx_rep+smalldelmx_rep-NTmx_rep./(2.*freqmx_rep)));
        exp_tm_smalldel_cyl = exp(-amSq_cyl.*dRes.*(tmmx_rep+smalldelmx_rep));
        exp_tm_cyl = exp(-amSq_cyl.*dRes.*(tmmx_rep));
        exp_smalldel_cyl = exp(-amSq_cyl.*dRes.*smalldelmx_rep);
        am_dRes_cyl = amSq_cyl.*dRes;
        am2_dRes2_cyl = amSq_cyl.^2.*dRes.^2;
        sgn = (-1).^NTmx_rep;

        % corresponds to the multiplication of the full oscillations with themselves
        sumterms_cyl11 = 2.*(-1 + exp_smalldel_NT_cyl + NTmx_rep.*(-1+exp_1_cyl) + am_dRes_cyl.*smalldelmx_rep)./am2_dRes2_cyl;

        % corresponds to the multiplication of the full oscillations with
        % other full oscillations from the same gradient
        sumterms_cyl11 = sumterms_cyl11 - 2.* sgn.*(1-exp_1_cyl).^2./am2_dRes2_cyl./(1+exp_1_cyl).^2.*...
            (exp_NT_cyl + sgn.*(NTmx_rep-1)+sgn.*NTmx_rep.*exp_1_cyl);

        % corresponds to the multiplication of the full oscillations with
        % other full oscillations from different gradients
        sumterms_cyl12 =  - 2.*exp_tm_smalldel_NT_cyl.*(1-exp_1_cyl).^2.*(1-sgn.*exp_NT_cyl).*(sgn-exp_NT_cyl)./...
            (1+exp_1_cyl).^2./am2_dRes2_cyl;

        sumterms1_cyl =(G1mx_rep.^2+G2mx_rep.^2).*sumterms_cyl11;   
        sumterms1_cyl = sumterms1_cyl+(G1mx_rep.*G2mx_rep.*CosPsimx_rep).*sumterms_cyl12; 
        
        sumterms1_cyl = Bmat_rep.*sumterms1_cyl;
        
        sumterms2_cyl = (G1mx_par_rep.^2+G2mx_par_rep.^2).*sumterms_cyl11; 
        sumterms2_cyl = sumterms2_cyl+ (G1mx_par_rep.*G2mx_par_rep).*sumterms_cyl12;  
        
        sumterms2_cyl = Bmat_rep.*sumterms2_cyl;

        s_cyl = sum(sumterms1_cyl,3)-sum(sumterms2_cyl,3);
        s_cyl = reshape(s_cyl,[l_q,l_a]);
        if(min(s_cyl)+1E-12<0)
            warning('Negative sums found in GPD sum.  Setting to zero.');
            s_cyl(find(s_cyl<0))=0;
            x;
        end
        
        logE = -0.5.*GAMMA^2.*s_cyl;
        ePerp = exp(logE); 

       
     else % the 2nd gradient is the negative of the first one reflected 
         
       error('Not yet implemented');

%       % Parallel component
%          %precomputing values that appear often in the expresions:
%         % rearrange all formulas to have only negative exponentials
%         exp_NT_plane = exp(-amSq_plane.*dRes.*NTmx_rep./(2.*freqmx_rep));
%         exp_1_plane = exp(-amSq_plane.*dRes./(2.*freqmx_rep));
%         exp_smalldel_NT_plane = exp(-amSq_plane.*dRes.*(smalldelmx_rep-NTmx_rep./(2.*freqmx_rep)));
%         exp_tm_smalldel_NT_plane = exp(-amSq_plane.*dRes.*(tmmx_rep-NTmx_rep./(2.*freqmx_rep)));
%         exp_tm_smalldel_plane = exp(-amSq_plane.*dRes.*tmmx_rep);
%         exp_tm_plane = exp(-amSq_plane.*dRes.*(tmmx_rep-smalldelmx_rep));
%         am_dRes_plane = amSq_plane.*dRes;
%         am2_dRes2_plane = amSq_plane.^2.*dRes.^2;
%         sgn = (-1).^NTmx_rep;
% 
%           % corresponds to the multiplication of the full oscillations with themselves
%         sumterms_plane = 4.*(-1 + exp_smalldel_NT_plane + NTmx_rep.*(-1+exp_1_plane) + am_dRes_plane.*smalldelmx_rep)./am2_dRes2_plane;
%         
%         % corresponds to the multiplication of the full oscillations with
%         % other full oscillations from the same gradient
%         sumterms_plane = sumterms_plane - 4.* sgn.*(1-exp_1_plane).^2./am2_dRes2_plane./(1+exp_1_plane).^2.*...
%             (exp_NT_plane + sgn.*(NTmx_rep-1)+sgn.*NTmx_rep.*exp_1_plane);
% 
%         % corresponds to the multiplication of the full oscillations with
%         % other full oscillations from different gradients
%         sumterms_plane = sumterms_plane + 2.*sgn.*exp_tm_smalldel_NT_plane.*exp_smalldel_NT_plane.*(1-exp_1_plane).^2.*(sgn.*exp_NT_plane-1).*(sgn-exp_NT_plane)./...
%             (1+exp_1_plane).^2./am2_dRes2_plane;
% 
%          % corresponds to the multiplication of the partial oscillations with themselves
%         sumterms_plane = sumterms_plane + 2.*exp_tm_plane./am2_dRes2_plane.*(exp_smalldel_NT_plane-1).*...
%             (1-exp_smalldel_NT_plane);
% 
%         % corresponds to the multiplication of the partial oscillations
%         % with full oscillations
% sumterms_plane = sumterms_plane -2./am2_dRes2_plane./(1+exp_1_plane).*(1-exp_smalldel_NT_plane).*(-1+exp_1_plane).*(sgn.*exp_NT_plane-1)+...
%       2./am2_dRes2_plane./(1+exp_1_plane).*(-1+exp_1_plane).*(sgn.*exp_NT_plane-1).*(1-exp_smalldel_NT_plane).*exp_tm_smalldel_NT_plane-...
%       2.*sgn./am2_dRes2_plane./(1+exp_1_plane).*(1-exp_1_plane).*(exp_smalldel_NT_plane-1).*(-exp_tm_smalldel_plane+sgn.*(exp_tm_smalldel_NT_plane))+...
%       2.*sgn./am2_dRes2_plane./(1+exp_1_plane).*(1-exp_1_plane).*(exp_smalldel_NT_plane-1).*(sgn-exp_NT_plane);
% 
%         sumterms_plane = Bmat_rep_plane.*sumterms_plane;
% 
%         testinds_plane = find(sumterms_plane(:,:,end)>0);
%         test_plane = sumterms_plane(testinds_plane,1)./sumterms_plane(testinds_plane,end);
% 
%         if(min(test_plane)<1E4)
%             warning('Ratio of largest to smallest terms in GPD model sum is <1E4.  May need more terms.');
%             x;
%         end
% 
%         s_plane = sum(sumterms_plane,3);
% 
%         s_plane = reshape(s_plane,[l_q,l_a]);
% 
%         if(min(s_plane)<0)
%             warning('Negative sums found in GPD sum.  Setting to zero.');
%             s_plane(find(s_plane<0))=0;
%             x;
%         end 
%         logEPar=-0.5.*GAMMA^2*GmxSq.*cosThetaSq_matrix.*s_plane;
%         ePar=exp(logEPar);
% 
%         % Perpendicular component
% 
%         %precomputing values that appear often in the expresions:
%         % rearrange all formulas to have only negative exponentials
%         exp_NT_cyl = exp(-amSq_cyl.*dRes.*NTmx_rep./(2.*freqmx_rep));
%         exp_1_cyl = exp(-amSq_cyl.*dRes./(2.*freqmx_rep));
%         exp_smalldel_NT_cyl = exp(-amSq_cyl.*dRes.*(smalldelmx_rep-NTmx_rep./(2.*freqmx_rep)));
%         exp_tm_smalldel_NT_cyl = exp(-amSq_cyl.*dRes.*(tmmx_rep-NTmx_rep./(2.*freqmx_rep)));
%         exp_tm_smalldel_cyl = exp(-amSq_cyl.*dRes.*tmmx_rep);
%         exp_tm_cyl = exp(-amSq_cyl.*dRes.*(tmmx_rep-smalldelmx_rep));
%         am_dRes_cyl = amSq_cyl.*dRes;
%         am2_dRes2_cyl = amSq_cyl.^2.*dRes.^2;
%         sgn = (-1).^NTmx_rep;
% 
%           % corresponds to the multiplication of the full oscillations with themselves
%         sumterms_cyl = 4.*(-1 + exp_smalldel_NT_cyl + NTmx_rep.*(-1+exp_1_cyl) + am_dRes_cyl.*smalldelmx_rep)./am2_dRes2_cyl;
%         
%         % corresponds to the multiplication of the full oscillations with
%         % other full oscillations from the same gradient
%         sumterms_cyl = sumterms_cyl - 4.* sgn.*(1-exp_1_cyl).^2./am2_dRes2_cyl./(1+exp_1_cyl).^2.*...
%             (exp_NT_cyl + sgn.*(NTmx_rep-1)+sgn.*NTmx_rep.*exp_1_cyl);
% 
%         % corresponds to the multiplication of the full oscillations with
%         % other full oscillations from different gradients
%         sumterms_cyl = sumterms_cyl + 2.*sgn.*exp_tm_smalldel_NT_cyl.*exp_smalldel_NT_cyl.*(1-exp_1_cyl).^2.*(sgn.*exp_NT_cyl-1).*(sgn-exp_NT_cyl)./...
%             (1+exp_1_cyl).^2./am2_dRes2_cyl;
% 
%          % corresponds to the multiplication of the partial oscillations with themselves
%         sumterms_cyl = sumterms_cyl + 2.*exp_tm_cyl./am2_dRes2_cyl.*(exp_smalldel_NT_cyl-1).*...
%             (1-exp_smalldel_NT_cyl);
% 
%         % corresponds to the multiplication of the partial oscillations
%         % with full oscillations
% sumterms_cyl = sumterms_cyl -2./am2_dRes2_cyl./(1+exp_1_cyl).*(1-exp_smalldel_NT_cyl).*(-1+exp_1_cyl).*(sgn.*exp_NT_cyl-1)+...
%       2./am2_dRes2_cyl./(1+exp_1_cyl).*(-1+exp_1_cyl).*(sgn.*exp_NT_cyl-1).*(1-exp_smalldel_NT_cyl).*exp_tm_smalldel_NT_cyl-...
%       2.*sgn./am2_dRes2_cyl./(1+exp_1_cyl).*(1-exp_1_cyl).*(exp_smalldel_NT_cyl-1).*(-exp_tm_smalldel_cyl+sgn.*(exp_tm_smalldel_NT_cyl))+...
%       2.*sgn./am2_dRes2_cyl./(1+exp_1_cyl).*(1-exp_1_cyl).*(exp_smalldel_NT_cyl-1).*(sgn-exp_NT_cyl);
% 
%         sumterms_cyl = Bmat_rep_cyl.*sumterms_cyl;
% 
%         testinds_cyl = find(sumterms_cyl(:,:,end)>0);
%         test_cyl = sumterms_cyl(testinds_cyl,1)./sumterms_cyl(testinds_cyl,end);
% 
%         if(min(test_cyl)<1E4)
%             warning('Ratio of largest to smallest terms in GPD model sum is <1E4.  May need more terms.');
%             x;
%         end
% 
%         s_cyl = sum(sumterms_cyl,3);
% 
%         s_cyl = reshape(s_cyl,[l_q,l_a]);
% 
%         if(min(s_cyl)<0)
%             warning('Negative sums found in GPD sum.  Setting to zero.');
%             s_cyl(find(s_cyl<0))=0;
%             x;
%         end 
% 
%         logE = -0.5.*GAMMA^2*GmxSq.*sinThetaSq_matrix.*s_cyl;
%         ePerp = exp(logE);
%         E_r_matrix = ePar.*ePerp; 
%         E_r=sum(E_r_matrix.*weight_matrix,2);
%         E=E_r;
     end    
E_r_matrix = ePar.*ePerp;

E_r=sum(E_r_matrix.*weight_matrix,2);

E=E_r;
 % Compute the Jacobian matrix; computed numerically
if(nargout>1)
    dx = 0.00001;
    J = zeros(length(E), length(x));
    if nargin < 3
         
        for i = 1:length(x); % compute the derivatives for all model parameters
            xpert = x;
            xpert(i) = xpert(i)*(1+dx);
            Epert = FiniteCylinder_GPD_DODE(xpert, protocol);
            dEdx = (Epert - E)/(xpert(i)*dx);
            J(:,i) = dEdx;
        end

       
    else  % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
      
        for i = 1:length(x_deriv);
            if x_deriv(i) ~= 0  
              
                xpert = x;
                xpert(i) = xpert(i)*(1+dx);
                Epert = FiniteCylinder_GPD_DODE(xpert, protocol);
                dEdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEdx;
            end
        end
    end   
 end

