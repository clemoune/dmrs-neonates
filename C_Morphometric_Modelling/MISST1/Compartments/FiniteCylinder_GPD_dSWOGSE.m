function [E,J]=FiniteCylinder_GPD_dSWOGSE(x, protocol, x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a finite cylinder compartment.
% 
% [E,J]=FiniteCylinder_GPD_dSWOGSE(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate consisting of parallel, 
% impermeable finite cylinders with a single radius and a diffusion protocol 
% specified in the input
% Substrate: parallel, impermeable finite cylinders with a single radius
% Diffusion pulse sequence: Double square wave oscillating gradients (dSWOGSE)
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
%       protocol.omega - gradient angular frequency, size [1 N]
%       protocol.tm - mixing time between the two gradient pairs, must be 
%       larger than protocol.smalldel, size [1 N]
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

G1 = protocol.G1';
G2 = protocol.G2';
smalldel = protocol.smalldel';
delta = protocol.delta';
omega = protocol.omega';
tm = protocol.tm';
grad_dirs1 = protocol.grad_dirs1;
grad_dirs2 = protocol.grad_dirs2;


% calculate fibre direction from the specified angles
fibredir = [sin(theta).*cos(phi); sin(theta).*sin(phi); cos(theta)];

GAMMA = 2.675987E8;


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

G2mx=repmat(-G2,[1,l_a]);
G2mx_par = G2mx.*cosTheta2mx;
G2mx_rep = G2mx(:); G2mx_rep = repmat(G2mx_rep ,[1 1 k_max]);
G2mx_par_rep = G2mx_par(:); G2mx_par_rep = repmat(G2mx_par_rep ,[1 1 k_max]);

CosPsi = sum(grad_dirs1.*grad_dirs2,2);
CosPsimx = repmat(CosPsi,[1,l_a]);
CosPsimx_rep = CosPsimx(:); CosPsimx_rep = repmat(CosPsimx_rep ,[1 1 k_max]);

%Geometric factor B
Bmat_rep = 2./amSq./(repmat(root_m,[l_q*l_a 1 1]).^2 -1);

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

if ~isfield(protocol,'phase')  
   
    % the number of integer half periods:
  
    NTmx_rep = floor(2.*smalldelmx_rep.*freqmx_rep+0.00000000001);
    
    if ~isfield(protocol,'mirror') || protocol.mirror == 0 % the 2nd gradient is the negative of the first one
              
       
        % Parallel component - planar restriction

        exp_NT_plane = exp(-lambda_mat_plane.*dRes.*NTmx_rep./(2.*freqmx_rep));
        exp_1_plane = exp(-lambda_mat_plane.*dRes./(2.*freqmx_rep));
        exp_smalldel_NT_plane  = exp(-lambda_mat_plane.*dRes.*(smalldelmx_rep-NTmx_rep./(2.*freqmx_rep)));
        exp_delta_NT_plane  = exp(-lambda_mat_plane.*dRes.*(deltamx_rep-NTmx_rep./(2.*freqmx_rep)));
        exp_delta_tm_NT_plane  = exp(-lambda_mat_plane.*dRes.*(deltamx_rep+tmmx_rep-NTmx_rep./(2.*freqmx_rep)));
        exp_tm_NT_plane  = exp(-lambda_mat_plane.*dRes.*(tmmx_rep-NTmx_rep./(2.*freqmx_rep)));
        exp_delta_plane  = exp(-lambda_mat_plane.*dRes.*deltamx_rep);
        exp_tm_plane  = exp(-lambda_mat_plane.*dRes.*tmmx_rep);
        exp_delta_smalldel_plane  = exp(-lambda_mat_plane.*dRes.*(deltamx_rep-smalldelmx_rep));
        exp_tm_smalldel_plane  = exp(-lambda_mat_plane.*dRes.*(tmmx_rep-smalldelmx_rep));
        exp_smalldel_plane  = exp(-lambda_mat_plane.*dRes.*smalldelmx_rep);
        am_dRes_plane  = lambda_mat_plane.*dRes;
        am2_dRes2_plane  = lambda_mat_plane.^2.*dRes.^2;
        
        % Perpendicular component

        %precomputing values that appear often in the expresions; rearrange all formulas to have only negative exponentials
        exp_NT = exp(-amSq.*dRes.*NTmx_rep./(2.*freqmx_rep));
        exp_1 = exp(-amSq.*dRes./(2.*freqmx_rep));
        exp_smalldel_NT = exp(-amSq.*dRes.*(smalldelmx_rep-NTmx_rep./(2.*freqmx_rep)));
        exp_delta_NT = exp(-amSq.*dRes.*(deltamx_rep-NTmx_rep./(2.*freqmx_rep)));
        exp_delta_tm_NT = exp(-amSq.*dRes.*(deltamx_rep+tmmx_rep-NTmx_rep./(2.*freqmx_rep)));
        exp_tm_NT = exp(-amSq.*dRes.*(tmmx_rep-NTmx_rep./(2.*freqmx_rep)));
        exp_delta = exp(-amSq.*dRes.*deltamx_rep);
        exp_tm = exp(-amSq.*dRes.*tmmx_rep);
        exp_delta_smalldel = exp(-amSq.*dRes.*(deltamx_rep-smalldelmx_rep));
        exp_tm_smalldel = exp(-amSq.*dRes.*(tmmx_rep-smalldelmx_rep));
        exp_smalldel = exp(-amSq.*dRes.*smalldelmx_rep);
        am_dRes = amSq.*dRes;
        am2_dRes2 = amSq.^2.*dRes.^2;
        sgn = (-1).^NTmx_rep;
        
         % corresponds to the multiplication of the full and partial oscillations with themselves
        sum11 = 4.*(-1 + exp_smalldel_NT + NTmx_rep.*(-1+exp_1) + am_dRes.*smalldelmx_rep)./am2_dRes2;
        
        % corresponds to the multiplication of the full oscillations with
        % other full oscillations from the same gradient
        sum11 = sum11- 4.* sgn.*(1-exp_1).^2./am2_dRes2./(1+exp_1).^2.*...
            (exp_NT + sgn.*(NTmx_rep-1)+sgn.*NTmx_rep.*exp_1); 
        % corresponds to the multiplication of the full oscillations with
        % other full oscillations from different gradients
        sum11 = sum11 + 2.*exp_delta_NT.*(1-exp_1).^2.*(1-sgn.*exp_NT).*(sgn-exp_NT)./...
            (1+exp_1).^2./am2_dRes2; % 1&2, 3&4
        % corresponds to the multiplication of the partial oscillations
        % from the first pulse with the incomplete oscillations from the
        % other pulses
        sum11 = sum11 + 2.*exp_delta_smalldel./am2_dRes2.*(exp_smalldel-exp_NT).*...
            (1-exp_smalldel_NT); % 1&2, 3&4
        % corresponds to the multiplication of the partial oscillations
        % with full oscillations                
        sum11 = sum11 + 2.*sgn.*(1-exp_1)./(1+exp_1)./am2_dRes2.*(1-exp_smalldel_NT).*...
              (-exp_delta.*exp_NT+sgn.*exp_delta_smalldel.*exp_NT+2.* exp_NT-exp_delta_smalldel-...
              2.*sgn+ sgn.*exp_delta); % 1&2, 3&4          
        
        
        % corresponds to the multiplication of the full oscillations with
        % other full oscillations from different gradients
        sum12 = - 2.*2.*exp_delta_tm_NT.*(1-exp_1).^2.*(1-sgn.*exp_NT).*(sgn-exp_NT)./...
            (1+exp_1).^2./am2_dRes2; % 1%3, 2&4
        sum12 = sum12 + 2.*exp_delta_tm_NT.*exp_delta.*(1-exp_1).^2.*(1-sgn.*exp_NT).*(sgn-exp_NT)./...
            (1+exp_1).^2./am2_dRes2; % 1&4
        sum12 = sum12 + 2.*exp_tm_NT.*(1-exp_1).^2.*(1-sgn.*exp_NT).*(sgn-exp_NT)./...
            (1+exp_1).^2./am2_dRes2; % 2&3
                % corresponds to the multiplication of the partial oscillations
        % from the first pulse with the incomplete oscillations from the
        % other pulses
        sum12 = sum12 - 2.*2.*exp_tm.*exp_delta_smalldel./am2_dRes2.*(exp_smalldel-exp_NT).*...
            (1-exp_smalldel_NT); % 1&3, 2&4
        sum12 = sum12 + 2.*exp_tm.*exp_delta.*exp_delta_smalldel./am2_dRes2.*(exp_smalldel-exp_NT).*...
            (1-exp_smalldel_NT); % 1&4
        sum12 = sum12 + 2.*exp_tm_smalldel./am2_dRes2.*(exp_smalldel-exp_NT).*...
            (1-exp_smalldel_NT); % 2&3
                % corresponds to the multiplication of the partial oscillations
        % with full oscillations   
         sum12 = sum12  - 2.*2.*sgn.*(1-exp_1)./(1+exp_1)./am2_dRes2.*(1-exp_smalldel_NT).*...
              (-exp_delta.*exp_tm.*exp_NT+sgn.*exp_delta_smalldel.*exp_tm.*exp_NT+2.* exp_NT-exp_delta_smalldel.*exp_tm-...
              2.*sgn+ sgn.*exp_delta.*exp_tm); % 1&3, 2&4
          
        sum12 = sum12  + 2.*sgn.*(1-exp_1)./(1+exp_1)./am2_dRes2.*(1-exp_smalldel_NT).*...
              (-exp_delta.*exp_delta.*exp_tm.*exp_NT+sgn.*exp_delta.*exp_delta_smalldel.*exp_tm.*exp_NT+2.* exp_NT-exp_delta_smalldel.*exp_delta.*exp_tm-...
              2.*sgn+ sgn.*exp_delta.*exp_delta.*exp_tm); % 1&4
        sum12 = sum12  + 2.*sgn.*(1-exp_1)./(1+exp_1)./am2_dRes2.*(1-exp_smalldel_NT).*...
              (-exp_tm.*exp_NT+sgn.*exp_tm_smalldel.*exp_NT+2.* exp_NT-exp_tm_smalldel-...
              2.*sgn+ sgn.*exp_tm); % 2&3  
          
          
        sumterms1 =(G1mx_rep.^2+G2mx_rep.^2).*sum11;   
        sumterms1 = sumterms1+(G1mx_rep.*G2mx_rep.*CosPsimx_rep).*sum12; 
        
        sumterms1 = Bmat_rep.*sumterms1;
        
        sumterms2 = (G1mx_par_rep.^2+G2mx_par_rep.^2).*sum11; 
        sumterms2 = sumterms2+ (G1mx_par_rep.*G2mx_par_rep).*sum12; 
        
        sumterms2 = Bmat_rep.*sumterms2;

       % planar geometry

        % corresponds to the multiplication of the full and partial oscillations with themselves
        sumterms2_plane = 4.*(G1mx_par_rep.^2+G2mx_par_rep.^2).*(-1 + exp_smalldel_NT_plane + NTmx_rep.*(-1+exp_1_plane) + am_dRes_plane.*smalldelmx_rep)./am2_dRes2_plane;

        % corresponds to the multiplication of the full oscillations with
        % other full oscillations from the same gradient
        sumterms2_plane = sumterms2_plane - 4.*(G1mx_par_rep.^2+G2mx_par_rep.^2).* sgn.*(1-exp_1_plane).^2./am2_dRes2_plane./(1+exp_1_plane).^2.*...
            (exp_NT_plane + sgn.*(NTmx_rep-1)+sgn.*NTmx_rep.*exp_1_plane);

        % corresponds to the multiplication of the full oscillations with
        % other full oscillations from different gradients
        sumterms2_plane = sumterms2_plane + 2.*(G1mx_par_rep.^2+G2mx_par_rep.^2).*exp_delta_NT_plane.*(1-exp_1_plane).^2.*(1-sgn.*exp_NT_plane).*(sgn-exp_NT_plane)./...
            (1+exp_1_plane).^2./am2_dRes2_plane; % 1&2, 3&4
        sumterms2_plane = sumterms2_plane - 2.*(2*G1mx_par_rep.*G2mx_par_rep).*exp_delta_tm_NT_plane.*(1-exp_1_plane).^2.*(1-sgn.*exp_NT_plane).*(sgn-exp_NT_plane)./...
            (1+exp_1_plane).^2./am2_dRes2_plane; % 1%3, 2&4
        sumterms2_plane = sumterms2_plane + 2.*G1mx_par_rep.*G2mx_par_rep.*exp_delta_tm_NT_plane.*exp_delta_plane.*(1-exp_1_plane).^2.*(1-sgn.*exp_NT_plane).*(sgn-exp_NT_plane)./...
            (1+exp_1_plane).^2./am2_dRes2_plane; % 1&4
        sumterms2_plane = sumterms2_plane + 2.*G1mx_par_rep.*G2mx_par_rep.*exp_tm_NT_plane.*(1-exp_1_plane).^2.*(1-sgn.*exp_NT_plane).*(sgn-exp_NT_plane)./...
            (1+exp_1_plane).^2./am2_dRes2_plane; % 2&3
        

        % corresponds to the multiplication of the partial oscillations
        % from the first pulse with the incomplete oscillations from the
        % other pulses
        sumterms2_plane = sumterms2_plane + 2.*(G1mx_par_rep.^2+G2mx_par_rep.^2).*exp_delta_smalldel_plane./am2_dRes2_plane.*(exp_smalldel_plane-exp_NT_plane).*...
            (1-exp_smalldel_NT_plane); % 1&2, 3&4
        sumterms2_plane = sumterms2_plane - 2.*(2*G1mx_par_rep.*G2mx_par_rep).*exp_tm_plane.*exp_delta_smalldel_plane./am2_dRes2_plane.*(exp_smalldel_plane-exp_NT_plane).*...
            (1-exp_smalldel_NT_plane); % 1&3, 2&4
        sumterms2_plane = sumterms2_plane + 2.*(G1mx_par_rep.*G2mx_par_rep).*exp_tm_plane.*exp_delta_plane.*exp_delta_smalldel_plane./am2_dRes2_plane.*(exp_smalldel_plane-exp_NT_plane).*...
            (1-exp_smalldel_NT_plane); % 1&4
        sumterms2_plane = sumterms2_plane + 2.*(G1mx_par_rep.*G2mx_par_rep).*exp_tm_smalldel_plane./am2_dRes2_plane.*(exp_smalldel_plane-exp_NT_plane).*...
            (1-exp_smalldel_NT_plane); % 2&3

        % corresponds to the multiplication of the partial oscillations
        % with full oscillations        
        
         sumterms2_plane = sumterms2_plane + 2.*(G1mx_par_rep.^2+G2mx_par_rep.^2).*sgn.*(1-exp_1_plane)./(1+exp_1_plane)./am2_dRes2_plane.*(1-exp_smalldel_NT_plane).*...
              (-exp_delta_plane.*exp_NT_plane+sgn.*exp_delta_smalldel_plane.*exp_NT_plane+2.* exp_NT_plane-exp_delta_smalldel_plane-...
              2.*sgn+ sgn.*exp_delta_plane); % 1&2, 3&4
         sumterms2_plane = sumterms2_plane  - 2.*(2*G1mx_par_rep.*G2mx_par_rep).*sgn.*(1-exp_1_plane)./(1+exp_1_plane)./am2_dRes2_plane.*(1-exp_smalldel_NT_plane).*...
              (-exp_delta_plane.*exp_tm_plane.*exp_NT_plane+sgn.*exp_delta_smalldel_plane.*exp_tm_plane.*exp_NT_plane+2.* exp_NT_plane-exp_delta_smalldel_plane.*exp_tm_plane-...
              2.*sgn+ sgn.*exp_delta_plane.*exp_tm_plane); % 1&3, 2&4
          
        sumterms2_plane = sumterms2_plane  + 2.*(G1mx_par_rep.*G2mx_par_rep).*sgn.*(1-exp_1_plane)./(1+exp_1_plane)./am2_dRes2_plane.*(1-exp_smalldel_NT_plane).*...
              (-exp_delta_plane.*exp_delta_plane.*exp_tm_plane.*exp_NT_plane+sgn.*exp_delta_plane.*exp_delta_smalldel_plane.*exp_tm_plane.*exp_NT_plane+2.* exp_NT_plane-exp_delta_smalldel_plane.*exp_delta_plane.*exp_tm_plane-...
              2.*sgn+ sgn.*exp_delta_plane.*exp_delta_plane.*exp_tm_plane); % 1&4
        sumterms2_plane = sumterms2_plane  + 2.*(G1mx_par_rep.*G2mx_par_rep).*sgn.*(1-exp_1_plane)./(1+exp_1_plane)./am2_dRes2_plane.*(1-exp_smalldel_NT_plane).*...
              (-exp_tm_plane.*exp_NT_plane+sgn.*exp_tm_smalldel_plane.*exp_NT_plane+2.* exp_NT_plane-exp_tm_smalldel_plane-...
              2.*sgn+ sgn.*exp_tm_plane); % 2&3  

        sumterms2_plane = Bmat_plane.*sumterms2_plane;

        testinds = find(sumterms2(:,:,end)>0);
        test = sumterms2(testinds,1)./sumterms2(testinds,end);
        if(min(abs(test))<1E2)
            warning('Ratio of largest to smallest terms in GPD model sum is <1E2.  May need more terms.');
            x;
        end

        s = sum(sumterms1,3)-sum(sumterms2,3);
        s = reshape(s,[l_q,l_a]);
        if(min(s)+1E-12<0)
            warning('Negative sums found in GPD sum.  Setting to zero.');
            s(find(s<0))=0;
            x;
        end

        s_plane = sum(sumterms2_plane,3);
        s_plane = reshape(s_plane,[l_q,l_a]);
       
        logE = -0.5.*GAMMA^2.*s;
        ePerp = exp(logE);
        logEpar = -0.5.*GAMMA^2.*s_plane;
        ePar = exp(logEpar);
        E_r_matrix = ePar.*ePerp;
        E_r=sum(E_r_matrix.*weight_matrix,2);
        E=E_r;

       
     else % the 2nd gradient is the negative of the first one reflected 
         
        % Parallel component - planar restriction

        %precomputing values that appear often in the expresions:
        % rearrange all formulas to have only negative exponentials
        exp_NT_plane = exp(-lambda_mat_plane.*dRes.*NTmx_rep./(2.*freqmx_rep));
        exp_1_plane = exp(-lambda_mat_plane.*dRes./(2.*freqmx_rep));
        exp_smalldel_NT_plane = exp(-lambda_mat_plane.*dRes.*(smalldelmx_rep-NTmx_rep./(2.*freqmx_rep)));
        exp_delta_NT_plane = exp(-lambda_mat_plane.*dRes.*(deltamx_rep-NTmx_rep./(2.*freqmx_rep)));
        exp_delta_tm_NT_plane = exp(-lambda_mat_plane.*dRes.*(deltamx_rep+tmmx_rep-NTmx_rep./(2.*freqmx_rep)));
        exp_tm_NT_plane = exp(-lambda_mat_plane.*dRes.*(tmmx_rep-NTmx_rep./(2.*freqmx_rep)));
        exp_delta_plane = exp(-lambda_mat_plane.*dRes.*deltamx_rep);
        exp_tm_plane = exp(-lambda_mat_plane.*dRes.*tmmx_rep);
        exp_tm_smalldel_plane = exp(-lambda_mat_plane.*dRes.*(tmmx_rep-smalldelmx_rep));
        exp_delta_smalldel_plane = exp(-lambda_mat_plane.*dRes.*(deltamx_rep-smalldelmx_rep));
        am_dRes_plane = lambda_mat_plane.*dRes;
        am2_dRes2_plane = lambda_mat_plane.^2.*dRes.^2;

        % Perpendicular component

        %precomputing values that appear often in the expresions:
        % rearrange all formulas to have only negative exponentials
        exp_NT = exp(-amSq.*dRes.*NTmx_rep./(2.*freqmx_rep));
        exp_1 = exp(-amSq.*dRes./(2.*freqmx_rep));
        exp_smalldel_NT = exp(-amSq.*dRes.*(smalldelmx_rep-NTmx_rep./(2.*freqmx_rep)));
        exp_delta_NT = exp(-amSq.*dRes.*(deltamx_rep-NTmx_rep./(2.*freqmx_rep)));
        exp_delta_tm_NT = exp(-amSq.*dRes.*(deltamx_rep+tmmx_rep-NTmx_rep./(2.*freqmx_rep)));
        exp_tm_NT = exp(-amSq.*dRes.*(tmmx_rep-NTmx_rep./(2.*freqmx_rep)));
        exp_delta = exp(-amSq.*dRes.*deltamx_rep);
        exp_tm = exp(-amSq.*dRes.*tmmx_rep);
        exp_tm_smalldel = exp(-amSq.*dRes.*(tmmx_rep-smalldelmx_rep));
        exp_delta_smalldel = exp(-amSq.*dRes.*(deltamx_rep-smalldelmx_rep));
        am_dRes = amSq.*dRes;
        am2_dRes2 = amSq.^2.*dRes.^2;
        sgn = (-1).^NTmx_rep;
        
        
       % corresponds to the multiplication of the full and incomplete oscillations with themselves
        sum11 = 4.*(-1 + exp_smalldel_NT + NTmx_rep.*(-1+exp_1) + am_dRes.*smalldelmx_rep)./am2_dRes2;
        
        % corresponds to the multiplication of the full oscillations with
        % other full oscillations from the same gradient
        sum11 = sum11 - 4.* sgn.*(1-exp_1).^2./am2_dRes2./(1+exp_1).^2.*...
            (exp_NT + sgn.*(NTmx_rep-1)+sgn.*NTmx_rep.*exp_1);
                % corresponds to the multiplication of the full oscillations with
        % other full oscillations from different gradients
        sum11 = sum11 + 2.*sgn.*exp_delta_NT.*exp_smalldel_NT.*(1-exp_1).^2.*(sgn.*exp_NT-1).*(sgn-exp_NT)./...
            (1+exp_1).^2./am2_dRes2; % 1&2, 3&4
                % corresponds to the multiplication of the partial oscillations
        % from the first pulse with the incomplete oscillations from the
        % other pulses
        sum11 = sum11 + 2.*exp_delta_smalldel./am2_dRes2.*(exp_smalldel_NT-1).*...
            (1-exp_smalldel_NT); % 1&2, 3&4
        % corresponds to the multiplication of the partial oscillations
        % with full oscillations
        sum11 = sum11 -2./am2_dRes2./(1+exp_1).*(1-exp_smalldel_NT).*(-1+exp_1).*(sgn.*exp_NT-1)+...
          2./am2_dRes2./(1+exp_1).*(-1+exp_1).*(sgn.*exp_NT-1).*(1-exp_smalldel_NT).*exp_delta_NT-...
          2.*sgn./am2_dRes2./(1+exp_1).*(1-exp_1).*(exp_smalldel_NT-1).*(-exp_delta+sgn.*(exp_delta_NT))+...
          2.*sgn./am2_dRes2./(1+exp_1).*(1-exp_1).*(exp_smalldel_NT-1).*(sgn-exp_NT); %1&2,3&4

        
        
        
        % corresponds to the multiplication of the full oscillations with
        % other full oscillations from different gradients
        sum12 =  - 2.*2.*sgn.*exp_delta_tm_NT.*exp_smalldel_NT.*(1-exp_1).^2.*(sgn.*exp_NT-1).*(sgn-exp_NT)./...
            (1+exp_1).^2./am2_dRes2; % 1&3, 2&4
        sum12 = sum12 + 2.*sgn.*exp_delta_tm_NT.*exp_delta.*exp_smalldel_NT.*(1-exp_1).^2.*(sgn.*exp_NT-1).*(sgn-exp_NT)./...
            (1+exp_1).^2./am2_dRes2; % 1&4, 
        sum12 = sum12 + 2.*sgn.*exp_tm_NT.*exp_smalldel_NT.*(1-exp_1).^2.*(sgn.*exp_NT-1).*(sgn-exp_NT)./...
            (1+exp_1).^2./am2_dRes2; % 2&3, 

        % corresponds to the multiplication of the partial oscillations
        % from the first pulse with the incomplete oscillations from the
        % other pulses
        sum12 = sum12 - 2.*2.*exp_delta_smalldel.*exp_tm./am2_dRes2.*(exp_smalldel_NT-1).*...
            (1-exp_smalldel_NT); % 1&3, 2&4
        sum12 = sum12 + 2.*exp_delta_smalldel.*exp_tm.*exp_delta./am2_dRes2.*(exp_smalldel_NT-1).*...
            (1-exp_smalldel_NT); % 1&4
        sum12 = sum12 + 2.*exp_tm_smalldel./am2_dRes2.*(exp_smalldel_NT-1).*...
            (1-exp_smalldel_NT); % 2&3

        % corresponds to the multiplication of the partial oscillations
        % with full oscillations
  
        sum12 = sum12 +2.*2./am2_dRes2./(1+exp_1).*(1-exp_smalldel_NT).*(-1+exp_1).*(sgn.*exp_NT-1)+...
      2./am2_dRes2./(1+exp_1).*(-1+exp_1).*(sgn.*exp_NT-1).*(1-exp_smalldel_NT).*exp_delta_NT.*exp_tm-...
      2.*sgn./am2_dRes2./(1+exp_1).*(1-exp_1).*(exp_smalldel_NT-1).*(-exp_delta.*exp_tm+sgn.*(exp_delta_NT.*exp_tm))+...
      2.*sgn./am2_dRes2./(1+exp_1).*(1-exp_1).*(exp_smalldel_NT-1).*(sgn-exp_NT); %1&2,3&4
  
        sum12 = sum12 -2./am2_dRes2./(1+exp_1).*(1-exp_smalldel_NT).*(-1+exp_1).*(sgn.*exp_NT-1)+...
      2./am2_dRes2./(1+exp_1).*(-1+exp_1).*(sgn.*exp_NT-1).*(1-exp_smalldel_NT).*exp_delta_NT.*exp_tm.*exp_delta-...
      2.*sgn./am2_dRes2./(1+exp_1).*(1-exp_1).*(exp_smalldel_NT-1).*(-exp_delta.*exp_delta.*exp_tm+sgn.*(exp_delta_NT.*exp_tm.*exp_delta))+...
      2.*sgn./am2_dRes2./(1+exp_1).*(1-exp_1).*(exp_smalldel_NT-1).*(sgn-exp_NT); %1&4
  
        sum12 = sum12 -2./am2_dRes2./(1+exp_1).*(1-exp_smalldel_NT).*(-1+exp_1).*(sgn.*exp_NT-1)+...
      2./am2_dRes2./(1+exp_1).*(-1+exp_1).*(sgn.*exp_NT-1).*(1-exp_smalldel_NT).*exp_tm_NT-...
      2.*sgn./am2_dRes2./(1+exp_1).*(1-exp_1).*(exp_smalldel_NT-1).*(-exp_tm+sgn.*(exp_tm_NT))+...
      2.*sgn./am2_dRes2./(1+exp_1).*(1-exp_1).*(exp_smalldel_NT-1).*(sgn-exp_NT); %2&3

  
        sumterms1 = (G1mx_rep.^2+G2mx_rep.^2).*sum11;
        sumterms1 = sumterms1 + (G1mx_rep.*G2mx_rep.*CosPsimx_rep).*sum12;
        sumterms1 = Bmat_rep.*sumterms1;
        
        sumterms2 = (G1mx_par_rep.^2+G2mx_par_rep.^2).*sum11; 
        sumterms2 = sumterms2+ (G1mx_par_rep.*G2mx_par_rep).*sum12;         
        sumterms2 = Bmat_rep.*sumterms2;     


        % planar restriction


          % corresponds to the multiplication of the full and incomplete oscillations with themselves
        sumterms2_plane = 4.*(G1mx_par_rep.^2+G2mx_par_rep.^2).*(-1 + exp_smalldel_NT_plane + NTmx_rep.*(-1+exp_1_plane) + am_dRes_plane.*smalldelmx_rep)./am2_dRes2_plane;
        
        % corresponds to the multiplication of the full oscillations with
        % other full oscillations from the same gradient
        sumterms2_plane = sumterms2_plane - 4.*(G1mx_par_rep.^2+G2mx_par_rep.^2).* sgn.*(1-exp_1_plane).^2./am2_dRes2_plane./(1+exp_1_plane).^2.*...
            (exp_NT_plane + sgn.*(NTmx_rep-1)+sgn.*NTmx_rep.*exp_1_plane);

        % corresponds to the multiplication of the full oscillations with
        % other full oscillations from different gradients
        sumterms2_plane = sumterms2_plane + 2.*(G1mx_par_rep.^2+G2mx_par_rep.^2).*sgn.*exp_delta_NT_plane.*exp_smalldel_NT_plane.*(1-exp_1_plane).^2.*(sgn.*exp_NT_plane-1).*(sgn-exp_NT_plane)./...
            (1+exp_1_plane).^2./am2_dRes2_plane; % 1&2, 3&4
        sumterms2_plane = sumterms2_plane - 2.*(2*G1mx_par_rep.*G2mx_par_rep).*sgn.*exp_delta_tm_NT_plane.*exp_smalldel_NT_plane.*(1-exp_1_plane).^2.*(sgn.*exp_NT_plane-1).*(sgn-exp_NT_plane)./...
            (1+exp_1_plane).^2./am2_dRes2_plane; % 1&3, 2&4
        sumterms2_plane = sumterms2_plane + 2.*(G1mx_par_rep.*G2mx_par_rep).*sgn.*exp_delta_tm_NT.*exp_delta.*exp_smalldel_NT_plane.*(1-exp_1_plane).^2.*(sgn.*exp_NT_plane-1).*(sgn-exp_NT_plane)./...
            (1+exp_1_plane).^2./am2_dRes2_plane; % 1&4, 
        sumterms2_plane = sumterms2_plane + 2.*(G1mx_par_rep.*G2mx_par_rep).*sgn.*exp_tm_NT_plane.*exp_smalldel_NT_plane.*(1-exp_1_plane).^2.*(sgn.*exp_NT_plane-1).*(sgn-exp_NT_plane)./...
            (1+exp_1_plane).^2./am2_dRes2_plane; % 2&3, 

        % corresponds to the multiplication of the partial oscillations
        % from the first pulse with the incomplete oscillations from the
        % other pulses
        sumterms2_plane = sumterms2_plane + 2.*(G1mx_par_rep.^2+G2mx_par_rep.^2).*exp_delta_smalldel_plane./am2_dRes2_plane.*(exp_smalldel_NT_plane-1).*...
            (1-exp_smalldel_NT_plane); % 1&2, 3&4
        sumterms2_plane = sumterms2_plane - 2.*(2*G1mx_par_rep.*G2mx_par_rep).*exp_delta_smalldel_plane.*exp_tm_plane./am2_dRes2_plane.*(exp_smalldel_NT_plane-1).*...
            (1-exp_smalldel_NT_plane); % 1&3, 2&4
        sumterms2_plane = sumterms2_plane + 2.*(G1mx_par_rep.*G2mx_par_rep).*exp_delta_smalldel_plane.*exp_tm_plane.*exp_delta_plane./am2_dRes2_plane.*(exp_smalldel_NT_plane-1).*...
            (1-exp_smalldel_NT_plane); % 1&4
        sumterms2_plane = sumterms2_plane + 2.*(G1mx_par_rep.*G2mx_par_rep).*exp_tm_smalldel_plane./am2_dRes2_plane.*(exp_smalldel_NT_plane-1).*...
            (1-exp_smalldel_NT_plane); % 2&3

        % corresponds to the multiplication of the partial oscillations
        % with full oscillations
        sumterms2_plane = sumterms2_plane -2.*(G1mx_par_rep.^2+G2mx_par_rep.^2)./am2_dRes2_plane./(1+exp_1_plane).*(1-exp_smalldel_NT_plane).*(-1+exp_1_plane).*(sgn.*exp_NT_plane-1)+...
      2./am2_dRes2_plane./(1+exp_1_plane).*(-1+exp_1_plane).*(sgn.*exp_NT_plane-1).*(1-exp_smalldel_NT_plane).*exp_delta_NT_plane-...
      2.*sgn./am2_dRes2_plane./(1+exp_1_plane).*(1-exp_1_plane).*(exp_smalldel_NT_plane-1).*(-exp_delta_plane+sgn.*(exp_delta_NT_plane))+...
      2.*sgn./am2_dRes2_plane./(1+exp_1_plane).*(1-exp_1_plane).*(exp_smalldel_NT_plane-1).*(sgn-exp_NT_plane); %1&2,3&4
  
        sumterms2_plane = sumterms2_plane +2.*(2*G1mx_par_rep.*G2mx_par_rep)./am2_dRes2_plane./(1+exp_1_plane).*(1-exp_smalldel_NT_plane).*(-1+exp_1_plane).*(sgn.*exp_NT_plane-1)+...
      2./am2_dRes2_plane./(1+exp_1_plane).*(-1+exp_1_plane).*(sgn.*exp_NT_plane-1).*(1-exp_smalldel_NT_plane).*exp_delta_NT_plane.*exp_tm_plane-...
      2.*sgn./am2_dRes2_plane./(1+exp_1_plane).*(1-exp_1_plane).*(exp_smalldel_NT_plane-1).*(-exp_delta_plane.*exp_tm_plane+sgn.*(exp_delta_NT_plane.*exp_tm_plane))+...
      2.*sgn./am2_dRes2_plane./(1+exp_1_plane).*(1-exp_1_plane).*(exp_smalldel_NT_plane-1).*(sgn-exp_NT_plane); %1&2,3&4
  
        sumterms2_plane = sumterms2_plane -2.*(G1mx_par_rep.*G2mx_par_rep)./am2_dRes2_plane./(1+exp_1_plane).*(1-exp_smalldel_NT_plane).*(-1+exp_1_plane).*(sgn.*exp_NT_plane-1)+...
      2./am2_dRes2_plane./(1+exp_1_plane).*(-1+exp_1_plane).*(sgn.*exp_NT_plane-1).*(1-exp_smalldel_NT_plane).*exp_delta_NT_plane.*exp_tm_plane.*exp_delta_plane-...
      2.*sgn./am2_dRes2_plane./(1+exp_1_plane).*(1-exp_1_plane).*(exp_smalldel_NT_plane-1).*(-exp_delta_plane.*exp_delta_plane.*exp_tm_plane+sgn.*(exp_delta_NT_plane.*exp_tm_plane.*exp_delta_plane))+...
      2.*sgn./am2_dRes2_plane./(1+exp_1_plane).*(1-exp_1_plane).*(exp_smalldel_NT_plane-1).*(sgn-exp_NT_plane); %1&4
  
  sumterms2_plane = sumterms2_plane -2.*(G1mx_par_rep.*G2mx_par_rep)./am2_dRes2_plane./(1+exp_1_plane).*(1-exp_smalldel_NT_plane).*(-1+exp_1_plane).*(sgn.*exp_NT_plane-1)+...
      2./am2_dRes2_plane./(1+exp_1_plane).*(-1+exp_1_plane).*(sgn.*exp_NT_plane-1).*(1-exp_smalldel_NT_plane).*exp_tm_NT_plane-...
      2.*sgn./am2_dRes2_plane./(1+exp_1_plane).*(1-exp_1_plane).*(exp_smalldel_NT_plane-1).*(-exp_tm_plane+sgn.*(exp_tm_NT_plane))+...
      2.*sgn./am2_dRes2_plane./(1+exp_1_plane).*(1-exp_1_plane).*(exp_smalldel_NT_plane-1).*(sgn-exp_NT_plane); %2&3

        sumterms2_plane = Bmat_plane.*sumterms2_plane;


        testinds = find(sumterms2(:,:,end)>0);
        test = sumterms2(testinds,1)./sumterms2(testinds,end);

        if(min(test)<1E4)
            warning('Ratio of largest to smallest terms in GPD model sum is <1E4.  May need more terms.');
            x;
        end

        s =  sum(sumterms1,3)-sum(sumterms2,3);
        s = reshape(s,[l_q,l_a]);

     
        s_plane = sum(sumterms2_plane,3);
        s_plane = reshape(s_plane,[l_q,l_a]);
        if(min(s)+1E-12<0)
            warning('Negative sums found in GPD sum.  Setting to zero.');
            s(find(s<0))=0;
            x;
        end
        logE = -0.5.*GAMMA^2.*s;
        ePerp = exp(logE);
        logEpar = -0.5.*GAMMA^2.*s_plane;
        ePar = exp(logEpar);
        E_r_matrix = ePar.*ePerp; 
        E_r=sum(E_r_matrix.*weight_matrix,2);
        E=E_r;
     end    
else  
    error('Not implemented yet');
    phase = protocol.phase';
    phase = mod(phase,2*pi);
    phase(phase>pi) = phase(phase>pi)-2*pi;
    phase(phase<0) = pi-abs(phase(phase<0));
    phdelay = phase ./(2 *pi()* niu);
  
    phdelaymx = repmat(phdelay,[1,l_a]);
    phdelaymx_rep =phdelaymx(:);
    phdelaymx_rep = repmat(phdelaymx_rep,[1 1 k_max]);

    % the number of integer half periods:
    NT = floor(2.*(smalldel-phdelay).*niu+0.00000000001); % needed for numerical precision (if NT is int, the floor rounds it down)
    NTmx_rep = floor(2.*(smalldelmx_rep-phdelaymx_rep).*freqmx_rep+0.00000000001);
    
    if ~isfield(protocol,'mirror') || protocol.mirror == 0 % the 2nd gradient is the negative of the first one
        sgn = (-1).^NT;
        % calculate b value
        bval = GAMMA.^2.*(G1mx_par.^2+G2mx_par.^2).*(1./3.*(smalldel-NT./2./niu-phdelay).^3 + (delta-smalldel).*...
        (sgn.*(smalldel- (0.5.*(1-sgn)+NT)./2./niu-phdelay)-phdelay).^2  +(phdelay.^3)./3 +...
        sgn.*((-1+sgn-2.*NT+4.*niu.*(smalldel-phdelay)).^3)./192./(niu.^3) +...   
         NT./96./(niu.^3).*(8+12.*NT.*(1+NT)- 24.*smalldel.*niu.*(1 + 2.*NT) +48.*smalldel.^2.*niu.^2+...
        48.*NT.*niu.*phdelay - 96.*niu.^2.*phdelay.*(smalldel-phdelay)) +...
        1./3.*(phdelay-sgn.*(-1+sgn-2.*NT+4.*niu.*(smalldel-phdelay))./4./niu).^3 +...     
        1/3.*sgn.*((-1+sgn+4.*niu.*phdelay).^3./64./niu.^3-(phdelay-sgn.*(-1+sgn-2.*NT+4.*niu.*(smalldel-phdelay))./4./niu).^3));
    
        bval_matrix=repmat(bval,[1,l_a]);

        % Parallel component
        logEPar=-bval_matrix.*dRes;
        ePar=exp(logEPar);
        
        % Perpendicular component

        %precomputing values that appear often in the expresions:
        % rearrange all formulas to have only negative exponentials
        exp_NT = exp(-amSq.*dRes.*NTmx_rep./(2.*freqmx_rep));
        exp_1 = exp(-amSq.*dRes./(2.*freqmx_rep));
        exp_smalldel_NT_phdelay = exp(-amSq.*dRes.*(smalldelmx_rep-NTmx_rep./(2.*freqmx_rep)-...
            phdelaymx_rep));
        exp_delta_NT_phdelay = exp(-amSq.*dRes.*(deltamx_rep-NTmx_rep./(2.*freqmx_rep)-...
            phdelaymx_rep));
        exp_tm_NT_phdelay = exp(-amSq.*dRes.*(tmmx_rep-NTmx_rep./(2.*freqmx_rep)-...
            phdelaymx_rep));
        exp_tm_phdelay = exp(-amSq.*dRes.*(tmmx_rep-phdelaymx_rep));
        exp_phdelay = exp(-amSq.*dRes.*phdelaymx_rep);
        exp_delta_phdelay = exp(-amSq.*dRes.*(deltamx_rep-phdelaymx_rep));
        exp_delta_NT = exp(-amSq.*dRes.*(deltamx_rep-NTmx_rep./(2.*freqmx_rep)));       
        exp_tm_NT = exp(-amSq.*dRes.*(tmmx_rep-NTmx_rep./(2.*freqmx_rep)));
        exp_delta = exp(-amSq.*dRes.*deltamx_rep);
        exp_tm = exp(-amSq.*dRes.*tmmx_rep);
        exp_tm_smalldel = exp(-amSq.*dRes.*(tmmx_rep-smalldelmx_rep));
        exp_delta_smalldel = exp(-amSq.*dRes.*(deltamx_rep-smalldelmx_rep));
        exp_smalldel = exp(-amSq.*dRes.*smalldelmx_rep);
        am_dRes = amSq.*dRes;
        am2_dRes2 = amSq.^2.*dRes.^2;

        sgn = (-1).^NTmx_rep;
        
       % corresponds to the multiplication of the full and incmoplete oscillations with themselves
        sum11 = 4.*(-1 + exp_smalldel_NT_phdelay + NTmx_rep.*(-1+exp_1) + am_dRes.*(smalldelmx_rep-phdelaymx_rep))./am2_dRes2;

        % corresponds to the multiplication of the partial oscillations due to non zero phase with themselves
        sum11 = sum11 + 4.* (am_dRes.*phdelaymx_rep+exp_phdelay-1)./am2_dRes2;

                % corresponds to the multiplication of the partial oscillations due
        % to non zero phase with the partial oscillations in the end of the
        % pulse
        sum11 = sum11 -( 2.*sgn.*(2.*exp_NT-exp_delta.*exp_NT-exp_delta_smalldel).*(1-exp_phdelay).*...
            (1-exp_smalldel_NT_phdelay)./am2_dRes2); % 1&2, 3&4
           % corresponds to the multiplication of the partial oscillations due
        % to non zero phase with the partial oscillations in the end of the
        % pulse
        sum11 = sum11 -( 2.*sgn.*(2.*exp_NT-exp_delta.*exp_NT-exp_delta_smalldel).*(1-exp_phdelay).*...
            (1-exp_smalldel_NT_phdelay)./am2_dRes2); % 1&2, 3&4
          % corresponds to the multiplication of the partial oscillations due
        % to non zero phase from the two gradients
        sum11 = sum11 - 2.*exp_delta_phdelay.*(1-exp_phdelay).^2./am2_dRes2; %1&2,3&4
        
        % corresponds to the multiplication of the full oscillations with
        % other full oscillations from the same gradient
        sum11 = sum11 - 4.* sgn.*(1-exp_1).^2./am2_dRes2./(1+exp_1).^2.*...
            (exp_NT + sgn.*(NTmx_rep-1)+sgn.*NTmx_rep.*exp_1);
        
        % corresponds to the multiplication of the partial oscillations with themselves
        sum11 = sum11 + 2.*exp_delta_smalldel./am2_dRes2.*(exp_smalldel-exp_NT.* exp_phdelay).*...
            (1-exp_smalldel_NT_phdelay); % 1&2, 3&4
        % corresponds to the multiplication of the partial oscillations
        % with full oscillations
          sum11 = sum11 - (2.*sgn.*(1-exp_1)./(1+exp_1)./am2_dRes2.*(1-exp_smalldel_NT_phdelay).*...
             (-sgn.*exp_delta_smalldel.*exp_phdelay.*exp_NT+exp_delta_smalldel.*exp_phdelay+exp_delta.*exp_NT-...
             2.*exp_NT-sgn.*exp_delta+2.*sgn));         % 1&2, 3&4
        
        
        
        % corresponds to the multiplication of the partial oscillations due
        % to non zero phase with full oscillations from different pulsese
       
        sum12 = -2*(2.*(-1+exp_1).*(sgn.*exp_NT-1).*(1-exp_phdelay)./am2_dRes2./(exp_1+1).*(exp_delta.*exp_tm-2) - ...
            2.*exp_delta_NT_phdelay.*exp_tm.*(1-exp_1).*(sgn-exp_NT).*(1-exp_phdelay)./am2_dRes2./(exp_1+1)); % 1&3, 2&4
        
        sum12 = sum12 +(2.*(-1+exp_1).*(sgn.*exp_NT-1).*(1-exp_phdelay)./am2_dRes2./(exp_1+1).*(exp_delta.*exp_delta.*exp_tm-2) - ...
            2.*exp_delta_NT_phdelay.*exp_tm.*exp_delta.*(1-exp_1).*(sgn-exp_NT).*(1-exp_phdelay)./am2_dRes2./(exp_1+1));%1&4

        sum12 = sum12 +(2.*(-1+exp_1).*(sgn.*exp_NT-1).*(1-exp_phdelay)./am2_dRes2./(exp_1+1).*(exp_tm-2) - ...
            2.*exp_tm_NT_phdelay.*(1-exp_1).*(sgn-exp_NT).*(1-exp_phdelay)./am2_dRes2./(exp_1+1));%1&4
        
        % corresponds to the multiplication of the partial oscillations due
        % to non zero phase with the partial oscillations in the end of the
        % pulse
        
        sum12 = sum12 + 2.*(2.*sgn.*(2.*exp_NT-exp_delta.*exp_tm.*exp_NT-exp_delta_smalldel.*exp_tm).*(1-exp_phdelay).*...
            (1-exp_smalldel_NT_phdelay)./am2_dRes2); % 1&3, 2&4
        
        sum12 = sum12 - (2.*sgn.*(2.*exp_NT-exp_delta.*exp_delta.*exp_tm.*exp_NT-exp_delta_smalldel.*exp_tm.*exp_delta).*(1-exp_phdelay).*...
            (1-exp_smalldel_NT_phdelay)./am2_dRes2); % 1&4
        
        sum12 = sum12 - (2.*sgn.*(2.*exp_NT-exp_tm.*exp_NT-exp_tm_smalldel).*(1-exp_phdelay).*...
            (1-exp_smalldel_NT_phdelay)./am2_dRes2); % 2&3

        % corresponds to the multiplication of the partial oscillations due
        % to non zero phase from the two gradients
   
        sum12 = sum12 + 2.*2.*exp_delta_phdelay.*exp_tm.*(1-exp_phdelay).^2./am2_dRes2; %1&3,2&4
        sum12 = sum12 - 2.*exp_delta_phdelay.*exp_tm.*exp_delta.*(1-exp_phdelay).^2./am2_dRes2; %1&4
        sum12 = sum12 - 2.*exp_tm_phdelay.*(1-exp_phdelay).^2./am2_dRes2; %2&3



        % corresponds to the multiplication of the full oscillations with
        % other full oscillations from different gradients 

        sum12 = sum12 - 2.*2.*exp_delta_NT.*exp_tm.*(1-exp_1).^2.*(1-sgn.*exp_NT).*(sgn-exp_NT)./...
            (1+exp_1).^2./am2_dRes2; %1&3, 2&4
        sum12 = sum12 + 2.*exp_delta_NT.*exp_tm.*exp_delta.*(1-exp_1).^2.*(1-sgn.*exp_NT).*(sgn-exp_NT)./...
            (1+exp_1).^2./am2_dRes2; %1&4
        sum12 = sum12 + 2.*exp_tm_NT.*(1-exp_1).^2.*(1-sgn.*exp_NT).*(sgn-exp_NT)./...
            (1+exp_1).^2./am2_dRes2; %2&3
      
         % corresponds to the multiplication of the partial oscillations with themselves

        sum12 = sum12 - 2.*2.*exp_delta_smalldel.*exp_tm./am2_dRes2.*(exp_smalldel-exp_NT.* exp_phdelay).*...
            (1-exp_smalldel_NT_phdelay); % 1&3, 2&4
        sum12 = sum12 + 2.*exp_delta_smalldel.*exp_tm.*exp_delta./am2_dRes2.*(exp_smalldel-exp_NT.* exp_phdelay).*...
            (1-exp_smalldel_NT_phdelay); % 1&4
        sum12 = sum12 + 2.*exp_tm_smalldel./am2_dRes2.*(exp_smalldel-exp_NT.* exp_phdelay).*...
            (1-exp_smalldel_NT_phdelay); % 2&3
     
         % corresponds to the multiplication of the partial oscillations
        % with full oscillations
         
         sum12 = sum12 + 2*(2.*sgn.*(1-exp_1)./(1+exp_1)./am2_dRes2.*(1-exp_smalldel_NT_phdelay).*...
             (-sgn.*exp_delta_smalldel.*exp_tm.*exp_phdelay.*exp_NT+exp_delta_smalldel.*exp_tm.*exp_phdelay+exp_delta.*exp_NT.*exp_tm-...
             2.*exp_NT-sgn.*exp_delta.*exp_tm+2.*sgn)); % 1&3, 2&4
         
         sum12 = sum12 - (2.*sgn.*(1-exp_1)./(1+exp_1)./am2_dRes2.*(1-exp_smalldel_NT_phdelay).*...
             (-sgn.*exp_delta_smalldel.*exp_tm.*exp_delta.*exp_phdelay.*exp_NT+exp_delta_smalldel.*exp_tm.*exp_delta.*exp_phdelay+exp_delta.*exp_NT.*exp_tm.*exp_delta-...
             2.*exp_NT-sgn.*exp_delta.*exp_delta+2.*sgn)); % 1&4
         
         sum12 = sum12 - (2.*sgn.*(1-exp_1)./(1+exp_1)./am2_dRes2.*(1-exp_smalldel_NT_phdelay).*...
             (-sgn.*exp_tm_smalldel.*exp_phdelay.*exp_NT+exp_tm_smalldel.*exp_phdelay+exp_tm.*exp_NT-...
             2.*exp_NT-sgn.*exp_tm+2.*sgn)); % 2&3

        sumterms1 = (G1mx_rep.^2+G2mx_rep.^2).*sum11;
        sumterms1 = sumterms1 + (G1mx_rep.*G2mx_rep.*CosPsimx_rep).*sum12;
        sumterms1 = Bmat_rep.*sumterms1;
        
        sumterms2 = (G1mx_par_rep.^2+G2mx_par_rep.^2).*sum11; 
        sumterms2 = sumterms2+ (G1mx_par_rep.*G2mx_par_rep).*sum12;         
        sumterms2 = Bmat_rep.*sumterms2;  
        
        testinds = find(sumterms2(:,:,end)>0);
        test = sumterms2(testinds,1)./sumterms2(testinds,end);
        if(min(test)<1E4)
            warning('Ratio of largest to smallest terms in GPD model sum is <1E4.  May need more terms.');
            x;
        end

        s = sum(sumterms1,3)-sum(sumterms2,3);


        s = reshape(s,[l_q,l_a]);

        if(min(s)<0)
            warning('Negative sums found in GPD sum.  Setting to zero.');
            s(find(s<0))=0;
            x;
        end

        logE = -0.5.*GAMMA^2.*s;
        ePerp = exp(logE);

        E_r_matrix = ePar.*ePerp;

        E_r=sum(E_r_matrix.*weight_matrix,2);

        E=E_r;

    else % the 2nd gradient is the negative of the first one reflected 
        sgn = (-1).^NT;
        bval = (G1mx_par.^2+G2mx_par.^2).*GAMMA.^2.*((delta-smalldel).*(sgn.*(smalldel-(0.5.*(1-sgn)+NT)./(2.*niu)-phdelay)-phdelay).^2 +...
            2.*phdelay.^3./3+ NT.*(1-6.*niu.*phdelay+12.*niu.^2.*phdelay.^2)./12./niu.^3+...
            sgn./3.*((phdelay-sgn./4./niu.*(sgn-1)).^3-2.*(phdelay-sgn./4./niu.*(-1+sgn-2.*NT+4.*niu.*(smalldel-phdelay))).^3)+...
            sgn./3.*((sgn-1+4.*niu.*phdelay).^3./64./niu.^3));
        bval_matrix=repmat(bval,[1,l_a]);

        % Parallel component
        logEPar=-bval_matrix.*dRes;
        ePar=exp(logEPar);

        % Perpendicular component

        %precomputing values that appear often in the expresions:
        % rearrange all formulas to have only negative exponentials
        exp_NT = exp(-amSq.*dRes.*NTmx_rep./(2.*freqmx_rep));
        exp_1 = exp(-amSq.*dRes./(2.*freqmx_rep));
        exp_smalldel_NT_phdelay = exp(-amSq.*dRes.*(smalldelmx_rep-NTmx_rep./(2.*freqmx_rep)-...
            phdelaymx_rep));
        exp_delta_NT_phdelay = exp(-amSq.*dRes.*(deltamx_rep-NTmx_rep./(2.*freqmx_rep)-...
            phdelaymx_rep));
        exp_tm_NT_phdelay = exp(-amSq.*dRes.*(tmmx_rep-NTmx_rep./(2.*freqmx_rep)-...
            phdelaymx_rep));
        exp_phdelay = exp(-amSq.*dRes.*phdelaymx_rep);
        exp_delta_phdelay = exp(-amSq.*dRes.*(deltamx_rep-phdelaymx_rep));
        exp_tm_phdelay = exp(-amSq.*dRes.*(tmmx_rep-phdelaymx_rep));
        exp_delta = exp(-amSq.*dRes.*deltamx_rep);
        exp_tm = exp(-amSq.*dRes.*tmmx_rep);
        exp_tm_smalldel = exp(-amSq.*dRes.*(tmmx_rep-smalldelmx_rep));    
       
        exp_smalldel_phdelay = exp(-amSq.*dRes.*(smalldelmx_rep-phdelaymx_rep));
        exp_delta_smalldel = exp(-amSq.*dRes.*(deltamx_rep-smalldelmx_rep));       
        am_dRes = amSq.*dRes;
        am2_dRes2 = amSq.^2.*dRes.^2;

        sgn = (-1).^NTmx_rep;
        
            % corresponds to the multiplication of the full oscillations with themselves
        sum11 = 4.*(G1mx_rep.^2+G2mx_rep.^2).*(-1 + exp_smalldel_NT_phdelay + NTmx_rep.*(-1+exp_1) + am_dRes.*(smalldelmx_rep-phdelaymx_rep))./am2_dRes2;

        % corresponds to the multiplication of the partial oscillations due to non zero phase with themselves
        sum11 = sum11 + 4.* (G1mx_rep.^2+G2mx_rep.^2).*(am_dRes.*phdelaymx_rep+exp_phdelay-1)./am2_dRes2;

         % corresponds to the multiplication of the partial oscillations due
        % to non zero phase with full oscillations
        sum11 = sum11 -(G1mx_rep.^2+G2mx_rep.^2).*(2.*(exp_1-1).*(sgn.*exp_NT-1).*(1-exp_phdelay)./am2_dRes2./(1+exp_1)-...
            2.*exp_delta_phdelay.*exp_smalldel_NT_phdelay.*(1-exp_1).*(sgn-exp_NT).*(exp_phdelay-1)./am2_dRes2./(1+exp_1) +...
            2.*sgn.*exp_delta_phdelay.*exp_smalldel_NT_phdelay.*(exp_1-1).*(sgn.*exp_NT-1).*(1-exp_phdelay)./am2_dRes2./(1+exp_1) -...
            2.*sgn.*(1-exp_1).*(sgn-exp_NT).*(1-exp_phdelay)./am2_dRes2./(1+exp_1)); % 1&2, 3&4
        
                 % corresponds to the multiplication of the partial oscillations due
        % to non zero phase with the partial oscillations in the end of the
        % pulse
        sum11 = sum11 + 4.*(G1mx_rep.^2+G2mx_rep.^2).*sgn./am2_dRes2.*(1-exp_phdelay).*(exp_delta_phdelay+exp_smalldel_phdelay-...
            exp_delta_NT_phdelay.*exp_smalldel_phdelay-exp_NT); %1&2, 3&4
        
           % corresponds to the multiplication of the partial oscillations due
        % to non zero phase from the two gradients
        sum11 = sum11 - 2.*(G1mx_rep.^2+G2mx_rep.^2).*exp_delta_phdelay.*exp_smalldel_phdelay.*(1-exp_phdelay).^2./am2_dRes2; % 1&2, 3&4
        
           % corresponds to the multiplication of the full oscillations with
        % other full oscillations from the same gradient 
        sum11 = sum11 - 4.*(G1mx_rep.^2+G2mx_rep.^2).* sgn.*(1-exp_1).^2./am2_dRes2./(1+exp_1).^2.*...
            (exp_NT + sgn.*(NTmx_rep-1)+sgn.*NTmx_rep.*exp_1);
        
          % corresponds to the multiplication of the full oscillations with
        % other full oscillations from different gradients  
        sum11 = sum11 + 2.*(G1mx_rep.^2+G2mx_rep.^2).*sgn.*exp_delta_NT_phdelay.*exp_smalldel_NT_phdelay.*(1-exp_1).^2.*...
            (sgn.*exp_NT-1).*(sgn-exp_NT)./am2_dRes2./(1+exp_1).^2; %1&2, 3&4
        
         % corresponds to the multiplication of the partial oscillations with themselves
        sum11 = sum11 + 2.*(G1mx_rep.^2+G2mx_rep.^2).*exp_delta_smalldel./am2_dRes2.*(exp_smalldel_NT_phdelay-1).*...
            (1-exp_smalldel_NT_phdelay); % 1&2, 3&4
        
            % corresponds to the multiplication of the partial oscillations
        % with full oscillations
          sum11 = sum11 + (G1mx_rep.^2+G2mx_rep.^2).*(2.*sgn.*(1-exp_1)./am2_dRes2./(1+exp_1).*(-exp_delta_NT_phdelay.*(sgn-exp_NT).*(-1+exp_smalldel_NT_phdelay) +...
              sgn.*(sgn.*exp_NT-1).*(1-exp_smalldel_NT_phdelay)+ ...
              (sgn-exp_NT).*(exp_smalldel_NT_phdelay-1)-...
              sgn.*exp_delta_NT_phdelay.*(sgn.*exp_NT-1).*(1-exp_smalldel_NT_phdelay))); %1&2, 3&4
        
    
          % corresponds to the multiplication of the partial oscillations due
        % to non zero phase with full oscillations
        
        sum12=  +(2*G1mx_rep.*G2mx_rep.*CosPsimx_rep).*(2.*(exp_1-1).*(sgn.*exp_NT-1).*(1-exp_phdelay)./am2_dRes2./(1+exp_1)-...
            2.*exp_delta_phdelay.*exp_tm.*exp_smalldel_NT_phdelay.*(1-exp_1).*(sgn-exp_NT).*(exp_phdelay-1)./am2_dRes2./(1+exp_1) +...
            2.*sgn.*exp_delta_phdelay.*exp_tm.*exp_smalldel_NT_phdelay.*(exp_1-1).*(sgn.*exp_NT-1).*(1-exp_phdelay)./am2_dRes2./(1+exp_1) -...
            2.*sgn.*(1-exp_1).*(sgn-exp_NT).*(1-exp_phdelay)./am2_dRes2./(1+exp_1)); % 1&3, 2&4
        
        sum12 = sum12 -(G1mx_rep.*G2mx_rep.*CosPsimx_rep).*(2.*(exp_1-1).*(sgn.*exp_NT-1).*(1-exp_phdelay)./am2_dRes2./(1+exp_1)-...
            2.*exp_delta_phdelay.*exp_delta.*exp_tm.*exp_smalldel_NT_phdelay.*(1-exp_1).*(sgn-exp_NT).*(exp_phdelay-1)./am2_dRes2./(1+exp_1) +...
            2.*sgn.*exp_delta_phdelay.*exp_delta.*exp_tm.*exp_smalldel_NT_phdelay.*(exp_1-1).*(sgn.*exp_NT-1).*(1-exp_phdelay)./am2_dRes2./(1+exp_1) -...
            2.*sgn.*(1-exp_1).*(sgn-exp_NT).*(1-exp_phdelay)./am2_dRes2./(1+exp_1)); % 1&4
        
        sum12 = sum12 -(G1mx_rep.*G2mx_rep.*CosPsimx_rep).*(2.*(exp_1-1).*(sgn.*exp_NT-1).*(1-exp_phdelay)./am2_dRes2./(1+exp_1)-...
            2.*exp_tm_phdelay.*exp_smalldel_NT_phdelay.*(1-exp_1).*(sgn-exp_NT).*(exp_phdelay-1)./am2_dRes2./(1+exp_1) +...
            2.*sgn.*exp_tm_phdelay.*exp_smalldel_NT_phdelay.*(exp_1-1).*(sgn.*exp_NT-1).*(1-exp_phdelay)./am2_dRes2./(1+exp_1) -...
            2.*sgn.*(1-exp_1).*(sgn-exp_NT).*(1-exp_phdelay)./am2_dRes2./(1+exp_1)); % 2&3

         % corresponds to the multiplication of the partial oscillations due
        % to non zero phase with the partial oscillations in the end of the
        % pulse
        sum12 = sum12 - 4.*(2*G1mx_rep.*G2mx_rep.*CosPsimx_rep).*sgn./am2_dRes2.*(1-exp_phdelay).*(exp_delta_phdelay.*exp_tm+exp_smalldel_phdelay-...
            exp_delta_NT_phdelay.*exp_tm.*exp_smalldel_phdelay-exp_NT); % 1&3 2&4
        sum12 = sum12 + 4.*(G1mx_rep.*G2mx_rep.*CosPsimx_rep).*sgn./am2_dRes2.*(1-exp_phdelay).*(exp_delta_phdelay.*exp_delta.*exp_tm+exp_smalldel_phdelay-...
            exp_delta_NT_phdelay.*exp_delta.*exp_tm.*exp_smalldel_phdelay-exp_NT); % 1&4
        sum12 = sum12 + 4.*(G1mx_rep.*G2mx_rep.*CosPsimx_rep).*sgn./am2_dRes2.*(1-exp_phdelay).*(exp_tm_phdelay+exp_smalldel_phdelay-...
            exp_tm_NT_phdelay.*exp_smalldel_phdelay-exp_NT); %2&3

        % corresponds to the multiplication of the partial oscillations due
        % to non zero phase from the two gradients
    
        sum12 = sum12 + 2.*(2*G1mx_rep.*G2mx_rep.*CosPsimx_rep).*exp_delta_phdelay.*exp_tm.*exp_smalldel_phdelay.*(1-exp_phdelay).^2./am2_dRes2; % 1&3, 2&4 
        sum12 = sum12 - 2.*(G1mx_rep.*G2mx_rep.*CosPsimx_rep).*exp_delta_phdelay.*exp_tm.*exp_delta.*exp_smalldel_phdelay.*(1-exp_phdelay).^2./am2_dRes2; % 1&4
        sum12 = sum12 - 2.*(G1mx_rep.*G2mx_rep.*CosPsimx_rep).*exp_tm_phdelay.*exp_smalldel_phdelay.*(1-exp_phdelay).^2./am2_dRes2; % 2&3

     


        % corresponds to the multiplication of the full oscillations with
        % other full oscillations from different gradients  
        sum12 = sum12 - 2.*(2*G1mx_rep.*G2mx_rep.*CosPsimx_rep).*sgn.*exp_delta_NT_phdelay.*exp_tm.*exp_smalldel_NT_phdelay.*(1-exp_1).^2.*...
            (sgn.*exp_NT-1).*(sgn-exp_NT)./am2_dRes2./(1+exp_1).^2; % 1&3 2&4
        sum12 = sum12 + 2.*(G1mx_rep.*G2mx_rep.*CosPsimx_rep).*sgn.*exp_delta_NT_phdelay.*exp_delta.*exp_tm.*exp_smalldel_NT_phdelay.*(1-exp_1).^2.*...
            (sgn.*exp_NT-1).*(sgn-exp_NT)./am2_dRes2./(1+exp_1).^2; % 1&4
        sum12 = sum12 + 2.*(G1mx_rep.*G2mx_rep.*CosPsimx_rep).*sgn.*exp_tm_NT_phdelay.*exp_smalldel_NT_phdelay.*(1-exp_1).^2.*...
            (sgn.*exp_NT-1).*(sgn-exp_NT)./am2_dRes2./(1+exp_1).^2; % 2&3


         % corresponds to the multiplication of the partial oscillations with themselves

        sum12 = sum12 - 2.*(2*G1mx_rep.*G2mx_rep.*CosPsimx_rep).*exp_delta_smalldel.*exp_tm./am2_dRes2.*(exp_smalldel_NT_phdelay-1).*...
            (1-exp_smalldel_NT_phdelay); % 1&3, 2&4
        sum12 = sum12 + 2.*(G1mx_rep.*G2mx_rep.*CosPsimx_rep).*exp_delta_smalldel.*exp_tm.*exp_delta./am2_dRes2.*(exp_smalldel_NT_phdelay-1).*...
            (1-exp_smalldel_NT_phdelay); % 1&4
        sum12 = sum12 + 2.*(G1mx_rep.*G2mx_rep.*CosPsimx_rep).*exp_tm_smalldel./am2_dRes2.*(exp_smalldel_NT_phdelay-1).*...
            (1-exp_smalldel_NT_phdelay); % 2&3

        % corresponds to the multiplication of the partial oscillations
        % with full oscillations
      
          sum12 = sum12 - (2*G1mx_rep.*G2mx_rep.*CosPsimx_rep).*(2.*sgn.*(1-exp_1)./am2_dRes2./(1+exp_1).*(-exp_delta_NT_phdelay.*exp_tm.*(sgn-exp_NT).*(-1+exp_smalldel_NT_phdelay) +...
              sgn.*(sgn.*exp_NT-1).*(1-exp_smalldel_NT_phdelay)+ ...
              (sgn-exp_NT).*(exp_smalldel_NT_phdelay-1)-...
              sgn.*exp_delta_NT_phdelay.*exp_tm.*(sgn.*exp_NT-1).*(1-exp_smalldel_NT_phdelay))); % 1&3, 2&4
          sum12 = sum12 + (G1mx_rep.*G2mx_rep.*CosPsimx_rep).*(2.*sgn.*(1-exp_1)./am2_dRes2./(1+exp_1).*(-exp_delta_NT_phdelay.*exp_tm.*exp_delta.*(sgn-exp_NT).*(-1+exp_smalldel_NT_phdelay) +...
              sgn.*(sgn.*exp_NT-1).*(1-exp_smalldel_NT_phdelay)+ ...
              (sgn-exp_NT).*(exp_smalldel_NT_phdelay-1)-...
              sgn.*exp_delta_NT_phdelay.*exp_tm.*exp_delta.*(sgn.*exp_NT-1).*(1-exp_smalldel_NT_phdelay)));%1 &4
          sum12 = sum12 + (G1mx_rep.*G2mx_rep.*CosPsimx_rep).*(2.*sgn.*(1-exp_1)./am2_dRes2./(1+exp_1).*(-exp_tm_NT_phdelay.*(sgn-exp_NT).*(-1+exp_smalldel_NT_phdelay) +...
              sgn.*(sgn.*exp_NT-1).*(1-exp_smalldel_NT_phdelay)+ ...
              (sgn-exp_NT).*(exp_smalldel_NT_phdelay-1)-...
              sgn.*exp_tm_NT_phdelay.*(sgn.*exp_NT-1).*(1-exp_smalldel_NT_phdelay))); % 2&3

        sumterms1 = (G1mx_rep.^2+G2mx_rep.^2).*sum11;
        sumterms1 = sumterms1 + (G1mx_rep.*G2mx_rep.*CosPsimx_rep).*sum12;
        sumterms1 = Bmat_rep.*sumterms1;
        
        sumterms2 = (G1mx_par_rep.^2+G2mx_par_rep.^2).*sum11; 
        sumterms2 = sumterms2+ (G1mx_par_rep.*G2mx_par_rep).*sum12;         
        sumterms2 = Bmat_rep.*sumterms2;    

        testinds = find(sumterms2(:,:,end)>0);
        test = sumterms2(testinds,1)./sumterms2(testinds,end);
        if(min(test)<1E4)
            warning('Ratio of largest to smallest terms in GPD model sum is <1E4.  May need more terms.');
            x;
        end

        s = sum(sumterms1,3)- sum(sumterms2,3);
        s = reshape(s,[l_q,l_a]);

        if(min(s)<0)
            warning('Negative sums found in GPD sum.  Setting to zero.');
            s(find(s<0))=0;
            x;
        end

        logE = -0.5.*GAMMA^2.*s;
        ePerp = exp(logE);

        E_r_matrix = ePar.*ePerp; 
        E_r=sum(E_r_matrix.*weight_matrix,2);
        E=E_r;

        
    end
end

 % Compute the Jacobian matrix; computed numerically
if(nargout>1)
    dx = 0.0001;
    J = zeros(length(E),length(x));
    if nargin < 3 
         
        for i = 1:length(x); % compute the derivatives for all model parameters
            xpert = x;
            xpert(i) = xpert(i)*(1+dx);
            Epert = FiniteCylinder_GPD_dSWOGSE(xpert, protocol);
            dEdx = (Epert - E)/(xpert(i)*dx);
            J(:,i) = dEdx;
        end
    else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
       
        for i = 1:length(x_deriv);
            if x_deriv(i) ~= 0  
            
                xpert = x;
                xpert(i) = xpert(i)*(1+dx);
                Epert = FiniteCylinder_GPD_dSWOGSE(xpert, protocol);
                dEdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEdx;
            end
        end
       

    end   
 end
