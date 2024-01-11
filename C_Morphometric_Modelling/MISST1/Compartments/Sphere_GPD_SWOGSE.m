function [E,J]=Sphere_GPD_SWOGSE(x, protocol, x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a sphere compartment.
% 
% [E,J]=Sphere_GPD_SWOGSE(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate consisting of impermeable 
% spheres with a single radius and a diffusion protocol specified in the input
% Substrate: impermeable spheres with a single radius
% Diffusion pulse sequence: Square wave oscillating gradient spin echo (SWOGSE)
% Signal approximation: Gaussian Phase Distribution (GPD)  
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% x - size 2 vector of model parameters in SI units for Sphere:
%       x(1) - free diffusivity of the material inside the sphere.
%       x(2) - sphere radius 
% protocol - structure which includes all the information related to the 
%        diffusion protocol and required to generate the signal.    
%       protocol.grad_dirs - is the gradient direction of the gradient for
%       each measurement. It has size [N 3] where N is the number of measurements.
%       protocol.G - gradient strength, size [1 N]
%       protocol.delta - pulse separation, size [1 N]
%       protocol.smalldel - pulse duration, size [1 N]
%       protocol.omega - gradient angular frequency, size [1 N]
%       protocol.roots_sphere- solutions to the Bessel function equation from 
%       function BesselJ_RootsSphere.m
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

roots = protocol.roots_sphere;

% Check the roots array is correct
if(abs(roots(1) - 2.0816)>0.0001)
    error('Looks like the roots array is wrong.  First value should be 1.8412, but is %f', roots(1));
end

dRes=x(1);
R=[x(2)]; 


G = protocol.G';
smalldel = protocol.smalldel';
delta = protocol.delta';
omega = protocol.omega';
grad_dirs = protocol.grad_dirs;


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
Bmat_rep = 2./amSq./(repmat(root_m,[l_q*l_a 1 1]).^2 -2);

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

if ~isfield(protocol,'phase')     
   
    % the number of integer half periods:
    NT = floor(2.*smalldel.*niu+0.00000000001); % needed for numerical precision (if NT is int, the floor rounds it down)
    NTmx_rep = floor(2.*smalldelmx_rep.*freqmx_rep+0.00000000001);
    
    if ~isfield(protocol,'mirror') || protocol.mirror == 0 % the 2nd gradient is the negative of the first one
              
     

        % Find restricted signal decay
       
        %precomputing values that appear often in the expresions; rearrange all formulas to have only negative exponentials
        exp_NT = exp(-amSq.*dRes.*NTmx_rep./(2.*freqmx_rep));
        exp_1 = exp(-amSq.*dRes./(2.*freqmx_rep));
        exp_smalldel_NT = exp(-amSq.*dRes.*(smalldelmx_rep-NTmx_rep./(2.*freqmx_rep)));
        exp_delta_NT = exp(-amSq.*dRes.*(deltamx_rep-NTmx_rep./(2.*freqmx_rep)));
        exp_delta = exp(-amSq.*dRes.*deltamx_rep);
        exp_delta_smalldel = exp(-amSq.*dRes.*(deltamx_rep-smalldelmx_rep));
        exp_smalldel = exp(-amSq.*dRes.*smalldelmx_rep);
        am_dRes = amSq.*dRes;
        am2_dRes2 = amSq.^2.*dRes.^2;
        sgn = (-1).^NTmx_rep;

        % corresponds to the multiplication of the full oscillations with themselves
        sumterms = 4.*(-1 + exp_smalldel_NT + NTmx_rep.*(-1+exp_1) + am_dRes.*smalldelmx_rep)./am2_dRes2;

        % corresponds to the multiplication of the full oscillations with
        % other full oscillations from the same gradient
        sumterms = sumterms - 4.* sgn.*(1-exp_1).^2./am2_dRes2./(1+exp_1).^2.*...
            (exp_NT + sgn.*(NTmx_rep-1)+sgn.*NTmx_rep.*exp_1);

        % corresponds to the multiplication of the full oscillations with
        % other full oscillations from different gradients
        sumterms = sumterms + 2.*exp_delta_NT.*(1-exp_1).^2.*(1-sgn.*exp_NT).*(sgn-exp_NT)./...
            (1+exp_1).^2./am2_dRes2;

        % corresponds to the multiplication of the partial oscillations with themselves
        sumterms = sumterms + 2.*exp_delta_smalldel./am2_dRes2.*(exp_smalldel-exp_NT).*...
            (1-exp_smalldel_NT);

        % corresponds to the multiplication of the partial oscillations
        % with full oscillations
          sumterms = sumterms + 2.*sgn.*(1-exp_1)./(1+exp_1)./am2_dRes2.*(1-exp_smalldel_NT).*...
              (-exp_delta.*exp_NT+sgn.*exp_delta_smalldel.*exp_NT+2.* exp_NT-exp_delta_smalldel-...
              2.*sgn+ sgn.*exp_delta);

        sumterms = Bmat_rep.*sumterms;

        testinds = find(sumterms(:,:,end)>0);
        test = sumterms(testinds,1)./sumterms(testinds,end);
        if(min(abs(test))<1E2)
            warning('Ratio of largest to smallest terms in GPD model sum is <1E2.  May need more terms.');
            x;
        end

        s = sum(sumterms,3);
        s = reshape(s,[l_q,l_a]);
        if(min(s)+1E-12<0)
            warning('Negative sums found in GPD sum.  Setting to zero.');
            s(find(s<0))=0;
            x;
        end
        logE = -0.5.*GAMMA^2*GmxSq.*s;
        eRes= exp(logE);
        E_r=sum(eRes,2);
        E=E_r;

       
     else % the 2nd gradient is the negative of the first one reflected         
    

        % Restricted signal

        %precomputing values that appear often in the expresions:
        % rearrange all formulas to have only negative exponentials
        exp_NT = exp(-amSq.*dRes.*NTmx_rep./(2.*freqmx_rep));
        exp_1 = exp(-amSq.*dRes./(2.*freqmx_rep));
        exp_smalldel_NT = exp(-amSq.*dRes.*(smalldelmx_rep-NTmx_rep./(2.*freqmx_rep)));
        exp_delta_NT = exp(-amSq.*dRes.*(deltamx_rep-NTmx_rep./(2.*freqmx_rep)));
        exp_delta = exp(-amSq.*dRes.*deltamx_rep);
        exp_delta_smalldel = exp(-amSq.*dRes.*(deltamx_rep-smalldelmx_rep));
        am_dRes = amSq.*dRes;
        am2_dRes2 = amSq.^2.*dRes.^2;
        sgn = (-1).^NTmx_rep;

          % corresponds to the multiplication of the full oscillations with themselves
        sumterms = 4.*(-1 + exp_smalldel_NT + NTmx_rep.*(-1+exp_1) + am_dRes.*smalldelmx_rep)./am2_dRes2;
        
        % corresponds to the multiplication of the full oscillations with
        % other full oscillations from the same gradient
        sumterms = sumterms - 4.* sgn.*(1-exp_1).^2./am2_dRes2./(1+exp_1).^2.*...
            (exp_NT + sgn.*(NTmx_rep-1)+sgn.*NTmx_rep.*exp_1);

        % corresponds to the multiplication of the full oscillations with
        % other full oscillations from different gradients
        sumterms = sumterms + 2.*sgn.*exp_delta_NT.*exp_smalldel_NT.*(1-exp_1).^2.*(sgn.*exp_NT-1).*(sgn-exp_NT)./...
            (1+exp_1).^2./am2_dRes2;

         % corresponds to the multiplication of the partial oscillations with themselves
        sumterms = sumterms + 2.*exp_delta_smalldel./am2_dRes2.*(exp_smalldel_NT-1).*...
            (1-exp_smalldel_NT);

        % corresponds to the multiplication of the partial oscillations
        % with full oscillations
    sumterms = sumterms -2./am2_dRes2./(1+exp_1).*(1-exp_smalldel_NT).*(-1+exp_1).*(sgn.*exp_NT-1)+...
          2./am2_dRes2./(1+exp_1).*(-1+exp_1).*(sgn.*exp_NT-1).*(1-exp_smalldel_NT).*exp_delta_NT-...
          2.*sgn./am2_dRes2./(1+exp_1).*(1-exp_1).*(exp_smalldel_NT-1).*(-exp_delta+sgn.*(exp_delta_NT))+...
          2.*sgn./am2_dRes2./(1+exp_1).*(1-exp_1).*(exp_smalldel_NT-1).*(sgn-exp_NT);

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

         logE = -0.5.*GAMMA^2*GmxSq.*s;
        eRes= exp(logE);
        E_r=sum(eRes,2);
        E=E_r;
     end    
else  
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
      
        
        % Restricted signal

        %precomputing values that appear often in the expresions:
        % rearrange all formulas to have only negative exponentials
        exp_NT = exp(-amSq.*dRes.*NTmx_rep./(2.*freqmx_rep));
        exp_1 = exp(-amSq.*dRes./(2.*freqmx_rep));
        exp_smalldel_NT_phdelay = exp(-amSq.*dRes.*(smalldelmx_rep-NTmx_rep./(2.*freqmx_rep)-...
            phdelaymx_rep));
        exp_delta_NT_phdelay = exp(-amSq.*dRes.*(deltamx_rep-NTmx_rep./(2.*freqmx_rep)-...
            phdelaymx_rep));
        exp_phdelay = exp(-amSq.*dRes.*phdelaymx_rep);
        exp_delta_phdelay = exp(-amSq.*dRes.*(deltamx_rep-phdelaymx_rep));
        exp_delta_NT = exp(-amSq.*dRes.*(deltamx_rep-NTmx_rep./(2.*freqmx_rep)));
        exp_delta = exp(-amSq.*dRes.*deltamx_rep);
        exp_delta_smalldel = exp(-amSq.*dRes.*(deltamx_rep-smalldelmx_rep));
        exp_smalldel = exp(-amSq.*dRes.*smalldelmx_rep);
        am_dRes = amSq.*dRes;
        am2_dRes2 = amSq.^2.*dRes.^2;

        sgn = (-1).^NTmx_rep;


        % corresponds to the multiplication of the full oscillations with themselves
        sumterms = 4.*(-1 + exp_smalldel_NT_phdelay + NTmx_rep.*(-1+exp_1) + am_dRes.*(smalldelmx_rep-phdelaymx_rep))./am2_dRes2;

        % corresponds to the multiplication of the partial oscillations due to non zero phase with themselves
        sumterms = sumterms + 4.* (am_dRes.*phdelaymx_rep+exp_phdelay-1)./am2_dRes2;

        % corresponds to the multiplication of the partial oscillations due
        % to non zero phase with full oscillations
        sumterms = sumterms +2.*(-1+exp_1).*(sgn.*exp_NT-1).*(1-exp_phdelay)./am2_dRes2./(exp_1+1).*(exp_delta-2) - ...
            2.*exp_delta_NT_phdelay.*(1-exp_1).*(sgn-exp_NT).*(1-exp_phdelay)./am2_dRes2./(exp_1+1);

        % corresponds to the multiplication of the partial oscillations due
        % to non zero phase with the partial oscillations in the end of the
        % pulse
        sumterms = sumterms - 2.*sgn.*(2.*exp_NT-exp_delta.*exp_NT-exp_delta_smalldel).*(1-exp_phdelay).*...
            (1-exp_smalldel_NT_phdelay)./am2_dRes2;

        % corresponds to the multiplication of the partial oscillations due
        % to non zero phase from the two gradients
        sumterms = sumterms - 2.*exp_delta_phdelay.*(1-exp_phdelay).^2./am2_dRes2;

        % corresponds to the multiplication of the full oscillations with
        % other full oscillations from the same gradient
        sumterms = sumterms - 4.* sgn.*(1-exp_1).^2./am2_dRes2./(1+exp_1).^2.*...
            (exp_NT + sgn.*(NTmx_rep-1)+sgn.*NTmx_rep.*exp_1);

        % corresponds to the multiplication of the full oscillations with
        % other full oscillations from different gradients 
        sumterms = sumterms + 2.*exp_delta_NT.*(1-exp_1).^2.*(1-sgn.*exp_NT).*(sgn-exp_NT)./...
            (1+exp_1).^2./am2_dRes2;
      
         % corresponds to the multiplication of the partial oscillations with themselves
        sumterms = sumterms + 2.*exp_delta_smalldel./am2_dRes2.*(exp_smalldel-exp_NT.* exp_phdelay).*...
            (1-exp_smalldel_NT_phdelay);
     
         % corresponds to the multiplication of the partial oscillations
        % with full oscillations
          sumterms = sumterms - 2.*sgn.*(1-exp_1)./(1+exp_1)./am2_dRes2.*(1-exp_smalldel_NT_phdelay).*...
             (-sgn.*exp_delta_smalldel.*exp_phdelay.*exp_NT+exp_delta_smalldel.*exp_phdelay+exp_delta.*exp_NT-...
             2.*exp_NT-sgn.*exp_delta+2.*sgn);

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

        logE = -0.5.*GAMMA^2*GmxSq.*s;
        eRes= exp(logE);
        E_r=sum(eRes,2);
        E=E_r;

    else % the 2nd gradient is the negative of the first one reflected 
      
        % Restricted signal

        %precomputing values that appear often in the expresions:
        % rearrange all formulas to have only negative exponentials
        exp_NT = exp(-amSq.*dRes.*NTmx_rep./(2.*freqmx_rep));
        exp_1 = exp(-amSq.*dRes./(2.*freqmx_rep));
        exp_smalldel_NT_phdelay = exp(-amSq.*dRes.*(smalldelmx_rep-NTmx_rep./(2.*freqmx_rep)-...
            phdelaymx_rep));
        exp_delta_NT_phdelay = exp(-amSq.*dRes.*(deltamx_rep-NTmx_rep./(2.*freqmx_rep)-...
            phdelaymx_rep));
        exp_phdelay = exp(-amSq.*dRes.*phdelaymx_rep);
        exp_delta_phdelay = exp(-amSq.*dRes.*(deltamx_rep-phdelaymx_rep));
        exp_smalldel_phdelay = exp(-amSq.*dRes.*(smalldelmx_rep-phdelaymx_rep));
        exp_delta_smalldel = exp(-amSq.*dRes.*(deltamx_rep-smalldelmx_rep));       
        am_dRes = amSq.*dRes;
        am2_dRes2 = amSq.^2.*dRes.^2;

        sgn = (-1).^NTmx_rep;


        % corresponds to the multiplication of the full oscillations with themselves
        sumterms = 4.*(-1 + exp_smalldel_NT_phdelay + NTmx_rep.*(-1+exp_1) + am_dRes.*(smalldelmx_rep-phdelaymx_rep))./am2_dRes2;

        % corresponds to the multiplication of the partial oscillations due to non zero phase with themselves
        sumterms = sumterms + 4.* (am_dRes.*phdelaymx_rep+exp_phdelay-1)./am2_dRes2;

         % corresponds to the multiplication of the partial oscillations due
        % to non zero phase with full oscillations
        sumterms = sumterms -2.*(exp_1-1).*(sgn.*exp_NT-1).*(1-exp_phdelay)./am2_dRes2./(1+exp_1)-...
            2.*exp_delta_phdelay.*exp_smalldel_NT_phdelay.*(1-exp_1).*(sgn-exp_NT).*(exp_phdelay-1)./am2_dRes2./(1+exp_1) +...
            2.*sgn.*exp_delta_phdelay.*exp_smalldel_NT_phdelay.*(exp_1-1).*(sgn.*exp_NT-1).*(1-exp_phdelay)./am2_dRes2./(1+exp_1) -...
            2.*sgn.*(1-exp_1).*(sgn-exp_NT).*(1-exp_phdelay)./am2_dRes2./(1+exp_1);

         % corresponds to the multiplication of the partial oscillations due
        % to non zero phase with the partial oscillations in the end of the
        % pulse
        sumterms = sumterms + 4.*sgn./am2_dRes2.*(1-exp_phdelay).*(exp_delta_phdelay+exp_smalldel_phdelay-...
            exp_delta_NT_phdelay.*exp_smalldel_phdelay-exp_NT);

        % corresponds to the multiplication of the partial oscillations due
        % to non zero phase from the two gradients
        sumterms = sumterms - 2.*exp_delta_phdelay.*exp_smalldel_phdelay.*(1-exp_phdelay).^2./am2_dRes2;

        % corresponds to the multiplication of the full oscillations with
        % other full oscillations from the same gradient 
        sumterms = sumterms - 4.* sgn.*(1-exp_1).^2./am2_dRes2./(1+exp_1).^2.*...
            (exp_NT + sgn.*(NTmx_rep-1)+sgn.*NTmx_rep.*exp_1);


        % corresponds to the multiplication of the full oscillations with
        % other full oscillations from different gradients  
        sumterms = sumterms + 2.*sgn.*exp_delta_NT_phdelay.*exp_smalldel_NT_phdelay.*(1-exp_1).^2.*...
            (sgn.*exp_NT-1).*(sgn-exp_NT)./am2_dRes2./(1+exp_1).^2;


         % corresponds to the multiplication of the partial oscillations with themselves
        sumterms = sumterms + 2.*exp_delta_smalldel./am2_dRes2.*(exp_smalldel_NT_phdelay-1).*...
            (1-exp_smalldel_NT_phdelay);

        % corresponds to the multiplication of the partial oscillations
        % with full oscillations
          sumterms = sumterms + 2.*sgn.*(1-exp_1)./am2_dRes2./(1+exp_1).*(-exp_delta_NT_phdelay.*(sgn-exp_NT).*(-1+exp_smalldel_NT_phdelay) +...
              sgn.*(sgn.*exp_NT-1).*(1-exp_smalldel_NT_phdelay)+ ...
              (sgn-exp_NT).*(exp_smalldel_NT_phdelay-1)-...
              sgn.*exp_delta_NT_phdelay.*(sgn.*exp_NT-1).*(1-exp_smalldel_NT_phdelay));

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

        logE = -0.5.*GAMMA^2*GmxSq.*s;
        eRes= exp(logE);
        E_r=sum(eRes,2);
        E=E_r;
        
    end
end

 % Compute the Jacobian matrix; computed numerically
if(nargout>1)
    dx = 0.00001;
    J = zeros(length(E), 2);
    if nargin < 3 
         
        for i = 1:2; % compute the derivatives for all model parameters
            xpert = x;
            xpert(i) = xpert(i)*(1+dx);
            Epert = Sphere_GPD_SWOGSE(xpert, protocol);
            dEdx = (Epert - E)/(xpert(i)*dx);
            J(:,i) = dEdx;
        end
    else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
       
        for i = 1:length(x_deriv);
            if x_deriv(i) ~= 0  
                
                xpert = x;
                xpert(i) = xpert(i)*(1+dx);
                Epert = Sphere_GPD_SWOGSE(xpert, protocol);
                dEdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEdx;
            end
        end
       

    end   
 end
