function [E,J]=CuboidSym_GPD_SWOGSE(x, protocol, x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a cuboid/ensemble of cuboids with 2 equal
% sides
% 
% [E,J]=CuboidSym_GPD_SWOGSE(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate consisting of a cuboid/ 
% an ensemble of cuboids with 2 equal sides and a diffusion protocol 
% specified in the input
% Substrate: cuboid / ensemble of cuboids with 2 equal sides
% Diffusion pulse sequence: Square wave oscillating gradient spin echo (SWOGSE)
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
%       protocol.roots_plane - roots of the diffusion equation for planar
%       geometry
%       protocol.cube_rotation - cell array which stores rotation matrices
%       describing the orientation of the cuboids. If it has more than one
%       matrix, than the average is computed; The default is the identity
%       matrix.
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

if ~isfield(protocol,'phase')  
   
    % the number of integer half periods:
    NT = floor(2.*smalldel.*niu+0.00000000001); % needed for numerical precision (if NT is int, the floor rounds it down)
    NTmx_rep = floor(2.*smalldelmx_rep.*freqmx_rep+0.00000000001);
    
     
    %precomputing values that appear often in the expresions; rearrange all formulas to have only negative exponentials
    exp_NT_x = exp(-lambda1_mat.*dRes.*NTmx_rep./(2.*freqmx_rep));
    exp_1_x = exp(-lambda1_mat.*dRes./(2.*freqmx_rep));
    exp_smalldel_NT_x = exp(-lambda1_mat.*dRes.*(smalldelmx_rep-NTmx_rep./(2.*freqmx_rep)));
    exp_delta_NT_x = exp(-lambda1_mat.*dRes.*(deltamx_rep-NTmx_rep./(2.*freqmx_rep)));
    exp_delta_x = exp(-lambda1_mat.*dRes.*deltamx_rep);
    exp_delta_smalldel_x = exp(-lambda1_mat.*dRes.*(deltamx_rep-smalldelmx_rep));
    exp_smalldel_x = exp(-lambda1_mat.*dRes.*smalldelmx_rep);
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
    am_dRes_y = lambda2_mat.*dRes;
    am2_dRes2_y = lambda2_mat.^2.*dRes.^2;

    exp_NT_z = exp(-lambda3_mat.*dRes.*NTmx_rep./(2.*freqmx_rep));
    exp_1_z = exp(-lambda3_mat.*dRes./(2.*freqmx_rep));
    exp_smalldel_NT_z = exp(-lambda3_mat.*dRes.*(smalldelmx_rep-NTmx_rep./(2.*freqmx_rep)));
    exp_delta_NT_z = exp(-lambda3_mat.*dRes.*(deltamx_rep-NTmx_rep./(2.*freqmx_rep)));
    exp_delta_z = exp(-lambda3_mat.*dRes.*deltamx_rep);
    exp_delta_smalldel_z = exp(-lambda3_mat.*dRes.*(deltamx_rep-smalldelmx_rep));
    exp_smalldel_z = exp(-lambda3_mat.*dRes.*smalldelmx_rep);
    am_dRes_z = lambda3_mat.*dRes;
    am2_dRes2_z = lambda3_mat.^2.*dRes.^2;
    
    if ~isfield(protocol,'mirror') || protocol.mirror == 0 % the 2nd gradient is the negative of the first one
         
          
        

        % corresponds to the multiplication of the full oscillations with themselves
        sumterms_x = 4.*(-1 + exp_smalldel_NT_x + NTmx_rep.*(-1+exp_1_x) + am_dRes_x.*smalldelmx_rep)./am2_dRes2_x;
        sumterms_y = 4.*(-1 + exp_smalldel_NT_y + NTmx_rep.*(-1+exp_1_y) + am_dRes_y.*smalldelmx_rep)./am2_dRes2_y;
        sumterms_z = 4.*(-1 + exp_smalldel_NT_z + NTmx_rep.*(-1+exp_1_z) + am_dRes_z.*smalldelmx_rep)./am2_dRes2_z;

        % corresponds to the multiplication of the full oscillations with
        % other full oscillations from the same gradient
        sumterms_x = sumterms_x - 4.* sgn.*(1-exp_1_x).^2./am2_dRes2_x./(1+exp_1_x).^2.*...
            (exp_NT_x + sgn.*(NTmx_rep-1)+sgn.*NTmx_rep.*exp_1_x);
         sumterms_y = sumterms_y - 4.* sgn.*(1-exp_1_y).^2./am2_dRes2_y./(1+exp_1_y).^2.*...
            (exp_NT_y + sgn.*(NTmx_rep-1)+sgn.*NTmx_rep.*exp_1_y);
         sumterms_z = sumterms_z - 4.* sgn.*(1-exp_1_z).^2./am2_dRes2_z./(1+exp_1_z).^2.*...
            (exp_NT_z + sgn.*(NTmx_rep-1)+sgn.*NTmx_rep.*exp_1_z);

        % corresponds to the multiplication of the full oscillations with
        % other full oscillations from different gradients
        sumterms_x = sumterms_x + 2.*exp_delta_NT_x.*(1-exp_1_x).^2.*(1-sgn.*exp_NT_x).*(sgn-exp_NT_x)./...
            (1+exp_1_x).^2./am2_dRes2_x;
        sumterms_y = sumterms_y + 2.*exp_delta_NT_y.*(1-exp_1_y).^2.*(1-sgn.*exp_NT_y).*(sgn-exp_NT_y)./...
            (1+exp_1_y).^2./am2_dRes2_y;
        sumterms_z = sumterms_z + 2.*exp_delta_NT_z.*(1-exp_1_z).^2.*(1-sgn.*exp_NT_z).*(sgn-exp_NT_z)./...
            (1+exp_1_z).^2./am2_dRes2_z;

        % corresponds to the multiplication of the partial oscillations with themselves
        sumterms_x = sumterms_x + 2.*exp_delta_smalldel_x./am2_dRes2_x.*(exp_smalldel_x-exp_NT_x).*...
            (1-exp_smalldel_NT_x);
        sumterms_y = sumterms_y + 2.*exp_delta_smalldel_y./am2_dRes2_y.*(exp_smalldel_y-exp_NT_y).*...
            (1-exp_smalldel_NT_y);
        sumterms_z = sumterms_z + 2.*exp_delta_smalldel_z./am2_dRes2_z.*(exp_smalldel_z-exp_NT_z).*...
            (1-exp_smalldel_NT_z);

        % corresponds to the multiplication of the partial oscillations
        % with full oscillations
          sumterms_x = sumterms_x + 2.*sgn.*(1-exp_1_x)./(1+exp_1_x)./am2_dRes2_x.*(1-exp_smalldel_NT_x).*...
              (-exp_delta_x.*exp_NT_x+sgn.*exp_delta_smalldel_x.*exp_NT_x+2.* exp_NT_x-exp_delta_smalldel_x-...
              2.*sgn+ sgn.*exp_delta_x);
          sumterms_y = sumterms_y + 2.*sgn.*(1-exp_1_y)./(1+exp_1_y)./am2_dRes2_y.*(1-exp_smalldel_NT_y).*...
              (-exp_delta_y.*exp_NT_y+sgn.*exp_delta_smalldel_y.*exp_NT_y+2.* exp_NT_y-exp_delta_smalldel_y-...
              2.*sgn+ sgn.*exp_delta_y);
          sumterms_z = sumterms_z + 2.*sgn.*(1-exp_1_z)./(1+exp_1_z)./am2_dRes2_z.*(1-exp_smalldel_NT_z).*...
              (-exp_delta_z.*exp_NT_z+sgn.*exp_delta_smalldel_z.*exp_NT_z+2.* exp_NT_z-exp_delta_smalldel_z-...
              2.*sgn+ sgn.*exp_delta_z);

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

       
     else % the 2nd gradient is the negative of the first one reflected 
         
       
%         % Perpendicular component
% 
%         %precomputing values that appear often in the expresions:
%         % rearrange all formulas to have only negative exponentials
%         exp_NT_x = exp(-lambda1_mat.*dRes.*NTmx_rep./(2.*freqmx_rep));
%         exp_1_x = exp(-lambda1_mat.*dRes./(2.*freqmx_rep));
%         exp_smalldel_NT_x = exp(-lambda1_mat.*dRes.*(smalldelmx_rep-NTmx_rep./(2.*freqmx_rep)));
%         exp_delta_NT_x = exp(-lambda1_mat.*dRes.*(deltamx_rep-NTmx_rep./(2.*freqmx_rep)));
%         exp_delta_x = exp(-lambda1_mat.*dRes.*deltamx_rep);
%         exp_delta_smalldel_x = exp(-lambda1_mat.*dRes.*(deltamx_rep-smalldelmx_rep));
%         am_dRes_x = lambda1_mat.*dRes;
%         am2_dRes2_x = lambda1_mat.^2.*dRes.^2;
%         sgn = (-1).^NTmx_rep;

          % corresponds to the multiplication of the full oscillations with themselves
        sumterms_x = 4.*(-1 + exp_smalldel_NT_x + NTmx_rep.*(-1+exp_1_x) + am_dRes_x.*smalldelmx_rep)./am2_dRes2_x;
        sumterms_y = 4.*(-1 + exp_smalldel_NT_y + NTmx_rep.*(-1+exp_1_y) + am_dRes_y.*smalldelmx_rep)./am2_dRes2_y;
        sumterms_z = 4.*(-1 + exp_smalldel_NT_z + NTmx_rep.*(-1+exp_1_z) + am_dRes_z.*smalldelmx_rep)./am2_dRes2_z;
        
        % corresponds to the multiplication of the full oscillations with
        % other full oscillations from the same gradient
        sumterms_x = sumterms_x - 4.* sgn.*(1-exp_1_x).^2./am2_dRes2_x./(1+exp_1_x).^2.*...
            (exp_NT_x + sgn.*(NTmx_rep-1)+sgn.*NTmx_rep.*exp_1_x);
        sumterms_y = sumterms_y - 4.* sgn.*(1-exp_1_y).^2./am2_dRes2_y./(1+exp_1_y).^2.*...
            (exp_NT_y + sgn.*(NTmx_rep-1)+sgn.*NTmx_rep.*exp_1_y);
        sumterms_z = sumterms_z - 4.* sgn.*(1-exp_1_z).^2./am2_dRes2_z./(1+exp_1_z).^2.*...
            (exp_NT_z + sgn.*(NTmx_rep-1)+sgn.*NTmx_rep.*exp_1_z);

        % corresponds to the multiplication of the full oscillations with
        % other full oscillations from different gradients
        sumterms_x = sumterms_x + 2.*sgn.*exp_delta_NT_x.*exp_smalldel_NT_x.*(1-exp_1_x).^2.*(sgn.*exp_NT_x-1).*(sgn-exp_NT_x)./...
            (1+exp_1_x).^2./am2_dRes2_x;
        sumterms_y = sumterms_y + 2.*sgn.*exp_delta_NT_y.*exp_smalldel_NT_y.*(1-exp_1_y).^2.*(sgn.*exp_NT_y-1).*(sgn-exp_NT_y)./...
            (1+exp_1_y).^2./am2_dRes2_y;
        sumterms_z = sumterms_z + 2.*sgn.*exp_delta_NT_z.*exp_smalldel_NT_z.*(1-exp_1_z).^2.*(sgn.*exp_NT_z-1).*(sgn-exp_NT_z)./...
            (1+exp_1_z).^2./am2_dRes2_z;

         % corresponds to the multiplication of the partial oscillations with themselves
        sumterms_x = sumterms_x + 2.*exp_delta_smalldel_x./am2_dRes2_x.*(exp_smalldel_NT_x-1).*...
            (1-exp_smalldel_NT_x);
        sumterms_y = sumterms_y + 2.*exp_delta_smalldel_y./am2_dRes2_y.*(exp_smalldel_NT_y-1).*...
            (1-exp_smalldel_NT_y);
        sumterms_z = sumterms_z + 2.*exp_delta_smalldel_z./am2_dRes2_z.*(exp_smalldel_NT_z-1).*...
            (1-exp_smalldel_NT_z);

        % corresponds to the multiplication of the partial oscillations
        % with full oscillations
        sumterms_x = sumterms_x -2./am2_dRes2_x./(1+exp_1_x).*(1-exp_smalldel_NT_x).*(-1+exp_1_x).*(sgn.*exp_NT_x-1)+...
      2./am2_dRes2_x./(1+exp_1_x).*(-1+exp_1_x).*(sgn.*exp_NT_x-1).*(1-exp_smalldel_NT_x).*exp_delta_NT_x-...
      2.*sgn./am2_dRes2_x./(1+exp_1_x).*(1-exp_1_x).*(exp_smalldel_NT_x-1).*(-exp_delta_x+sgn.*(exp_delta_NT_x))+...
      2.*sgn./am2_dRes2_x./(1+exp_1_x).*(1-exp_1_x).*(exp_smalldel_NT_x-1).*(sgn-exp_NT_x);
     sumterms_y = sumterms_y -2./am2_dRes2_y./(1+exp_1_y).*(1-exp_smalldel_NT_y).*(-1+exp_1_y).*(sgn.*exp_NT_y-1)+...
      2./am2_dRes2_y./(1+exp_1_y).*(-1+exp_1_y).*(sgn.*exp_NT_y-1).*(1-exp_smalldel_NT_y).*exp_delta_NT_y-...
      2.*sgn./am2_dRes2_y./(1+exp_1_y).*(1-exp_1_y).*(exp_smalldel_NT_y-1).*(-exp_delta_y+sgn.*(exp_delta_NT_y))+...
      2.*sgn./am2_dRes2_y./(1+exp_1_y).*(1-exp_1_y).*(exp_smalldel_NT_y-1).*(sgn-exp_NT_y);
     sumterms_z = sumterms_z -2./am2_dRes2_z./(1+exp_1_z).*(1-exp_smalldel_NT_z).*(-1+exp_1_z).*(sgn.*exp_NT_z-1)+...
      2./am2_dRes2_z./(1+exp_1_z).*(-1+exp_1_z).*(sgn.*exp_NT_z-1).*(1-exp_smalldel_NT_z).*exp_delta_NT_z-...
      2.*sgn./am2_dRes2_z./(1+exp_1_z).*(1-exp_1_z).*(exp_smalldel_NT_z-1).*(-exp_delta_z+sgn.*(exp_delta_NT_z))+...
      2.*sgn./am2_dRes2_z./(1+exp_1_z).*(1-exp_1_z).*(exp_smalldel_NT_z-1).*(sgn-exp_NT_z);
        

        sumterms_x = B1_mat.*sumterms_x;
        sumterms_y = B2_mat.*sumterms_y;
        sumterms_z = B3_mat.*sumterms_z;

        testinds = find(sumterms_x(:,:,end)>0);
        test = sumterms_x(testinds,1)./sumterms_x(testinds,end);

        if(min(test)<1E4)
            warning('Ratio of largest to smallest terms in GPD model sum is <1E4.  May need more terms.');
            x;
        end

        s_x = sum(sumterms_x,3);
        s_y = sum(sumterms_y,3);
        s_z = sum(sumterms_z,3);
        s_x = reshape(s_x,[l_q,l_a]);
        s_y = reshape(s_y,[l_q,l_a]);
        s_z = reshape(s_z,[l_q,l_a]);

        if(min(s_x)<0)
            warning('Negative sums found in GPD sum.  Setting to zero.');
            s_x(find(s_x<0))=0;
            x;
        end 

        logEx = - 0.5.*GAMMA.^2.*Gxmx.^2.*s_x;
        logEy = - 0.5.*GAMMA.^2.*Gymx.^2.*s_y;
        logEz = - 0.5.*GAMMA.^2.*Gzmx.^2.*s_z;
        
        logE = logEx+logEy+logEz;
        E_r_matrix = exp(logE);
        E_r=sum(E_r_matrix,2);
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
    
     %precomputing values that appear often in the expresions:
    % rearrange all formulas to have only negative exponentials
    exp_NT_x = exp(-lambda1_mat.*dRes.*NTmx_rep./(2.*freqmx_rep));
    exp_1_x = exp(-lambda1_mat.*dRes./(2.*freqmx_rep));
    exp_smalldel_NT_phdelay_x = exp(-lambda1_mat.*dRes.*(smalldelmx_rep-NTmx_rep./(2.*freqmx_rep)-...
        phdelaymx_rep));
    exp_delta_NT_phdelay_x = exp(-lambda1_mat.*dRes.*(deltamx_rep-NTmx_rep./(2.*freqmx_rep)-...
        phdelaymx_rep));
    exp_phdelay_x = exp(-lambda1_mat.*dRes.*phdelaymx_rep);
    exp_delta_phdelay_x = exp(-lambda1_mat.*dRes.*(deltamx_rep-phdelaymx_rep));
    exp_smalldel_phdelay_x = exp(-lambda1_mat.*dRes.*(smalldelmx_rep-phdelaymx_rep));
    exp_delta_NT_x = exp(-lambda1_mat.*dRes.*(deltamx_rep-NTmx_rep./(2.*freqmx_rep)));
    exp_delta_x = exp(-lambda1_mat.*dRes.*deltamx_rep);
    exp_delta_smalldel_x = exp(-lambda1_mat.*dRes.*(deltamx_rep-smalldelmx_rep));
    exp_smalldel_x = exp(-lambda1_mat.*dRes.*smalldelmx_rep);
    am_dRes_x = lambda1_mat.*dRes;
    am2_dRes2_x = lambda1_mat.^2.*dRes.^2;
    
     exp_NT_y = exp(-lambda2_mat.*dRes.*NTmx_rep./(2.*freqmx_rep));
    exp_1_y = exp(-lambda2_mat.*dRes./(2.*freqmx_rep));
    exp_smalldel_NT_phdelay_y = exp(-lambda2_mat.*dRes.*(smalldelmx_rep-NTmx_rep./(2.*freqmx_rep)-...
        phdelaymx_rep));
    exp_delta_NT_phdelay_y = exp(-lambda2_mat.*dRes.*(deltamx_rep-NTmx_rep./(2.*freqmx_rep)-...
        phdelaymx_rep));
    exp_phdelay_y = exp(-lambda2_mat.*dRes.*phdelaymx_rep);
    exp_delta_phdelay_y = exp(-lambda2_mat.*dRes.*(deltamx_rep-phdelaymx_rep));
    exp_smalldel_phdelay_y = exp(-lambda2_mat.*dRes.*(smalldelmx_rep-phdelaymx_rep));
    exp_delta_NT_y = exp(-lambda2_mat.*dRes.*(deltamx_rep-NTmx_rep./(2.*freqmx_rep)));
    exp_delta_y = exp(-lambda2_mat.*dRes.*deltamx_rep);
    exp_delta_smalldel_y = exp(-lambda2_mat.*dRes.*(deltamx_rep-smalldelmx_rep));
    exp_smalldel_y = exp(-lambda2_mat.*dRes.*smalldelmx_rep);
    am_dRes_y = lambda2_mat.*dRes;
    am2_dRes2_y = lambda2_mat.^2.*dRes.^2;
    
     exp_NT_z = exp(-lambda3_mat.*dRes.*NTmx_rep./(2.*freqmx_rep));
    exp_1_z = exp(-lambda3_mat.*dRes./(2.*freqmx_rep));
    exp_smalldel_NT_phdelay_z = exp(-lambda3_mat.*dRes.*(smalldelmx_rep-NTmx_rep./(2.*freqmx_rep)-...
        phdelaymx_rep));
    exp_delta_NT_phdelay_z = exp(-lambda3_mat.*dRes.*(deltamx_rep-NTmx_rep./(2.*freqmx_rep)-...
        phdelaymx_rep));
    exp_phdelay_z = exp(-lambda3_mat.*dRes.*phdelaymx_rep);
    exp_delta_phdelay_z = exp(-lambda3_mat.*dRes.*(deltamx_rep-phdelaymx_rep));
    exp_smalldel_phdelay_z = exp(-lambda3_mat.*dRes.*(smalldelmx_rep-phdelaymx_rep));
    exp_delta_NT_z = exp(-lambda3_mat.*dRes.*(deltamx_rep-NTmx_rep./(2.*freqmx_rep)));
    exp_delta_z = exp(-lambda3_mat.*dRes.*deltamx_rep);
    exp_delta_smalldel_z = exp(-lambda3_mat.*dRes.*(deltamx_rep-smalldelmx_rep));
    exp_smalldel_z = exp(-lambda3_mat.*dRes.*smalldelmx_rep);
    am_dRes_z = lambda3_mat.*dRes;
    am2_dRes2_z = lambda3_mat.^2.*dRes.^2;



    sgn = (-1).^NTmx_rep;

    
    if ~isfield(protocol,'mirror') || protocol.mirror == 0 % the 2nd gradient is the negative of the first one
               
        % restricted signal
        % restricted signal
        % corresponds to the multiplication of the full oscillations with themselves
        sumterms_x = 4.*(-1 + exp_smalldel_NT_phdelay_x + NTmx_rep.*(-1+exp_1_x) + am_dRes_x.*(smalldelmx_rep-phdelaymx_rep))./am2_dRes2_x;
        sumterms_y = 4.*(-1 + exp_smalldel_NT_phdelay_y + NTmx_rep.*(-1+exp_1_y) + am_dRes_y.*(smalldelmx_rep-phdelaymx_rep))./am2_dRes2_y;
        sumterms_z = 4.*(-1 + exp_smalldel_NT_phdelay_z + NTmx_rep.*(-1+exp_1_z) + am_dRes_z.*(smalldelmx_rep-phdelaymx_rep))./am2_dRes2_z;

        % corresponds to the multiplication of the partial oscillations due to non zero phase with themselves
        sumterms_x = sumterms_x + 4.* (am_dRes_x.*phdelaymx_rep+exp_phdelay_x-1)./am2_dRes2_x;
        sumterms_y = sumterms_y + 4.* (am_dRes_y.*phdelaymx_rep+exp_phdelay_y-1)./am2_dRes2_y;
        sumterms_z = sumterms_z + 4.* (am_dRes_z.*phdelaymx_rep+exp_phdelay_z-1)./am2_dRes2_z;

        % corresponds to the multiplication of the partial oscillations due
        % to non zero phase with full oscillations
        sumterms_x = sumterms_x +2.*(-1+exp_1_x).*(sgn.*exp_NT_x-1).*(1-exp_phdelay_x)./am2_dRes2_x./(exp_1_x+1).*(exp_delta_x-2) - ...
            2.*exp_delta_NT_phdelay_x.*(1-exp_1_x).*(sgn-exp_NT_x).*(1-exp_phdelay_x)./am2_dRes2_x./(exp_1_x+1);
        sumterms_y = sumterms_y +2.*(-1+exp_1_y).*(sgn.*exp_NT_y-1).*(1-exp_phdelay_y)./am2_dRes2_y./(exp_1_y+1).*(exp_delta_y-2) - ...
            2.*exp_delta_NT_phdelay_y.*(1-exp_1_y).*(sgn-exp_NT_y).*(1-exp_phdelay_y)./am2_dRes2_y./(exp_1_y+1);
        sumterms_z = sumterms_z +2.*(-1+exp_1_z).*(sgn.*exp_NT_z-1).*(1-exp_phdelay_z)./am2_dRes2_z./(exp_1_z+1).*(exp_delta_z-2) - ...
            2.*exp_delta_NT_phdelay_z.*(1-exp_1_z).*(sgn-exp_NT_z).*(1-exp_phdelay_z)./am2_dRes2_z./(exp_1_z+1);

        % corresponds to the multiplication of the partial oscillations due
        % to non zero phase with the partial oscillations in the end of the
        % pulse
        sumterms_x = sumterms_x - 2.*sgn.*(2.*exp_NT_x-exp_delta_x.*exp_NT_x-exp_delta_smalldel_x).*(1-exp_phdelay_x).*...
            (1-exp_smalldel_NT_phdelay_x)./am2_dRes2_x;
        sumterms_y = sumterms_y - 2.*sgn.*(2.*exp_NT_y-exp_delta_y.*exp_NT_y-exp_delta_smalldel_y).*(1-exp_phdelay_y).*...
            (1-exp_smalldel_NT_phdelay_y)./am2_dRes2_y;
        sumterms_z = sumterms_z - 2.*sgn.*(2.*exp_NT_z-exp_delta_z.*exp_NT_z-exp_delta_smalldel_z).*(1-exp_phdelay_z).*...
            (1-exp_smalldel_NT_phdelay_z)./am2_dRes2_z;

        % corresponds to the multiplication of the partial oscillations due
        % to non zero phase from the two gradients
        sumterms_x = sumterms_x - 2.*exp_delta_phdelay_x.*(1-exp_phdelay_x).^2./am2_dRes2_x;
        sumterms_y = sumterms_y - 2.*exp_delta_phdelay_y.*(1-exp_phdelay_y).^2./am2_dRes2_y;
        sumterms_z = sumterms_z - 2.*exp_delta_phdelay_z.*(1-exp_phdelay_z).^2./am2_dRes2_z;

        % corresponds to the multiplication of the full oscillations with
        % other full oscillations from the same gradient
        sumterms_x = sumterms_x - 4.* sgn.*(1-exp_1_x).^2./am2_dRes2_x./(1+exp_1_x).^2.*...
            (exp_NT_x + sgn.*(NTmx_rep-1)+sgn.*NTmx_rep.*exp_1_x);
        sumterms_y = sumterms_y - 4.* sgn.*(1-exp_1_y).^2./am2_dRes2_y./(1+exp_1_y).^2.*...
            (exp_NT_y + sgn.*(NTmx_rep-1)+sgn.*NTmx_rep.*exp_1_y);
        sumterms_z = sumterms_z - 4.* sgn.*(1-exp_1_z).^2./am2_dRes2_z./(1+exp_1_z).^2.*...
            (exp_NT_z + sgn.*(NTmx_rep-1)+sgn.*NTmx_rep.*exp_1_z);

        % corresponds to the multiplication of the full oscillations with
        % other full oscillations from different gradients 
        sumterms_x = sumterms_x + 2.*exp_delta_NT_x.*(1-exp_1_x).^2.*(1-sgn.*exp_NT_x).*(sgn-exp_NT_x)./...
            (1+exp_1_x).^2./am2_dRes2_x;
        sumterms_y = sumterms_y + 2.*exp_delta_NT_y.*(1-exp_1_y).^2.*(1-sgn.*exp_NT_y).*(sgn-exp_NT_y)./...
            (1+exp_1_y).^2./am2_dRes2_y;
        sumterms_z = sumterms_z + 2.*exp_delta_NT_z.*(1-exp_1_z).^2.*(1-sgn.*exp_NT_z).*(sgn-exp_NT_z)./...
            (1+exp_1_z).^2./am2_dRes2_z;
      
         % corresponds to the multiplication of the partial oscillations with themselves
        sumterms_x = sumterms_x + 2.*exp_delta_smalldel_x./am2_dRes2_x.*(exp_smalldel_x-exp_NT_x.* exp_phdelay_x).*...
            (1-exp_smalldel_NT_phdelay_x);
        sumterms_y = sumterms_y + 2.*exp_delta_smalldel_y./am2_dRes2_y.*(exp_smalldel_y-exp_NT_y.* exp_phdelay_y).*...
            (1-exp_smalldel_NT_phdelay_y);
        sumterms_z = sumterms_z + 2.*exp_delta_smalldel_z./am2_dRes2_z.*(exp_smalldel_z-exp_NT_z.* exp_phdelay_z).*...
            (1-exp_smalldel_NT_phdelay_z);
     
         % corresponds to the multiplication of the partial oscillations
        % with full oscillations
          sumterms_x = sumterms_x - 2.*sgn.*(1-exp_1_x)./(1+exp_1_x)./am2_dRes2_x.*(1-exp_smalldel_NT_phdelay_x).*...
             (-sgn.*exp_delta_smalldel_x.*exp_phdelay_x.*exp_NT_x+exp_delta_smalldel_x.*exp_phdelay_x+exp_delta_x.*exp_NT_x-...
             2.*exp_NT_x-sgn.*exp_delta_x+2.*sgn);
         sumterms_y = sumterms_y - 2.*sgn.*(1-exp_1_y)./(1+exp_1_y)./am2_dRes2_y.*(1-exp_smalldel_NT_phdelay_y).*...
             (-sgn.*exp_delta_smalldel_y.*exp_phdelay_y.*exp_NT_y+exp_delta_smalldel_y.*exp_phdelay_y+exp_delta_y.*exp_NT_y-...
             2.*exp_NT_y-sgn.*exp_delta_y+2.*sgn);
         sumterms_z = sumterms_z - 2.*sgn.*(1-exp_1_z)./(1+exp_1_z)./am2_dRes2_z.*(1-exp_smalldel_NT_phdelay_z).*...
             (-sgn.*exp_delta_smalldel_z.*exp_phdelay_z.*exp_NT_z+exp_delta_smalldel_z.*exp_phdelay_z+exp_delta_z.*exp_NT_z-...
             2.*exp_NT_z-sgn.*exp_delta_z+2.*sgn);

        sumterms_x = B1_mat.*sumterms_x;
        sumterms_y = B2_mat.*sumterms_y;
        sumterms_z = B3_mat.*sumterms_z;

        testinds = find(sumterms_x(:,:,end)>0);
        test = sumterms_x(testinds,1)./sumterms_x(testinds,end);
        if(min(test)<1E4)
            warning('Ratio of largest to smallest terms in GPD model sum is <1E4.  May need more terms.');
            x;
        end

        s_x = sum(sumterms_x,3);
        s_x = reshape(s_x,[l_q,l_a]);
        s_y = sum(sumterms_y,3);
        s_y = reshape(s_y,[l_q,l_a]);
        s_z = sum(sumterms_z,3);
        s_z = reshape(s_z,[l_q,l_a]);

        if(min(s_x)<0)
            warning('Negative sums found in GPD sum.  Setting to zero.');
            s_x(find(s_x<0))=0;
            x;
        end

        logEx = - 0.5.*GAMMA.^2.*Gxmx.^2.*s_x;
        logEy = - 0.5.*GAMMA.^2.*Gymx.^2.*s_y;
        logEz = - 0.5.*GAMMA.^2.*Gzmx.^2.*s_z;
        
        logE = logEx+logEy+logEz;
        E_r_matrix = exp(logE);
        E_r=sum(E_r_matrix,2);
        E=E_r;

    else % the 2nd gradient is the negative of the first one reflected 
      
         % corresponds to the multiplication of the full oscillations with themselves
        sumterms_x = 4.*(-1 + exp_smalldel_NT_phdelay_x + NTmx_rep.*(-1+exp_1_x) + am_dRes_x.*(smalldelmx_rep-phdelaymx_rep))./am2_dRes2_x;
        sumterms_y = 4.*(-1 + exp_smalldel_NT_phdelay_y + NTmx_rep.*(-1+exp_1_y) + am_dRes_y.*(smalldelmx_rep-phdelaymx_rep))./am2_dRes2_y;
        sumterms_z = 4.*(-1 + exp_smalldel_NT_phdelay_z + NTmx_rep.*(-1+exp_1_z) + am_dRes_z.*(smalldelmx_rep-phdelaymx_rep))./am2_dRes2_z;

        % corresponds to the multiplication of the partial oscillations due to non zero phase with themselves
        sumterms_x = sumterms_x + 4.* (am_dRes_x.*phdelaymx_rep+exp_phdelay_x-1)./am2_dRes2_x;
        sumterms_y = sumterms_y + 4.* (am_dRes_y.*phdelaymx_rep+exp_phdelay_y-1)./am2_dRes2_y;
        sumterms_z = sumterms_z + 4.* (am_dRes_z.*phdelaymx_rep+exp_phdelay_z-1)./am2_dRes2_z;

         % corresponds to the multiplication of the partial oscillations due
        % to non zero phase with full oscillations
        sumterms_x = sumterms_x -2.*(exp_1_x-1).*(sgn.*exp_NT_x-1).*(1-exp_phdelay_x)./am2_dRes2_x./(1+exp_1_x)-...
            2.*exp_delta_phdelay_x.*exp_smalldel_NT_phdelay_x.*(1-exp_1_x).*(sgn-exp_NT_x).*(exp_phdelay_x-1)./am2_dRes2_x./(1+exp_1_x) +...
            2.*sgn.*exp_delta_phdelay_x.*exp_smalldel_NT_phdelay_x.*(exp_1_x-1).*(sgn.*exp_NT_x-1).*(1-exp_phdelay_x)./am2_dRes2_x./(1+exp_1_x) -...
            2.*sgn.*(1-exp_1_x).*(sgn-exp_NT_x).*(1-exp_phdelay_x)./am2_dRes2_x./(1+exp_1_x);
        sumterms_y = sumterms_y -2.*(exp_1_y-1).*(sgn.*exp_NT_y-1).*(1-exp_phdelay_y)./am2_dRes2_y./(1+exp_1_y)-...
            2.*exp_delta_phdelay_y.*exp_smalldel_NT_phdelay_y.*(1-exp_1_y).*(sgn-exp_NT_y).*(exp_phdelay_y-1)./am2_dRes2_y./(1+exp_1_y) +...
            2.*sgn.*exp_delta_phdelay_y.*exp_smalldel_NT_phdelay_y.*(exp_1_y-1).*(sgn.*exp_NT_y-1).*(1-exp_phdelay_y)./am2_dRes2_y./(1+exp_1_y) -...
            2.*sgn.*(1-exp_1_y).*(sgn-exp_NT_y).*(1-exp_phdelay_y)./am2_dRes2_y./(1+exp_1_y);

        sumterms_z = sumterms_z -2.*(exp_1_z-1).*(sgn.*exp_NT_z-1).*(1-exp_phdelay_z)./am2_dRes2_z./(1+exp_1_z)-...
            2.*exp_delta_phdelay_z.*exp_smalldel_NT_phdelay_z.*(1-exp_1_z).*(sgn-exp_NT_z).*(exp_phdelay_z-1)./am2_dRes2_z./(1+exp_1_z) +...
            2.*sgn.*exp_delta_phdelay_z.*exp_smalldel_NT_phdelay_z.*(exp_1_z-1).*(sgn.*exp_NT_z-1).*(1-exp_phdelay_z)./am2_dRes2_z./(1+exp_1_z) -...
            2.*sgn.*(1-exp_1_z).*(sgn-exp_NT_z).*(1-exp_phdelay_z)./am2_dRes2_z./(1+exp_1_z);


         % corresponds to the multiplication of the partial oscillations due
        % to non zero phase with the partial oscillations in the end of the
        % pulse
        sumterms_x = sumterms_x + 4.*sgn./am2_dRes2_x.*(1-exp_phdelay_x).*(exp_delta_phdelay_x+exp_smalldel_phdelay_x-...
            exp_delta_NT_phdelay_x.*exp_smalldel_phdelay_x-exp_NT_x);
        sumterms_y = sumterms_y + 4.*sgn./am2_dRes2_y.*(1-exp_phdelay_y).*(exp_delta_phdelay_y+exp_smalldel_phdelay_y-...
            exp_delta_NT_phdelay_y.*exp_smalldel_phdelay_y-exp_NT_y);
        sumterms_z = sumterms_z + 4.*sgn./am2_dRes2_z.*(1-exp_phdelay_z).*(exp_delta_phdelay_z+exp_smalldel_phdelay_z-...
            exp_delta_NT_phdelay_z.*exp_smalldel_phdelay_z-exp_NT_z);

        % corresponds to the multiplication of the partial oscillations due
        % to non zero phase from the two gradients
        sumterms_x = sumterms_x - 2.*exp_delta_phdelay_x.*exp_smalldel_phdelay_x.*(1-exp_phdelay_x).^2./am2_dRes2_x;
        sumterms_y = sumterms_y - 2.*exp_delta_phdelay_y.*exp_smalldel_phdelay_y.*(1-exp_phdelay_y).^2./am2_dRes2_y;
        sumterms_z = sumterms_z - 2.*exp_delta_phdelay_z.*exp_smalldel_phdelay_z.*(1-exp_phdelay_z).^2./am2_dRes2_z;

        % corresponds to the multiplication of the full oscillations with
        % other full oscillations from the same gradient 
        sumterms_x = sumterms_x - 4.* sgn.*(1-exp_1_x).^2./am2_dRes2_x./(1+exp_1_x).^2.*...
            (exp_NT_x + sgn.*(NTmx_rep-1)+sgn.*NTmx_rep.*exp_1_x);
        sumterms_y = sumterms_y - 4.* sgn.*(1-exp_1_y).^2./am2_dRes2_y./(1+exp_1_y).^2.*...
            (exp_NT_y + sgn.*(NTmx_rep-1)+sgn.*NTmx_rep.*exp_1_y);
        sumterms_z = sumterms_z - 4.* sgn.*(1-exp_1_z).^2./am2_dRes2_z./(1+exp_1_z).^2.*...
            (exp_NT_z + sgn.*(NTmx_rep-1)+sgn.*NTmx_rep.*exp_1_z);


        % corresponds to the multiplication of the full oscillations with
        % other full oscillations from different gradients  
        sumterms_x = sumterms_x + 2.*sgn.*exp_delta_NT_phdelay_x.*exp_smalldel_NT_phdelay_x.*(1-exp_1_x).^2.*...
            (sgn.*exp_NT_x-1).*(sgn-exp_NT_x)./am2_dRes2_x./(1+exp_1_x).^2;
        sumterms_y = sumterms_y + 2.*sgn.*exp_delta_NT_phdelay_y.*exp_smalldel_NT_phdelay_y.*(1-exp_1_y).^2.*...
            (sgn.*exp_NT_y-1).*(sgn-exp_NT_y)./am2_dRes2_y./(1+exp_1_y).^2;
        sumterms_z = sumterms_z + 2.*sgn.*exp_delta_NT_phdelay_z.*exp_smalldel_NT_phdelay_z.*(1-exp_1_z).^2.*...
            (sgn.*exp_NT_z-1).*(sgn-exp_NT_z)./am2_dRes2_z./(1+exp_1_z).^2;


         % corresponds to the multiplication of the partial oscillations with themselves
        sumterms_x = sumterms_x + 2.*exp_delta_smalldel_x./am2_dRes2_x.*(exp_smalldel_NT_phdelay_x-1).*...
            (1-exp_smalldel_NT_phdelay_x);
        sumterms_y = sumterms_y + 2.*exp_delta_smalldel_y./am2_dRes2_y.*(exp_smalldel_NT_phdelay_y-1).*...
            (1-exp_smalldel_NT_phdelay_y);
        sumterms_z = sumterms_z + 2.*exp_delta_smalldel_z./am2_dRes2_z.*(exp_smalldel_NT_phdelay_z-1).*...
            (1-exp_smalldel_NT_phdelay_z);

        % corresponds to the multiplication of the partial oscillations
        % with full oscillations
          sumterms_x = sumterms_x + 2.*sgn.*(1-exp_1_x)./am2_dRes2_x./(1+exp_1_x).*(-exp_delta_NT_phdelay_x.*(sgn-exp_NT_x).*(-1+exp_smalldel_NT_phdelay_x) +...
              sgn.*(sgn.*exp_NT_x-1).*(1-exp_smalldel_NT_phdelay_x)+ ...
              (sgn-exp_NT_x).*(exp_smalldel_NT_phdelay_x-1)-...
              sgn.*exp_delta_NT_phdelay_x.*(sgn.*exp_NT_x-1).*(1-exp_smalldel_NT_phdelay_x));
          sumterms_y = sumterms_y + 2.*sgn.*(1-exp_1_y)./am2_dRes2_y./(1+exp_1_y).*(-exp_delta_NT_phdelay_y.*(sgn-exp_NT_y).*(-1+exp_smalldel_NT_phdelay_y) +...
              sgn.*(sgn.*exp_NT_y-1).*(1-exp_smalldel_NT_phdelay_y)+ ...
              (sgn-exp_NT_y).*(exp_smalldel_NT_phdelay_y-1)-...
              sgn.*exp_delta_NT_phdelay_y.*(sgn.*exp_NT_y-1).*(1-exp_smalldel_NT_phdelay_y));
          sumterms_z = sumterms_z + 2.*sgn.*(1-exp_1_z)./am2_dRes2_z./(1+exp_1_z).*(-exp_delta_NT_phdelay_z.*(sgn-exp_NT_z).*(-1+exp_smalldel_NT_phdelay_z) +...
              sgn.*(sgn.*exp_NT_z-1).*(1-exp_smalldel_NT_phdelay_z)+ ...
              (sgn-exp_NT_z).*(exp_smalldel_NT_phdelay_z-1)-...
              sgn.*exp_delta_NT_phdelay_z.*(sgn.*exp_NT_z-1).*(1-exp_smalldel_NT_phdelay_z));

          sumterms_x = B1_mat.*sumterms_x;
        sumterms_y = B2_mat.*sumterms_y;
        sumterms_z = B3_mat.*sumterms_z;

        testinds = find(sumterms_x(:,:,end)>0);
        test = sumterms_x(testinds,1)./sumterms_x(testinds,end);
        if(min(test)<1E4)
            warning('Ratio of largest to smallest terms in GPD model sum is <1E4.  May need more terms.');
            x;
        end

        s_x = sum(sumterms_x,3);
        s_x = reshape(s_x,[l_q,l_a]);
        s_y = sum(sumterms_y,3);
        s_y = reshape(s_y,[l_q,l_a]);
        s_z = sum(sumterms_z,3);
        s_z = reshape(s_z,[l_q,l_a]);

        if(min(s_x)<0)
            warning('Negative sums found in GPD sum.  Setting to zero.');
            s_x(find(s_x<0))=0;
            x;
        end

        logEx = - 0.5.*GAMMA.^2.*Gxmx.^2.*s_x;
        logEy = - 0.5.*GAMMA.^2.*Gymx.^2.*s_y;
        logEz = - 0.5.*GAMMA.^2.*Gzmx.^2.*s_z;
        
        logE = logEx+logEy+logEz;
        E_r_matrix = exp(logE);
        E_r=sum(E_r_matrix,2);
        E=E_r;

        
    end
end

 % Compute the Jacobian matrix; computed numerically
if(nargout>1)
    dx = 0.00001;
     J = zeros(length(E), 3);
    if nargin < 3 
         
        for i = 1:3; % compute the derivatives for all model parameters
            xpert = x;
            xpert(i) = xpert(i)*(1+dx);
            Epert = CuboidSym_GPD_SWOGSE(xpert, protocol);
            dEdx = (Epert - E)/(xpert(i)*dx);
            J(:,i) = dEdx;
        end
    else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
       
        for i = 1:length(x_deriv);
            if x_deriv(i) ~= 0  
               
                xpert = x;
                xpert(i) = xpert(i)*(1+dx);
                Epert = CuboidSym_GPD_SWOGSE(xpert, protocol);
                dEdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEdx;
            end
        end
       

    end   
 end