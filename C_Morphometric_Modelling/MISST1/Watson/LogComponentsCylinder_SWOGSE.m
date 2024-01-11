function [Lperp,Lpar]=LogComponentsCylinder_SWOGSE(x, protocol)
% Substrate: Parallel, impermeable cylinders with one radius
% Pulse sequence: square wave oscillating gradient spin echo with angular
% frequency omega (and phase  - optional)
% Signal approximation: Gaussian phase distribution.
%
% [Lperp,Lpar]=LogComponentsCylinder_SWOGSE(x, protocol)
% returns the Log of the signal when the gradient is parallel or perpendicular to the fibre direction 
% needed for watson and bingham distributions
%
% x is the list of model parameters in SI units:
% x(1) is the free diffusivity of the material inside the cylinders.
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
%
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
              
        % Caculate bvalue
        bval = GAMMA.^2.*G.^2.*(1./(48.*niu.^3)).*((2.*NT.^3 +... 
             3.*NT.^2.*(1 + (-1).^NT - 4.*niu.*smalldel) - 4.*niu.^2.*smalldel.^2.*...
             (-3 + 3.*(-1).^NT + 4.*niu.*smalldel) +3.*delta.*niu.*(-1 + (-1).^NT - ...
             2.*NT + 4.*niu.*smalldel).^2 + NT.*(1 + 3.*(-1).^NT - 12.*niu.*smalldel + ...
             24.*niu.^2.*smalldel.^2)));
        bval_matrix=repmat(bval,[1,l_a]);

        % Find restricted signal decay
        % Parallel component
        Lpar=-bval_matrix.*dRes;
       
        
        % Perpendicular component

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

        % corresponds to the multiplication of the full and incomplete oscillations with themselves
        sumterms = 4.*(-1 + exp_smalldel_NT + NTmx_rep.*(-1+exp_1) + am_dRes.*smalldelmx_rep)./am2_dRes2;

        % corresponds to the multiplication of the full oscillations with
        % other full oscillations from the same gradient
        sumterms = sumterms - 4.* sgn.*(1-exp_1).^2./am2_dRes2./(1+exp_1).^2.*...
            (exp_NT + sgn.*(NTmx_rep-1)+sgn.*NTmx_rep.*exp_1);

        % corresponds to the multiplication of the full oscillations with
        % other full oscillations from different gradients
        sumterms = sumterms + 2.*exp_delta_NT.*(1-exp_1).^2.*(1-sgn.*exp_NT).*(sgn-exp_NT)./...
            (1+exp_1).^2./am2_dRes2;

        % corresponds to the multiplication of the partial oscillations
        % from the first pulse with the incomplete oscillations from the
        % 2nd pulse
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
        Lperp = -0.5.*GAMMA^2*GmxSq.*s;
        Lperp(G==0) = 0;
        

       
     else % the 2nd gradient is the negative of the first one reflected 
         
         % calculate b value
        bval = GAMMA.^2.*G.^2.*((delta-smalldel).*(smalldel-(0.5.*(1-(-1).^NT)+NT)./(2.*niu)).^2 + NT./(12.*niu.^3) +...
            1./(192.*niu.^3).*(-(-1+(-1).^NT).^3+(-1+(-1).^NT-2.*NT+4.*smalldel.*niu).^3) - ...
            1./(96.*niu.^3).*(NT-2.*smalldel.*niu).*(3.*(-1+(-1).^NT).^2+4.*NT.^2+12.*smalldel.*niu.*(-1+(-1).^NT)+...
            16.*smalldel.^2.*niu.^2-2.*NT.*(-3+3.*(-1).^NT+8.*smalldel.*niu)));
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
        exp_smalldel = exp(-amSq.*dRes.*smalldelmx_rep);
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

        Lperp = -0.5.*GAMMA^2*GmxSq.*s;
        Lperp(G==0) = 0;
 
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
        sgn = (-1).^NT;
        % calculate b value
        bval = GAMMA.^2.*G.^2.*(1./3.*(smalldel-NT./2./niu-phdelay).^3 + (delta-smalldel).*...
        (sgn.*(smalldel- (0.5.*(1-sgn)+NT)./2./niu-phdelay)-phdelay).^2  +(phdelay.^3)./3 +...
        sgn.*((-1+sgn-2.*NT+4.*niu.*(smalldel-phdelay)).^3)./192./(niu.^3) +...   
         NT./96./(niu.^3).*(8+12.*NT.*(1+NT)- 24.*smalldel.*niu.*(1 + 2.*NT) +48.*smalldel.^2.*niu.^2+...
        48.*NT.*niu.*phdelay - 96.*niu.^2.*phdelay.*(smalldel-phdelay)) +...
        1./3.*(phdelay-sgn.*(-1+sgn-2.*NT+4.*niu.*(smalldel-phdelay))./4./niu).^3 +...     
        1/3.*sgn.*((-1+sgn+4.*niu.*phdelay).^3./64./niu.^3-(phdelay-sgn.*(-1+sgn-2.*NT+4.*niu.*(smalldel-phdelay))./4./niu).^3));
    
        bval_matrix=repmat(bval,[1,l_a]);

        % Parallel component
        Lpar=-bval_matrix.*dRes;
     
        
        % Perpendicular component

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

        Lperp = -0.5.*GAMMA^2*GmxSq.*s;
        Lperp(G==0) = 0;
   

    else % the 2nd gradient is the negative of the first one reflected 
        sgn = (-1).^NT;
        bval = G.^2.*GAMMA.^2.*((delta-smalldel).*(sgn.*(smalldel-(0.5.*(1-sgn)+NT)./(2.*niu)-phdelay)-phdelay).^2 +...
            2.*phdelay.^3./3+ NT.*(1-6.*niu.*phdelay+12.*niu.^2.*phdelay.^2)./12./niu.^3+...
            sgn./3.*((phdelay-sgn./4./niu.*(sgn-1)).^3-2.*(phdelay-sgn./4./niu.*(-1+sgn-2.*NT+4.*niu.*(smalldel-phdelay))).^3)+...
            sgn./3.*((sgn-1+4.*niu.*phdelay).^3./64./niu.^3));
        bval_matrix=repmat(bval,[1,l_a]);

        % Parallel component
        Lpar=-bval_matrix.*dRes;
       
        % Perpendicular component

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

        Lperp = -0.5.*GAMMA^2*GmxSq.*s;
        Lperp(G==0) = 0;

        
    end
end

