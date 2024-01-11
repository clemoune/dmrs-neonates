function [E J]=Zeppelin(x,protocol,x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating the Zeppelin compartment.
% 
% [E,J]=Zeppelin(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a cyllindrically symmetric diffusion tensor 
% and a diffusion protocol specified in the input
% Substrate: Cyllinderically symmetric diffusion tensor (hindered diffusion)
% Diffusion pulse sequence: wrapper for various sequences
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% x - size 4 vector of model parameters in SI units for Zeppelin:
%       x(1) - parallel diffusivity
%       x(2) - perpendicular diffusivity (hindered)
%       x(3) - polar angle theta in spherical coordinates desbribing the
%       principal direction
%       x(4) - azimuthal angle phi in spherical coordinates describing the
% principal direction
% protocol - structure which includes all the information related to the 
%        diffusion protocol and required to generate the signal. 
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

% Model parameters
if iscell(x) % tis will be the case for models that vary diffusivity with frequency
    dPar=cell2mat(x(1));
    dPerp=cell2mat(x(2)); 
    theta = cell2mat(x(3));
    phi = cell2mat(x(4));
else
    dPar = x(1);
    dPerp = x(2);
    theta = x(3);
    phi = x(4);
end

% calculate fibre direction from the specified angles
fibredir = [sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];



GAMMA = 2.675987E8;
if strcmp(protocol.pulseseq,'PGSE') || strcmp(protocol.pulseseq,'SWOGSE') || strcmp(protocol.pulseseq,'TWOGSE') || ...
        strcmp(protocol.pulseseq,'DSE') || strcmp(protocol.pulseseq,'STEAM')% sequences with one gradient orientation
        grad_dirs = protocol.grad_dirs;
        % Angles between gradient directions and fibre direction.
        cosTheta = grad_dirs*fibredir;
        cosThetaSq = cosTheta.^2;

        Bval = GetBvalues(protocol)';
        Bval_par = Bval.*cosThetaSq;
        Bval_perp = Bval.*(1-cosThetaSq);
        ePar=exp(-Bval_par.*dPar); 
        ePerp = exp(-Bval_perp.*dPerp);
        E = ePar.*ePerp;
elseif strcmp(protocol.pulseseq,'FullSTEAM')
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
        % Angles between gradient directions and fibre direction.
        cosTheta = grad_dirs*fibredir;
        % Components parallel to fibre
        GD_CompPar = G.*cosTheta;
        GC_CompPar = cG*fibredir;
        GR_CompPar = rG*fibredir;
        
        fibdirrep = repmat(fibredir',[size(grad_dirs,1) 1]);
        comp1_dirs = grad_dirs - repmat(sum(grad_dirs.*fibdirrep,2),[1 3]).*fibdirrep;
        % Find any measurements where the gradient and fibre directions
        % are parallel
        fgpar = find(abs(cosTheta) == 1);
        comp1_dirs(fgpar, :) = repmat([fibredir(3) fibredir(3) fibredir(1)+fibredir(2)], [length(fgpar) 1]);
        comp1_dirs = comp1_dirs./repmat(sqrt(sum(comp1_dirs.^2, 2)), [1 3]);
        GD_Comp1 = G.*sqrt(1-cosTheta.^2);
        GC_Comp1 = sum(comp1_dirs.*cG,2);
        GR_Comp1 = sum(comp1_dirs.*rG,2);

        % Components along a perpendicular direction also perpendicular
        % to the diffusion gradient.
        comp2_dirs = cross(grad_dirs, fibdirrep);
        comp2_dirs(fgpar, :) = cross(grad_dirs(fgpar,:), comp1_dirs(fgpar,:));
        comp2_dirs = comp2_dirs./repmat(sqrt(sum(comp2_dirs.^2, 2)), [1 3]);
        GC_Comp2 = sum(comp2_dirs.*cG,2);
        GR_Comp2 = sum(comp2_dirs.*rG,2);
        GD_Comp2 = zeros(size(GC_Comp2));
        
        tdd = gap1 + gap2 + TM + 2*sdelc + 2*smalldel/3 + 2*sdelr;
        tcc = TM + 2*sdelc/3 + 2*sdelr;
        trr = TM + 2*sdelr/3;
        tdc = TM + sdelc + 2*sdelr;
        tdr = TM + sdelr;
        tcr = TM + sdelr;
        
        bvalpar = GAMMA^2*(GD_CompPar.^2.*smalldel.^2.*tdd + GC_CompPar.^2.*sdelc.^2.*tcc + ...
    GR_CompPar.^2.*sdelr.^2.*trr + 2*GD_CompPar.*GC_CompPar.*smalldel.*sdelc.*tdc + ...
    2*GD_CompPar.*GR_CompPar.*smalldel.*sdelr.*tdr + 2*GC_CompPar.*GR_CompPar.*sdelc.*sdelr.*tcr);
    
        bvalperp1 = GAMMA^2*(GD_Comp1.^2.*smalldel.^2.*tdd + GC_Comp1.^2.*sdelc.^2.*tcc + ...
    GR_Comp1.^2.*sdelr.^2.*trr + 2*GD_Comp1.*GC_Comp1.*smalldel.*sdelc.*tdc + ...
    2*GD_Comp1.*GR_Comp1.*smalldel.*sdelr.*tdr + 2*GC_Comp1.*GR_Comp1.*sdelc.*sdelr.*tcr);

        bvalperp2 = GAMMA^2*(GD_Comp2.^2.*smalldel.^2.*tdd + GC_Comp2.^2.*sdelc.^2.*tcc + ...
    GR_Comp2.^2.*sdelr.^2.*trr + 2*GD_Comp2.*GC_Comp2.*smalldel.*sdelc.*tdc + ...
    2*GD_Comp2.*GR_Comp2.*smalldel.*sdelr.*tdr + 2*GC_Comp2.*GR_Comp2.*sdelc.*sdelr.*tcr);        

E = exp(-bvalpar.*dPar).*exp(-bvalperp1.*dPerp).*exp(-bvalperp2.*dPerp);
elseif strcmp(protocol.pulseseq,'dPGSE')
        b1 = GAMMA.^2.*protocol.G1.^2.*protocol.smalldel.^2.*(protocol.delta-protocol.smalldel/3);    
        b2 = GAMMA.^2.*protocol.G2.^2.*protocol.smalldel.^2.*(protocol.delta-protocol.smalldel/3);
        grad_dirs1 = protocol.grad_dirs1;
        grad_dirs2 = protocol.grad_dirs2;

        cosTheta1 = grad_dirs1*fibredir;
        cosThetaSq1 = cosTheta1.^2;

        cosTheta2 = grad_dirs2*fibredir;
        cosThetaSq2 = cosTheta2.^2;
        
        Bval_par = b1'.*cosThetaSq1+b2'.*cosThetaSq2;
        Bval_perp = b1'.*(1-cosThetaSq1)+b2'.*(1-cosThetaSq2);
        
        ePar=exp(-Bval_par.*dPar); 
        ePerp = exp(-Bval_perp.*dPerp);
        E = ePar.*ePerp;
elseif  strcmp(protocol.pulseseq,'dSWOGSE') 
    G1 = protocol.G1;
    G2 = protocol.G2;
    omega = protocol.omega;    
    niu = omega/(2*pi());
    smalldel = protocol.smalldel;
    delta = protocol.delta;
    
    if ~isfield(protocol,'mirror') || protocol.mirror == 0
        if ~isfield(protocol,'phase') 
            NT = floor(2.*smalldel.*niu+0.00000000001);
            b1 = GAMMA.^2.*G1.^2.*(1./(48.*niu.^3)).*((2.*NT.^3 +... 
             3.*NT.^2.*(1 + (-1).^NT - 4.*niu.*smalldel) - 4.*niu.^2.*smalldel.^2.*...
             (-3 + 3.*(-1).^NT + 4.*niu.*smalldel) +3.*delta.*niu.*(-1 + (-1).^NT - ...
             2.*NT + 4.*niu.*smalldel).^2 + NT.*(1 + 3.*(-1).^NT - 12.*niu.*smalldel + ...
               24.*niu.^2.*smalldel.^2)));
            b2 = GAMMA.^2.*G2.^2.*(1./(48.*niu.^3)).*((2.*NT.^3 +... 
             3.*NT.^2.*(1 + (-1).^NT - 4.*niu.*smalldel) - 4.*niu.^2.*smalldel.^2.*...
             (-3 + 3.*(-1).^NT + 4.*niu.*smalldel) +3.*delta.*niu.*(-1 + (-1).^NT - ...
             2.*NT + 4.*niu.*smalldel).^2 + NT.*(1 + 3.*(-1).^NT - 12.*niu.*smalldel + ...
               24.*niu.^2.*smalldel.^2)));
        else
             phase = protocol.phase;
             for i = 1:length(omega)
                if omega(i)<pi/smalldel(i);
                    omega(i) = pi/smalldel(i);
                    phase(i) = 0;
                end
             end
            phase = mod(phase,2*pi);
            phase(phase>pi) = phase(phase>pi)-2*pi;
            phase(phase<0) = pi-abs(phase(phase<0));
            phdelay = phase ./(2 *pi()* niu);

            NT = floor(2.*(smalldel-phdelay).*niu+0.00000000001);
            sgn = (-1).^NT;

            b1 = GAMMA.^2.*G1.^2.*(1./3.*(smalldel-NT./2./niu-phdelay).^3 + (delta-smalldel).*...
                (sgn.*(smalldel- (0.5.*(1-sgn)+NT)./2./niu-phdelay)-phdelay).^2  +(phdelay.^3)./3 +...
                sgn.*((-1+sgn-2.*NT+4.*niu.*(smalldel-phdelay)).^3)./192./(niu.^3) +...   
                 NT./96./(niu.^3).*(8+12.*NT.*(1+NT)- 24.*smalldel.*niu.*(1 + 2.*NT) +48.*smalldel.^2.*niu.^2+...
                48.*NT.*niu.*phdelay - 96.*niu.^2.*phdelay.*(smalldel-phdelay)) +...
                1./3.*(phdelay-sgn.*(-1+sgn-2.*NT+4.*niu.*(smalldel-phdelay))./4./niu).^3 +...     
                1/3.*sgn.*((-1+sgn+4.*niu.*phdelay).^3./64./niu.^3-(phdelay-sgn.*(-1+sgn-2.*NT+4.*niu.*(smalldel-phdelay))./4./niu).^3));
            b2 = GAMMA.^2.*G2.^2.*(1./3.*(smalldel-NT./2./niu-phdelay).^3 + (delta-smalldel).*...
                (sgn.*(smalldel- (0.5.*(1-sgn)+NT)./2./niu-phdelay)-phdelay).^2  +(phdelay.^3)./3 +...
                sgn.*((-1+sgn-2.*NT+4.*niu.*(smalldel-phdelay)).^3)./192./(niu.^3) +...   
                 NT./96./(niu.^3).*(8+12.*NT.*(1+NT)- 24.*smalldel.*niu.*(1 + 2.*NT) +48.*smalldel.^2.*niu.^2+...
                48.*NT.*niu.*phdelay - 96.*niu.^2.*phdelay.*(smalldel-phdelay)) +...
                1./3.*(phdelay-sgn.*(-1+sgn-2.*NT+4.*niu.*(smalldel-phdelay))./4./niu).^3 +...     
                1/3.*sgn.*((-1+sgn+4.*niu.*phdelay).^3./64./niu.^3-(phdelay-sgn.*(-1+sgn-2.*NT+4.*niu.*(smalldel-phdelay))./4./niu).^3));

        end
    else
         if ~isfield(protocol,'phase') 
             NT = floor(2.*smalldel.*niu+0.00000000001);

            b1 = GAMMA.^2.*G1.^2.*((delta-smalldel).*(smalldel-(0.5.*(1-(-1).^NT)+NT)./(2.*niu)).^2 + NT./(12.*niu.^3) +...
                1./(192.*niu.^3).*(-(-1+(-1).^NT).^3+(-1+(-1).^NT-2.*NT+4.*smalldel.*niu).^3) - ...
                1./(96.*niu.^3).*(NT-2.*smalldel.*niu).*(3.*(-1+(-1).^NT).^2+4.*NT.^2+12.*smalldel.*niu.*(-1+(-1).^NT)+...
                16.*smalldel.^2.*niu.^2-2.*NT.*(-3+3.*(-1).^NT+8.*smalldel.*niu)));
            b2 = GAMMA.^2.*G2.^2.*((delta-smalldel).*(smalldel-(0.5.*(1-(-1).^NT)+NT)./(2.*niu)).^2 + NT./(12.*niu.^3) +...
                1./(192.*niu.^3).*(-(-1+(-1).^NT).^3+(-1+(-1).^NT-2.*NT+4.*smalldel.*niu).^3) - ...
                1./(96.*niu.^3).*(NT-2.*smalldel.*niu).*(3.*(-1+(-1).^NT).^2+4.*NT.^2+12.*smalldel.*niu.*(-1+(-1).^NT)+...
                16.*smalldel.^2.*niu.^2-2.*NT.*(-3+3.*(-1).^NT+8.*smalldel.*niu)));
         else
              phase = protocol.phase;
             for i = 1:length(omega)
                if omega(i)<pi/smalldel(i);
                    omega(i) = pi/smalldel(i);
                    phase(i) = 0;
                end
             end
            phase = mod(phase,2*pi);
            phase(phase>pi) = phase(phase>pi)-2*pi;
            phase(phase<0) = pi-abs(phase(phase<0));
            phdelay = phase ./(2 *pi()* niu);

            NT = floor(2.*(smalldel-phdelay).*niu+0.00000000001);
            sgn = (-1).^NT;

            b1 = G1.^2.*GAMMA.^2.*((delta-smalldel).*(sgn.*(smalldel-(0.5.*(1-sgn)+NT)./(2.*niu)-phdelay)-phdelay).^2 +...
                2.*phdelay.^3./3+ NT.*(1-6.*niu.*phdelay+12.*niu.^2.*phdelay.^2)./12./niu.^3+...
                sgn./3.*((phdelay-sgn./4./niu.*(sgn-1)).^3-2.*(phdelay-sgn./4./niu.*(-1+sgn-2.*NT+4.*niu.*(smalldel-phdelay))).^3)+...
                sgn./3.*((sgn-1+4.*niu.*phdelay).^3./64./niu.^3));
            b2 = G2.^2.*GAMMA.^2.*((delta-smalldel).*(sgn.*(smalldel-(0.5.*(1-sgn)+NT)./(2.*niu)-phdelay)-phdelay).^2 +...
                2.*phdelay.^3./3+ NT.*(1-6.*niu.*phdelay+12.*niu.^2.*phdelay.^2)./12./niu.^3+...
                sgn./3.*((phdelay-sgn./4./niu.*(sgn-1)).^3-2.*(phdelay-sgn./4./niu.*(-1+sgn-2.*NT+4.*niu.*(smalldel-phdelay))).^3)+...
                sgn./3.*((sgn-1+4.*niu.*phdelay).^3./64./niu.^3));
         end
    end
    grad_dirs1 = protocol.grad_dirs1;
    grad_dirs2 = protocol.grad_dirs2;
    
    cosTheta1 = grad_dirs1*fibredir;
    cosThetaSq1 = cosTheta1.^2;
    
    cosTheta2 = grad_dirs2*fibredir;
    cosThetaSq2 = cosTheta2.^2;
    
    Bval_par = b1'.*cosThetaSq1+b2'.*cosThetaSq2;
    Bval_perp = b1'.*(1-cosThetaSq1)+b2'.*(1-cosThetaSq2);
    ePar=exp(-Bval_par.*dPar);   
    ePerp=exp(-Bval_perp.*dPerp);   
    E = ePar.*ePerp;
elseif strcmp(protocol.pulseseq,'isoPGSE') 
    b1 = GAMMA.^2.*2*protocol.G.^2.*(protocol.smalldel./6).^2.*(2*protocol.smalldel/6/3);
    
    grad_dirs1 = protocol.grad_dirs1;
    grad_dirs2 = protocol.grad_dirs2;
    grad_dirs3 = protocol.grad_dirs3;
    
    cosTheta1 = grad_dirs1*fibredir;
    cosThetaSq1 = cosTheta1.^2;
    
    cosTheta2 = grad_dirs2*fibredir;
    cosThetaSq2 = cosTheta2.^2;
    
    cosTheta3 = grad_dirs3*fibredir;
    cosThetaSq3 = cosTheta3.^2;
    
    Bval_par = b1'.*cosThetaSq1+b1'.*cosThetaSq2+b1'.*cosThetaSq3;
    Bval_perp = b1'.*(1-cosThetaSq1)+b1'.*(1-cosThetaSq2)+b1'.*(1-cosThetaSq3);    
    ePar=exp(-Bval_par.*dPar);   
    ePerp=exp(-Bval_perp.*dPerp);   
    E = ePar.*ePerp;
elseif strcmp(protocol.pulseseq,'DODE')

     if ~isfield(protocol,'phase') || max(protocol.phase) == 0      
        b1 = GAMMA.^2.*protocol.G1.^2./2.*protocol.smalldel.^3./6./protocol.Nosc.^2;
         b2 = GAMMA.^2.*protocol.G2.^2./2.*protocol.smalldel.^3./6./protocol.Nosc.^2;
     elseif max(protocol.phase) == pi/2 || max(protocol.phase) == -pi/2   
         b1 = GAMMA.^2.*protocol.G1.^2./2.*protocol.smalldel.^3./24./protocol.Nosc.^2;
         b2 = GAMMA.^2.*protocol.G2.^2./2.*protocol.smalldel.^3./24./protocol.Nosc.^2;
         
     else error('DODE defined only for phase= 0 or pi/2')
     end    
        grad_dirs1 = protocol.grad_dirs1;
        grad_dirs2 = protocol.grad_dirs2;

        cosTheta1 = grad_dirs1*fibredir;
        cosThetaSq1 = cosTheta1.^2;

        cosTheta2 = grad_dirs2*fibredir;
        cosThetaSq2 = cosTheta2.^2;
        
        Bval_par = b1'.*cosThetaSq1+b2'.*cosThetaSq2;
        Bval_perp = b1'.*(1-cosThetaSq1)+b2'.*(1-cosThetaSq2);
        
        ePar=exp(-Bval_par.*dPar); 
        ePerp = exp(-Bval_perp.*dPerp);
        E = ePar.*ePerp;
else

    tau=protocol.tau;
    wf = wave_form(protocol);
    G_dot_fibre = wf(:,1:3:end)*fibredir(1)+ wf(:,2:3:end)*fibredir(2) + wf(:,3:3:end)*fibredir(3); 
    % x and y directions in he system of the principal axes
     v = [cos(theta)*cos(phi)^2+sin(phi)^2; -(1-cos(theta))*sin(phi)*cos(phi); -sin(theta)*cos(phi)]; % vectors for the new x and y directions; see Ozarslan 2010
     w = [ -(1-cos(theta))*sin(phi)*cos(phi); cos(theta)*sin(phi)^2+cos(phi)^2; -sin(theta)*sin(phi)];
     G_dot_x =  wf(:,1:3:end)*v(1)+ wf(:,2:3:end)*v(2) +wf(:,3:3:end)*v(3) ;  
     G_dot_y =  wf(:,1:3:end)*w(1)+ wf(:,2:3:end)*w(2) +wf(:,3:3:end)*w(3) ;  
    Fpar = cumsum(G_dot_fibre.*tau,2);
    Fx =  cumsum(G_dot_x.*tau,2);
    Fy =  cumsum(G_dot_y.*tau,2);    
   
    Bval_par=sum(Fpar.^2*tau,2);      
    Bval_par=GAMMA^2*Bval_par;     
    Bval_perp=sum((Fx.^2+Fy.^2)*tau,2);
    Bval_perp=GAMMA^2*Bval_perp;
    ePar=exp(-Bval_par.*dPar);  
    ePerp=exp(-Bval_perp.*dPerp); 
    E = ePar.*ePerp;
end
     

% Compute the Jacobian matrix
if(nargout>1)
   dx = 0.00001;
   J = zeros(length(E), 4);
    if nargin < 3 
         
         dEddPar = -Bval_par.*ePar.*ePerp;
         J(:,1) = dEddPar; 
         dEddPerp = -Bval_perp.*ePar.*ePerp;
         J(:,2) = dEddPerp;
         for i = 3:4
            xpert = x;
            xpert(i) = xpert(i)*(1+dx);
            Epert = Zeppelin(xpert, protocol);
            dEdx = (Epert - E)/(xpert(i)*dx);
            J(:,i) = dEdx;
         end         
        
    else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
       
        for i = 1:length(x_deriv);
            if x_deriv(i) ~= 0  
                
                xpert = x;
                xpert(i) = xpert(i)*(1+dx);
                Epert = Zeppelin(xpert, protocol);
                dEdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEdx;
            end
        end
        
       
        
    end   
   
end