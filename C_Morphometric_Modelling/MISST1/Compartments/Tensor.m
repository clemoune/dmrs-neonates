function [E J]=Tensor(x,protocol,x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating the Tensor compartment.
% 
% [E,J]=Tensor(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a diffusion tensor substrate and a 
% diffusion protocol specified in the input
% Substrate: Diffusion tensor
% Diffusion pulse sequence: wrapper for various sequences
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% x - size 6 vector of model parameters in SI units for Stick:
%       x(1) - diffusivity of the material along the first direction
%       x(2) - diffusivity of the material along the second direction
%       x(3) - diffusivity of the material along the third direction
%       x(4) - angle from the z direction of the main direction of the
% tensor
%       x(5) - azymuthal angle of the main direction of the
% tensor
%       x(6) - angle of rotation of the second direction 
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
d1 = x(1);
d2 = x(2);
d3 = x(3);
theta = x(4);
phi = x(5);
psi = x(6);


GAMMA = 2.675987E8; % This is what is used throughout Wuzi.
%GAMMA = 2.6751525E8; % This is the latest best estimate of GAMMA (used in Camino)


% calculate main direction from the specified angles
n1 = [sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];

% rotate by 90 degrees in xz plane
n2_zrot = [sin(theta+pi/2)*cos(phi); sin(theta+pi/2)*sin(phi); cos(theta+pi/2)];

n2 = RotateVec(n2_zrot,n1,psi);

n3 = cross(n1,n2);


if strcmp(protocol.pulseseq,'PGSE') || strcmp(protocol.pulseseq,'SWOGSE') || strcmp(protocol.pulseseq,'TWOGSE') || ...
        strcmp(protocol.pulseseq,'DSE') || strcmp(protocol.pulseseq,'STEAM')% sequences with one gradient orientation
    grad_dirs = protocol.grad_dirs;
    % Angles between gradient directions and fibre direction.
    cosTheta1 = grad_dirs*n1;
    cosTheta1Sq = cosTheta1.^2;
    cosTheta2 = grad_dirs*n2;
    cosTheta2Sq = cosTheta2.^2;
    cosTheta3 = grad_dirs*n3;
    cosTheta3Sq = cosTheta3.^2;
    
    Bval = GetBvalues(protocol)';
    B1 = Bval.*cosTheta1Sq;
    B2 = Bval.*cosTheta2Sq;
    B3 = Bval.*cosTheta3Sq;
    e1 = exp(-B1.*d1); 
    e2 = exp(-B2.*d2);
    e3 = exp(-B3.*d3);
    E = e1.*e2.*e3;
elseif strcmp(protocol.pulseseq,'dPGSE')
    b1 = GAMMA.^2.*protocol.G1.^2.*protocol.smalldel.^2.*(protocol.delta-protocol.smalldel/3);    
    b2 = GAMMA.^2.*protocol.G2.^2.*protocol.smalldel.^2.*(protocol.delta-protocol.smalldel/3);
    
    grad_dirs1 = protocol.grad_dirs1;
    grad_dirs2 = protocol.grad_dirs2;
    
    cosTheta11 = grad_dirs1*n1;
    cosTheta21 = grad_dirs2*n1;
    cosTheta11Sq = cosTheta11.^2;
    cosTheta21Sq = cosTheta21.^2;
    
    cosTheta12 = grad_dirs1*n2;
    cosTheta22 = grad_dirs2*n2;
    cosTheta12Sq = cosTheta12.^2;
    cosTheta22Sq = cosTheta22.^2;
    
    cosTheta13 = grad_dirs1*n3;
    cosTheta23 = grad_dirs2*n3;
    cosTheta13Sq = cosTheta13.^2;
    cosTheta23Sq = cosTheta23.^2;
    

    B1 = b1'.*cosTheta11Sq+b2'.*cosTheta21Sq;
    B2 = b1'.*cosTheta12Sq+b2'.*cosTheta22Sq;
    B3 = b1'.*cosTheta13Sq+b2'.*cosTheta23Sq;
    e1 = exp(-B1.*d1); 
    e2 = exp(-B2.*d2);
    e3 = exp(-B3.*d3);
    E = e1.*e2.*e3;
elseif strcmp(protocol.pulseseq,'dSWOGSE') 
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

    
     cosTheta11 = grad_dirs1*n1;
    cosTheta21 = grad_dirs2*n1;
    cosTheta11Sq = cosTheta11.^2;
    cosTheta21Sq = cosTheta21.^2;
    
    cosTheta12 = grad_dirs1*n2;
    cosTheta22 = grad_dirs2*n2;
    cosTheta12Sq = cosTheta12.^2;
    cosTheta22Sq = cosTheta22.^2;
    
    cosTheta13 = grad_dirs1*n3;
    cosTheta23 = grad_dirs2*n3;
    cosTheta13Sq = cosTheta13.^2;
    cosTheta23Sq = cosTheta23.^2;
    

    B1 = b1'.*cosTheta11Sq+b2'.*cosTheta21Sq;
    B2 = b1'.*cosTheta12Sq+b2'.*cosTheta22Sq;
    B3 = b1'.*cosTheta13Sq+b2'.*cosTheta23Sq;
    e1 = exp(-B1.*d1); 
    e2 = exp(-B2.*d2);
    e3 = exp(-B3.*d3);
    E = e1.*e2.*e3;
elseif strcmp(protocol.pulseseq,'isoPGSE') 
    b1 = GAMMA.^2.*2*protocol.G.^2.*(protocol.smalldel./6).^2.*(2*protocol.smalldel/6/3);
    
    grad_dirs1 = protocol.grad_dirs1;
    grad_dirs2 = protocol.grad_dirs2;
    grad_dirs3 = protocol.grad_dirs3;
    
     cosTheta11 = grad_dirs1*n1;
    cosTheta21 = grad_dirs2*n1;
    cosTheta31 = grad_dirs3*n1;
    cosTheta11Sq = cosTheta11.^2;
    cosTheta21Sq = cosTheta21.^2;
    cosTheta31Sq = cosTheta31.^2;
    
    cosTheta12 = grad_dirs1*n2;
    cosTheta22 = grad_dirs2*n2;
    cosTheta32 = grad_dirs3*n2;
    cosTheta12Sq = cosTheta12.^2;
    cosTheta22Sq = cosTheta22.^2;
    cosTheta32Sq = cosTheta32.^2;
    
    cosTheta13 = grad_dirs1*n3;
    cosTheta23 = grad_dirs2*n3;
    cosTheta33 = grad_dirs3*n3;
    cosTheta13Sq = cosTheta13.^2;
    cosTheta23Sq = cosTheta23.^2;
    cosTheta33Sq = cosTheta33.^2;
    

    B1 = b1'.*cosTheta11Sq+b1'.*cosTheta21Sq +b1'.*cosTheta31Sq;
    B2 = b1'.*cosTheta12Sq+b1'.*cosTheta22Sq +b1'.*cosTheta32Sq;
    B3 = b1'.*cosTheta13Sq+b1'.*cosTheta23Sq +b1'.*cosTheta33Sq;
    e1 = exp(-B1.*d1); 
    e2 = exp(-B2.*d2);
    e3 = exp(-B3.*d3);
    E = e1.*e2.*e3;
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
    
    cosTheta11 = grad_dirs1*n1;
    cosTheta21 = grad_dirs2*n1;
    cosTheta11Sq = cosTheta11.^2;
    cosTheta21Sq = cosTheta21.^2;
    
    cosTheta12 = grad_dirs1*n2;
    cosTheta22 = grad_dirs2*n2;
    cosTheta12Sq = cosTheta12.^2;
    cosTheta22Sq = cosTheta22.^2;
    
    cosTheta13 = grad_dirs1*n3;
    cosTheta23 = grad_dirs2*n3;
    cosTheta13Sq = cosTheta13.^2;
    cosTheta23Sq = cosTheta23.^2;
    

    B1 = b1'.*cosTheta11Sq+b2'.*cosTheta21Sq;
    B2 = b1'.*cosTheta12Sq+b2'.*cosTheta22Sq;
    B3 = b1'.*cosTheta13Sq+b2'.*cosTheta23Sq;
    e1 = exp(-B1.*d1); 
    e2 = exp(-B2.*d2);
    e3 = exp(-B3.*d3);
    E = e1.*e2.*e3;
    
else

    tau=protocol.tau;  
    wf=wave_form(protocol);

    % dot product between gradient direction and main direction of the tensor
    F1 = cumsum((wf(:,1:3:end).*n1(1)+wf(:,2:3:end).*n1(2)+wf(:,3:3:end).*n1(3)).*tau,2); 
    % dot product between gradient direction and second direction of the tensor
    F2 = cumsum((wf(:,1:3:end).*n2(1)+wf(:,2:3:end).*n2(2)+wf(:,3:3:end).*n2(3)).*tau,2); 
    % dot product between gradient direction and third direction of the tensor
    F3 = cumsum((wf(:,1:3:end).*n3(1)+wf(:,2:3:end).*n3(2)+wf(:,3:3:end).*n3(3)).*tau,2);

    B1=GAMMA^2*sum((F1.^2)*tau,2);   
    B2=GAMMA^2*sum((F2.^2)*tau,2); 
    B3=GAMMA^2*sum((F3.^2)*tau,2); 


    E=exp(-(B1.*d1 + B2*d2 +B3*d3) );
end

% Compute the Jacobian matrix
if(nargout>1)
   dx = 0.00001;
   J = zeros(length(E), 4);
    if nargin < 3 
         
         dEdd1 = -B1.*E;
         dEdd2 = -B2.*E;
         dEdd3 = -B3.*E;
         J(:,1) = dEdd1;
         J(:,2) = dEdd2;
         J(:,3) = dEdd3;
         for i = 4:6
            xpert = x;
            xpert(i) = xpert(i)*(1+dx);
            Epert = Tensor(xpert, protocol);
            dEdx = (Epert - E)/(xpert(i)*dx);
            J(:,i) = dEdx;
         end         
        
    else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
    
        for i = 1:length(x_deriv);
            if x_deriv(i) ~= 0                  
                xpert = x;
                xpert(i) = xpert(i)*(1+dx);
                Epert = Tensor(xpert, protocol);
                dEdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEdx;
            end
        end
        
       
        
    end   
   
end
end
function vec2 = RotateVec(vec1,axis,alpha) % rotates vector vec1 around axis by alpha
    axis = axis./norm(axis);
    % Rodriguez formula
    vec2 = cos(alpha)*vec1+sin(alpha)*cross(vec1,axis) + (1-cos(alpha))*dot(vec1,axis)*axis;
   
end