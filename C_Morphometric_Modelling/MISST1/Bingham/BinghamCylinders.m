function [E,J]=BinghamCylinders(x, protocol,x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a BinghamCylinders compartment.
%
% [E,J]=BinghamCylinders(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate of cylinders with a Bingham 
% distribution of orientations and a diffusion protocol specified in the input 
%
% Substrate: Cylinders with a Bingham distribution of orientations
% Diffusion pulse sequence: wrapper for various sequences
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% x - size 7 vector of model parameters in SI units for AstroCylinders:
%       x(1) - free diffusivity of the material inside the cylinders.
%       x(2) - radius of the cylinders.
%       x(3) - kappa: concentration parameter of Bingham distribution
%       x(4) - beta: concentration parameter of Bingham distribution
%       x(5) - theta: polar angle in spherical coordinates desbribing the fibre
% direction
%       x(6) - phi: azimuthal angle in spherical coordinates describing the
% fibre direction
%       x(7) - psi: 3rd euler angle of the Bingham's distribution 
% protocol - structure which includes all the information related to the 
%       diffusion protocol and required to generate the signal. 
%       The basic information for a PGSE sequence is the following:
%       protocol.grad_dirs - is the gradient direction for each measurement.
%           It has size [N 3] where N is the number of measurements.
%       protocol.G - gradient strength, size [1 N]
%       protocol.delta - pulse separation, size [1 N]
%       protocol.smalldel - pulse duration, size [1 N]
%       Other diffusion sequences might have additional fields.
% x_deriv - a vector of 0s and 1s with the same size as x, indicating
%       which parameters are considered in the Jacobian;
%
%------------------------------------------------------------------------
% This file is part of the camino.m toolbox.
% Copyright (c) 2015, UCL Microstructure Imaging Group (MIG), All rights reserved.
% Distributed under the Modified BSD Licence (see: LICENSE.pdf).
%
% Authors:   
%   Gary Hui Zhang (gary.zhang@ucl.ac.uk)
%   Maira Tariq (maira.tariq.11@ucl.ac.uk)
% 
    

if length(x) ~= 7
    error('the model should have exactly 7 parameters');
end

d=x(1);
R=x(2);
kappa=x(3);
beta = x(4);
theta = x(5); phi = x(6); psi = x(7);
fanningdir = computeFanningOrientation(theta, phi, psi);
fibredir = GetFibreOrientation('BinghamCylinders', x);
% the direction normal to the fanning plane
normaldir=cross(fanningdir,fibredir);

grad_dirs = protocol.grad_dirs;

l_q = size(protocol.grad_dirs,1);

if strcmp(protocol.pulseseq,'PGSE') || strcmp(protocol.pulseseq,'OGSE') || strcmp(protocol.pulseseq,'STEAM') ||...
    strcmp(protocol.pulseseq,'SWOGSE') || strcmp(protocol.pulseseq,'TWOGSE')
    
     fun = str2func(['LogComponentsCylinder_' protocol.pulseseq]);
    
  
    % Parallel & perpendicular components of a cylinder
    x_comp = [d R];
    
    [Lperp,Lpar] = fun(x_comp,protocol); % compute perpendicular component (cylinder along z, gradient along x)  
    
    
    ePerp = exp(Lperp);
    
    % Compute the Legendre weighted signal
    
    Lpmp = Lpar - Lperp;
    
       
     E = zeros(size(grad_dirs,1),1);
    BinghamNC = 1/BinghamNormalizationCoeff(kappa,beta);
    for i=1:length(E)
        % Project the gradient directions in the local frame of the Bingham
        % distribution: without loss of generality, we choose
        % the dominant direction is along z, fanning is along y
        matrix = [0 0 0; 0 beta 0; 0 0 kappa];
        projG = [grad_dirs(i,:)*normaldir grad_dirs(i,:)*fanningdir grad_dirs(i,:)*fibredir]';
        matrix = matrix + Lpmp(i)*projG*projG';
        [eigV, eigD] = eig(matrix);
        E(i) = BinghamNormalizationCoeff(eigD(1,1), eigD(2,2), eigD(3,3))*BinghamNC;
    end

    % tag on other factors
    E = E.*ePerp;
    
       
    
    % Compute the Jacobian matrix; computed numerically
    if(nargout>1)
        dx = 0.00001;
        J = zeros(length(E), length(x));
        if nargin < 3
            for i = 1:length(x); % compute the derivatives for all model parameters
                xpert = x;
                xpert(i) = xpert(i)*(1+dx);
                Epert = BinghamCylinders(xpert, protocol);
                dEdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEdx;
            end
        else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
            
            for i = 1:length(x_deriv);
                if x_deriv(i) ~= 0
                    xpert = x;
                    xpert(i) = xpert(i)*(1+dx);
                    Epert = BinghamCylinders(xpert, protocol);
                    dEdx = (Epert - E)/(xpert(i)*dx);
                    J(:,i) = dEdx;
                end
            end
        end
    end
else
    % sequences that use the MM method for the signal calculation
    error('Bingham not yet implemented for sequences with MM')
%     fun = str2func(['Cylinder_' protocol.pulseseq]);
%     Ndir = protocol.Ndir;
%     kappa = x(3);
%     theta = x(4);
%     phi = x(5);
%     fibre_dirs = WatsonSamples(Ndir,theta,phi,kappa)';  
%     [az incl] = cart2sph(fibre_dirs(1,:),fibre_dirs(2,:),fibre_dirs(3,:)); % azimuth and inclinations
% 
%     if isfield(protocol,'smalldel')
%         M = size(protocol.smalldel,2);       
%     else
%          M = size(protocol.G,1);       
%     end
%     E0 = zeros(M*2,Ndir);
%     tmp = protocol.complex;
%     protocol.complex = 'complex';
%     for i = 1:Ndir
%         params = [d R pi/2-incl(i) az(i)];
%         E0(:,i) = fun(params,protocol);
%     end
%     protocol.complex = tmp;
%     if strcmp(protocol.complex,'complex')
%         E = mean(E0,2);
%     elseif strcmp(protocol.complex,'real')
%         E = mean(E0(1:M,:),2);
%     elseif strcmp(protocol.complex,'abs')
%         E = abs(mean(E0(1:M,:)+1i*E0(M+1:end,:),2));
%     else error unknown protocol.complex
%     end
%     
%     if(nargout>1)
%     J = zeros(length(E),5);
%     dx = protocol.pert;
%     if nargin < 3          
%         for i = 1:5               
%             xpert = x;
%                 if i<=2
%                 protocol.diff=i;
%                 else
%                 protocol.diff = 0; % derivatives with respect to fibre direction
%                 end
%             xpert(i) = xpert(i)*(1+dx);    
%             Epert = WatsonCylinders(xpert,protocol);
%             dEtdx = (Epert - E)/(xpert(i)*dx);
%             J(:,i) = dEtdx;
%         
%         end
%     else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian      
%         for i = 1:5  
%         if x_deriv(i) ~= 0                
%                 xpert = x;
%                     if i<=2
%                     protocol.diff=i;
%                     else
%                     protocol.diff = 0; % derivatives with respect to fibre direction
%                     end
%                 xpert(i) = xpert(i)*(1+dx);    
%                 Epert = WatsonCylinders(xpert,protocol);
%                 dEtdx = (Epert - E)/(xpert(i)*dx);
%                 J(:,i) = dEtdx;
%         end
%         end
%                 
%     end   
%     protocol.diff=0;   
%     end

    
    
end

