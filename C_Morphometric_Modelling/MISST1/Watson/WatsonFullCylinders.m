function [E,J]=WatsonFullCylinders(x, protocol,x_deriv)

% [E,J]=WatsonCylinders(x, protocol,x_deriv)
% returns the measurements E according to the model and the Jacobian J of the
% measurements with respect to the parameters.
%
% Substrate: Impermeable sticks in an empty background.
% Orientation distribution: Watson's distribution 
% Signal approximation: Gaussian phase distribution.
%
% x is the list of model parameters in SI units:
% x(1) is the diffusivity of the material inside the cylinders.
% x(2) is the radius of the cylinders.
% x(3) is the concentration parameter of the Watson's distribution
% x(4) is the angle the polar angle theta in spherical coordinates desbribing the fibre
% direction
% x(5) is the azimuthal angle phi in spherical coordinates describing the
% fibre direction

% protocol includes all the information required to generate the signal
% protocol.grad_dirs is the gradient direction for each measurement.  It has size [N
% 3] where N is the number of measurements.
%
% protocol.G, protocol.delta and protocol.smalldel are the gradient strength, pulse separation and
% pulse length of each measurement in the protocol.  Each has
% size [N 1].
%
%
% protocol.roots_cyl contains solutions to the Bessel function equation from function
% BesselJ_RootsCyl.
%
% $Id$
%
% author: Gary Hui Zhang (gary.zhang@ucl.ac.uk)
% Modified by: Maira Tariq (maira.tariq.11@ucl.ac.uk)
%
% $Id$
% model parameters are: di, rad, kappa, theta, phi
% $Id$

if length(x) ~= 5
    error('the model should have exactly five parameters');
end

d=x(1);
R=x(2);
kappa=x(3);
grad_dirs = protocol.grad_dirs;
l_q = size(protocol.grad_dirs,1);

if ~isfield(protocol,'A_cyl') % sequences that are computed using the GPD approxiation
    
      fun = str2func(['LogComponentsCylinder_' protocol.pulseseq]);
    
  
    % Parallel & perpendicular components of a cylinder
    x_comp = [d R];
    
    [Lperp,Lpar] = fun(x_comp,protocol); % compute perpendicular component (cylinder along z, gradient along x)  
    
    
    ePerp = exp(Lperp);
    
    fibredir = GetFibreOrientation('WatsonFullCylinders', x);
    
    % Compute the Legendre weighted signal
    
    
    Lpmp = Lpar - Lperp;
    E = zeros(size(grad_dirs,1),1);
    watsonNC = WatsonNormalizationCoeff(kappa);
    for i=1:length(E)
        % integrate the contribution from all orientation
        E(i) = dblquad(@(z,phi)WatsonCylNeumanComponent(z, phi, Lpmp(i), kappa, watsonNC, grad_dirs(i,:)*fibredir), -1, 1, 0, pi/2, 1E-5);
    end

    % tag on other factors
    E = (1./pi)*E.*ePerp;
    
    
    % Compute the Jacobian matrix; computed numerically
    if(nargout>1)
        dx = 0.00001;
        J = zeros(length(E), length(x));
        if nargin < 3
            for i = 1:length(x); % compute the derivatives for all model parameters
                xpert = x;
                xpert(i) = xpert(i)*(1+dx);
                Epert = WatsonFullCylinders(xpert, protocol);
                dEdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEdx;
            end
        else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
            
            for i = 1:length(x_deriv);
                if x_deriv(i) ~= 0
                    xpert = x;
                    xpert(i) = xpert(i)*(1+dx);
                    Epert = WatsonFullCylinders(xpert, protocol);
                    dEdx = (Epert - E)/(xpert(i)*dx);
                    J(:,i) = dEdx;
                end
            end
        end
    end
else
    % sequences that use the MM method for the signal calculation
    fun = str2func(['Cylinder_' protocol.pulseseq]);
    Ndir = protocol.Ndir;
    kappa = x(3);
    theta = x(4);
    phi = x(5);
    fibre_dirs = WatsonSamples(Ndir,theta,phi,kappa)';  
    [az incl] = cart2sph(fibre_dirs(1,:),fibre_dirs(2,:),fibre_dirs(3,:)); % azimuth and inclinations

    if isfield(protocol,'smalldel')
        M = size(protocol.smalldel,2);       
    else
         M = size(protocol.G,1);       
    end
    E0 = zeros(M*2,Ndir);
    tmp = protocol.complex;
    protocol.complex = 'complex';
    for i = 1:Ndir
        params = [d R pi/2-incl(i) az(i)];
        E0(:,i) = fun(params,protocol);
    end
    protocol.complex = tmp;
    if strcmp(protocol.complex,'complex')
        E = mean(E0,2);
    elseif strcmp(protocol.complex,'real')
        E = mean(E0(1:M,:),2);
    elseif strcmp(protocol.complex,'abs')
        E = abs(mean(E0(1:M,:)+1i*E0(M+1:end,:),2));
    else error unknown protocol.complex
    end
    
    if(nargout>1)
    J = zeros(length(E),5);
    dx = protocol.pert;
    if nargin < 3          
        for i = 1:5               
            xpert = x;
                if i<=2
                protocol.diff=i;
                else
                protocol.diff = 0; % derivatives with respect to fibre direction
                end
            xpert(i) = xpert(i)*(1+dx);    
            Epert = WatsonFullCylinders(xpert,protocol);
            dEtdx = (Epert - E)/(xpert(i)*dx);
            J(:,i) = dEtdx;
        
        end
    else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian      
        for i = 1:5  
        if x_deriv(i) ~= 0                
                xpert = x;
                    if i<=2
                    protocol.diff=i;
                    else
                    protocol.diff = 0; % derivatives with respect to fibre direction
                    end
                xpert(i) = xpert(i)*(1+dx);    
                Epert = WatsonFullCylinders(xpert,protocol);
                dEtdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEtdx;
        end
        end
                
    end   
    protocol.diff=0;   
    end

    
    
end

