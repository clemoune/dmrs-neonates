function [E,J]=WatsonFiniteCylinders_oldJac(x, protocol,x_deriv)

% [E,J]=WatsonCylinders(x, protocol,x_deriv)
% returns the measurements E according to the model and the Jacobian J of the
% measurements with respect to the parameters.
%
% Substrate: Impermeable sticks in an empty background.
% Orientation distribution: Watson's distribution with SH approximation
% Signal approximation: Gaussian phase distribution.
%
% x is the list of model parameters in SI units:
% x(1) is the diffusivity of the material inside the cylinders.
% x(2) is the radius of the cylinders.
% x(3) is the length of the cylinders
% x(4) is the concentration parameter of the Watson's distribution
% x(5) is the angle the polar angle theta in spherical coordinates desbribing the fibre
% direction
% x(6) is the azimuthal angle phi in spherical coordinates describing the
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
% Modified by: Andrada Ianus (a.ianus.11@ucl.ac.uk)
%
% $Id$
% model parameters are: di, rad, kappa, theta, phi
% $Id$

if length(x) ~= 6
    error('the model should have exactly five parameters');
end

d=x(1);
R=x(2);
lcyl = x(3);
kappa=x(4);



if isfield(protocol,'approx')
    approx = protocol.approx;
else
    approx = 'GPD';
end

% if ~isfield(protocol,'numint')
%     numint = 0;
% else
%     if protocol.numint == 0
%         numint = 0;
%     else
%         numint =1;
%     end
% end


if strcmp(protocol.pulseseq,'PGSE') || strcmp(protocol.pulseseq,'SWOGSE') || strcmp(protocol.pulseseq,'TWOGSE') || ...
        strcmp(protocol.pulseseq,'DSE') || strcmp(protocol.pulseseq,'STEAM')% sequences with 1 gradient direction
    
    l_q = size(protocol.grad_dirs,1);
    fun = str2func(['FiniteCylinder_' approx '_' protocol.pulseseq]);
    
    protocol1 = protocol; %  make protocol with gradient gradient oriented along x direction
    protocol1.grad_dirs = repmat([1 0 0],length(protocol.smalldel),1);
    
    % Parallel & perpendicular components of a cylinder
    x_comp = [d R lcyl];
    
    Lperp = log(fun([x_comp 0 0],protocol1)); % compute perpendicular component (cylinder along z, gradient along x)
    Lpar = log(fun([x_comp pi/2 0],protocol1)); % compute parallel component (cylinder along x, gradient along x)
    
    
    ePerp = exp(Lperp);
    
    % Compute the Legendre weighted signal
    
    Lpmp = Lperp - Lpar;
    
    if nargout > 1
        [lgi, J_lgi] = LegendreGaussianIntegral(Lpmp, 6);
    else
        lgi = LegendreGaussianIntegral(Lpmp, 6);
    end
    
    % Compute the spherical harmonic coefficients of the Watson's distribution
    if nargout > 1
        [coeff, J_coeff] = WatsonSHCoeff(kappa);
    else
        coeff = WatsonSHCoeff(kappa);
    end
    coeffMatrix = repmat(coeff, [l_q, 1]);
    
    % Compute the dot product between the symmetry axis of the Watson's distribution
    % and the gradient direction
    %
    % For numerical reasons, cosTheta might not always be between -1 and 1
    % Due to round off errors, individual gradient vectors in grad_dirs and the
    % fibredir are never exactly normal.  When a gradient vector and fibredir are
    % essentially parallel, their dot product can fall outside of -1 and 1.
    %
    % BUT we need make sure it does, otherwise the legendre function call below
    % will FAIL and abort the calculation!!!
    %
    
    fibredir = GetFibreOrientation('WatsonFiniteCylinders', x);
    
    cosTheta = protocol.grad_dirs*fibredir;
    badCosTheta = find(abs(cosTheta)>1);
    cosTheta(badCosTheta) = cosTheta(badCosTheta)./abs(cosTheta(badCosTheta));
    
    % Compute the SH values at cosTheta
    sh = zeros(size(coeff));
    shMatrix = repmat(sh, [l_q, 1]);
    for i = 1:7
        shMatrix(:,i) = sqrt((i - .75)/pi);
        % legendre function returns coefficients of all m from 0 to l
        % we only need the coefficient corresponding to m = 0
        % WARNING: make sure to input ROW vector as variables!!!
        % cosTheta is expected to be a COLUMN vector.
        tmp = legendre(2*i - 2, cosTheta');
        tmp = tmp';
        shMatrix(:,i) = shMatrix(:,i) .* tmp(:,1);
    end
    
    E = sum(lgi.*coeffMatrix.*shMatrix, 2);
    % with the SH approximation, there will be no guarantee that E will be positive
    % but we need to make sure it does!!! replace the negative values with 10% of
    % the smallest positive values
    E(find(E<=0)) = min(E(find(E>0)))*0.1;
    E = 0.5*E.*ePerp;
    
    
    % Compute the Jacobian matrix; computed numerically
    if(nargout>1)
        dx = 0.00001;
        J = zeros(length(E), length(x));
        if nargin < 3
            for i = 1:length(x); % compute the derivatives for all model parameters
                xpert = x;
                xpert(i) = xpert(i)*(1+dx);
                Epert = WatsonFiniteCylinders_oldJac(xpert, protocol);
                dEdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEdx;
            end
        else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
            
            for i = 1:length(x_deriv);
                if x_deriv(i) ~= 0
                    xpert = x;
                    xpert(i) = xpert(i)*(1+dx);
                    Epert = WatsonFiniteCylinders_oldJac(xpert, protocol);
                    dEdx = (Epert - E)/(xpert(i)*dx);
                    J(:,i) = dEdx;
                end
            end
        end
    end
% elseif strcmp(protocol.pulseseq,'dPGSE') && numint ~= 0
%      fun = str2func(['FiniteCylinderComponents_' protocol.pulseseq]);
%      kappa = x(4);   
%      
%      maindir = GetFibreOrientation('WatsonFiniteCylinders', x);
%      % change coordinate system: maindirection along z
%      v = cross(maindir,[0 0 1]');
%      s = norm(v);
%      if s < 1E-4
%          Rot = eye(3);
%      else
%         c = dot(maindir,[0 0 1]');
%         vmat = [0 - v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0 ];
%         Rot = eye(3)+vmat+vmat^2*(1-c)/s^2;
%      end     
%      
%      E = zeros(size(protocol.grad_dirs1,1),1);
%      watsonNC = WatsonNormalizationCoeff(kappa);
%      [C, A11, A22, A12] = fun(x,protocol);
%      for i = 1:length(E)
%          trans_gd1 = Rot*protocol.grad_dirs1(i,:)';
%          trans_gd2 = Rot*protocol.grad_dirs2(i,:)';    
%           E(i) = exp(C(i)).*dblquad(@(theta,phi)WatsonComponent_dPGSE(theta, phi, A11(i), A22(i),A12(i), kappa, watsonNC, trans_gd1,trans_gd2), 0, pi, 0, 2*pi, 1E-5);
%      end

else
    
  
    % sequences that use the MM method for the signal calculation
    fun = str2func(['FiniteCylinder_' approx '_' protocol.pulseseq]);
    if isfield(protocol,'Ndir')
        Ndir = protocol.Ndir;
    else
        Ndir = 100;
    end
    
    if isfield(protocol,'smalldel')
        M = size(protocol.smalldel,2);       
    elseif isfield(protocol,'smalldel1')
        M = size(protocol.smalldel1,2);   
    else
         M = size(protocol.G,1);       
    end
    kappa = x(4);
    theta = x(5);
    phi = x(6);
%     fibre_dirs = WatsonSamples(Ndir,theta,phi,kappa)';  
    fibre_dirs = ReadCaminoElecPS(sprintf('PointSets/Elec%03i.txt',Ndir))';  
    [az incl] = cart2sph(fibre_dirs(1,:),fibre_dirs(2,:),fibre_dirs(3,:)); % azimuth and inclinations

    maindir = GetFibreOrientation('WatsonFiniteCylinders', x); 

    WatsonWeights = hypergeom(1/2,3/2,kappa)^(-1).*exp(kappa.*(maindir'*fibre_dirs).^2);
    WatsonWeights = repmat(WatsonWeights,M,1);
    
   

        if strcmp(approx,'MM')
            E0 = zeros(M*2,Ndir);
            tmp = protocol.complex;
            protocol.complex = 'complex';
            for i = 1:Ndir
                params = [d R lcyl pi/2-incl(i) az(i)];
                E0(:,i) = fun(params,protocol);
            end
            protocol.complex = tmp;
            if strcmp(protocol.complex,'complex')
                E = mean(E0.*repmat(WatsonWeights,2,1),2);
            elseif strcmp(protocol.complex,'real')
                E = mean(E0(1:M,:).*WatsonWeights,2);
            elseif strcmp(protocol.complex,'abs')
                E = abs(mean(E0(1:M,:).*WatsonWeights+1i*E0(M+1:end,:).*WatsonWeights,2));
            else error unknown protocol.complex
            end
        else
          E0 = zeros(M,Ndir);
            for i = 1:Ndir
                params = [d R lcyl pi/2-incl(i) az(i)];
                E0(:,i) = fun(params,protocol);
            end
            E = mean(E0.*WatsonWeights,2);
        end
    
   if nargout >1 % calculate jacobians
        
        J = zeros(length(E),6);
        
        if isfield(protocol,'pert') 
            dx = protocol.pert;
        else
            dx = 1E-5;
        end      
    
        if nargin < 3               
        
        
        for i = 1:5               
            xpert = x;
                if i<=2
                protocol.diff=i;
                else
                protocol.diff = 0; % derivatives with respect to fibre direction
                end
            xpert(i) = xpert(i)*(1+dx);    
            Epert = WatsonFiniteCylinders_oldJac(xpert,protocol);
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
                Epert = WatsonFiniteCylinders_oldJac(xpert,protocol);
                dEtdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEtdx;
        end
        end
                
    end   
    protocol.diff=0;   
     end    
    
end

