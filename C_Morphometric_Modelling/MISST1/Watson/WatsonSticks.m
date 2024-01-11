function [E,J]=WatsonSticks(x, protocol,x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a WatsonSticks compartment.
%
% [E,J]=WatsonSticks(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate of sticks with a Watson 
% distribution of orientations and a diffusion protocol specified in the input 
%
% Substrate: Sticks with a Watson distribution of orientations; there is
% not diffusion perpendicular to the axis of each stick
% Diffusion pulse sequence: wrapper for various sequences
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% x - size 4 vector of model parameters in SI units for AstroCylinders:
%       x(1) - free diffusivity of the material along the sticks.
%       x(2) - concentration parameter of Watson distribution
%       x(3) is the angle the polar angle theta in spherical coordinates desbribing the fibre
% direction
%       x(4) is the azimuthal angle phi in spherical coordinates describing the
% fibre direction
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
%   Andrada Ianus (a.ianus.11@ucl.ac.uk)
% 

if length(x) ~= 4
    error('the model should have exactly four parameters');
end

d=x(1);
kappa=x(2);


    if  strcmp(protocol.pulseseq,'PGSE') || strcmp(protocol.pulseseq,'SWOGSE') || strcmp(protocol.pulseseq,'TWOGSE') || ...
            strcmp(protocol.pulseseq,'DSE') || strcmp(protocol.pulseseq,'STEAM')% sequences with 1 gradient direction

            fun = str2func('Stick');
            
            l_q = size(protocol.grad_dirs,1);
    
            protocol1 = protocol; %  make protocol with gradient gradient oriented along x direction
            protocol1.grad_dirs = repmat([1 0 0],length(protocol.smalldel),1);

            % Parallel components of the stick
            x_comp = [d];

            Lpar = log(fun([x_comp pi/2 0],protocol1)); % compute parallel component (stick along x, gradient along x)

            % Since R=0, no need to do any calculation
            Lperp = zeros(size(protocol.G,1), 1);

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

            fibredir = GetFibreOrientation('WatsonSticks', x);

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
    else
        % sequences that use the MM method for the signal calculation
        fun = str2func('Stick');    
        Ndir = 100;
   

        fibre_dirs = ReadCaminoElecPS(sprintf('PointSets/Elec%03i.txt',Ndir))';  
        [az incl] = cart2sph(fibre_dirs(1,:),fibre_dirs(2,:),fibre_dirs(3,:)); % azimuth and inclinations

        maindir = GetFibreOrientation('WatsonSticks', x);

        WatsonWeights = hypergeom(1/2,3/2,kappa)^(-1).*exp(kappa.*(maindir'*fibre_dirs).^2);   
        WatsonWeights = WatsonWeights./sum(WatsonWeights);
        theta = pi/2-incl;
        phi = az;
        Et = [];
        for i = 1:Ndir
            xnew = [x(1) theta(i) phi(i)];
            Ei = fun(xnew,protocol);
            Et = [Et Ei];
        end
        E = sum(Et.*repmat(WatsonWeights,size(Et,1),1),2);
        
    end
    
    
    % Compute the Jacobian matrix; computed numerically
    if(nargout>1)
        dx = 0.00001;
        J = zeros(length(E), length(x));
        if nargin < 3
            for i = 1:length(x); % compute the derivatives for all model parameters
                xpert = x;
                xpert(i) = xpert(i)*(1+dx);
                Epert = WatsonSticks(xpert, protocol);
                dEdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEdx;
            end
        else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
            
            for i = 1:length(x_deriv);
                if x_deriv(i) ~= 0
                    xpert = x;
                    xpert(i) = xpert(i)*(1+dx);
                    Epert = WatsonSticks(xpert, protocol);
                    dEdx = (Epert - E)/(xpert(i)*dx);
                    J(:,i) = dEdx;
                end
            end
        end
    end














% 
% d=x(1);
% kappa=x(2);
% 
% l_q = size(protocol.grad_dirs,1);
% 
% if ~isfield(protocol,'A_cyl') % sequences that are computed using the GPD approxiation
%     
%     fun = str2func(['Stick_' protocol.pulseseq]);
%     
%     protocol1 = protocol; %  make protocol with gradient gradient oriented along x direction
%     protocol1.grad_dirs = repmat([1 0 0],length(protocol.smalldel),1);
%     
%     % Parallel & perpendicular components of a cylinder
%     x_comp = [d];
% 
%     Lpar = log(fun([x_comp pi/2 0],protocol1)); % compute parallel component (cylinder along x, gradient along x)
%     
%     % Since R=0, no need to do any calculation
%     Lperp = zeros(size(protocol.G,1), 1);
% 
%     ePerp = exp(Lperp);
%     
%     % Compute the Legendre weighted signal
%     
%     Lpmp = Lperp - Lpar;
%     
%     if nargout > 1
%         [lgi, J_lgi] = LegendreGaussianIntegral(Lpmp, 6);
%     else
%         lgi = LegendreGaussianIntegral(Lpmp, 6);
%     end
%     
%     % Compute the spherical harmonic coefficients of the Watson's distribution
%     if nargout > 1
%         [coeff, J_coeff] = WatsonSHCoeff(kappa);
%     else
%         coeff = WatsonSHCoeff(kappa);
%     end
%     coeffMatrix = repmat(coeff, [l_q, 1]);
%     
%     % Compute the dot product between the symmetry axis of the Watson's distribution
%     % and the gradient direction
%     %
%     % For numerical reasons, cosTheta might not always be between -1 and 1
%     % Due to round off errors, individual gradient vectors in grad_dirs and the
%     % fibredir are never exactly normal.  When a gradient vector and fibredir are
%     % essentially parallel, their dot product can fall outside of -1 and 1.
%     %
%     % BUT we need make sure it does, otherwise the legendre function call below
%     % will FAIL and abort the calculation!!!
%     %
%     
%     fibredir = GetFibreOrientation('WatsonSticks', x);
%     
%     cosTheta = protocol.grad_dirs*fibredir;
%     badCosTheta = find(abs(cosTheta)>1);
%     cosTheta(badCosTheta) = cosTheta(badCosTheta)./abs(cosTheta(badCosTheta));
%     
%     % Compute the SH values at cosTheta
%     sh = zeros(size(coeff));
%     shMatrix = repmat(sh, [l_q, 1]);
%     for i = 1:7
%         shMatrix(:,i) = sqrt((i - .75)/pi);
%         % legendre function returns coefficients of all m from 0 to l
%         % we only need the coefficient corresponding to m = 0
%         % WARNING: make sure to input ROW vector as variables!!!
%         % cosTheta is expected to be a COLUMN vector.
%         tmp = legendre(2*i - 2, cosTheta');
%         tmp = tmp';
%         shMatrix(:,i) = shMatrix(:,i) .* tmp(:,1);
%     end
%     
%     E = sum(lgi.*coeffMatrix.*shMatrix, 2);
%     % with the SH approximation, there will be no guarantee that E will be positive
%     % but we need to make sure it does!!! replace the negative values with 10% of
%     % the smallest positive values
%     E(find(E<=0)) = min(E(find(E>0)))*0.1;
%     E = 0.5*E.*ePerp;
%     
%     
%     % Compute the Jacobian matrix; computed numerically
%     if(nargout>1)
%         dx = 0.00001;
%         J = zeros(length(E), length(x));
%         if nargin < 3
%             for i = 1:length(x); % compute the derivatives for all model parameters
%                 xpert = x;
%                 xpert(i) = xpert(i)*(1+dx);
%                 Epert = WatsonSticks(xpert, protocol);
%                 dEdx = (Epert - E)/(xpert(i)*dx);
%                 J(:,i) = dEdx;
%             end
%         else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
%             
%             for i = 1:length(x_deriv);
%                 if x_deriv(i) ~= 0
%                     xpert = x;
%                     xpert(i) = xpert(i)*(1+dx);
%                     Epert = WatsonSticks(xpert, protocol);
%                     dEdx = (Epert - E)/(xpert(i)*dx);
%                     J(:,i) = dEdx;
%                 end
%             end
%         end
%     end
% else
%     % sequences that use the MM method for the signal calculation
%     error Computation of singal not implemented for sequences that use the MM method for the signal calculation
% end

