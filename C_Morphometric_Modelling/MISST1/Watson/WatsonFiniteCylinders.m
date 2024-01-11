function [E,J]=WatsonFiniteCylinders(x, protocol,x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a WatsonFiniteCylinders compartment.
%
% [E,J]=WatsonFiniteCylinders(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate of finite cylinders with a Watson 
% distribution of orientations and a diffusion protocol specified in the input 
%
% Substrate: Finite cylinders with a Watson distribution of orientations
% Diffusion pulse sequence: wrapper for various sequences
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% x - size 6 vector of model parameters in SI units for AstroCylinders:
%       x(1) - free diffusivity of the material inside the cylinders.
%       x(2) - radius of the cylinders.
%       x(3) - eccentricity (ratio between cylinder length and diameter - 
% the same for all sizes)
%       x(4) - concentration parameter of Watson distribution
%       x(5) is the angle the polar angle theta in spherical coordinates desbribing the fibre
% direction
%       x(6) is the azimuthal angle phi in spherical coordinates describing the
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

if length(x) ~= 6
    error('the model should have exactly 6 parameters');
end

d=x(1);
R=x(2);
ecc = x(3);
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

if strcmp(approx,'GPD')
    if strcmp(protocol.pulseseq,'PGSE') || strcmp(protocol.pulseseq,'SWOGSE') || strcmp(protocol.pulseseq,'TWOGSE') || ...
        strcmp(protocol.pulseseq,'DSE') || strcmp(protocol.pulseseq,'STEAM')% sequences with 1 gradient direction
    
        l_q = size(protocol.grad_dirs,1);
        fun = str2func(['FiniteCylinder_' approx '_' protocol.pulseseq]);

        protocol1 = protocol; %  make protocol with gradient gradient oriented along x direction
        protocol1.grad_dirs = repmat([1 0 0],length(protocol.smalldel),1);

        % Parallel & perpendicular components of a cylinder
        x_comp = [d R ecc];

        Lperp = log(fun([x_comp 0 0],protocol1)); % compute perpendicular component (cylinder along z, gradient along x)
        Lpar = log(fun([x_comp pi/2 0],protocol1)); % compute parallel component (cylinder along x, gradient along x)


        ePerp = exp(Lperp);

        % Compute the Legendre weighted signal
        G = protocol.G;
        Lpmp = zeros(size(G));     
        Lpmp(G>0) = Lperp(G>0)-Lpar(G>0);
        Lpmp(Lpmp<0) = 0;

        if nargout > 1
            [lgi, J_lgi] = LegendreGaussianIntegral(Lpmp', 6);
        else
            lgi = LegendreGaussianIntegral(Lpmp', 6);
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
        E(E<=0) = min(E(E>0))*0.1;
        E = 0.5*E.*ePerp;
    else
          % sequences that use the MM method for the signal calculation
        fun = str2func(['FiniteCylinder_' approx '_' protocol.pulseseq]);
        if isfield(protocol,'Ndir')
            Ndir = protocol.Ndir;
        else
            Ndir = 100;
        end
        kappa = x(4);

        fibre_dirs = ReadCaminoElecPS(sprintf('PointSets/Elec%03i.txt',Ndir))';  
        [az incl] = cart2sph(fibre_dirs(1,:),fibre_dirs(2,:),fibre_dirs(3,:)); % azimuth and inclinations

        maindir = GetFibreOrientation('WatsonFiniteCylinders', x); 

        WatsonWeights = hypergeom(1/2,3/2,kappa)^(-1).*exp(kappa.*(maindir'*fibre_dirs).^2);
    
        
         WatsonWeights = WatsonWeights./sum(WatsonWeights);
         theta = pi/2-incl;
        phi = az;
        xnew{1} = x(1);
        xnew{2} = x(2);
        xnew{3} = x(3);
        xnew{4} = theta;
        xnew{5} = phi;
        protocol.Nweight = WatsonWeights;
        E = fun(xnew,protocol);

        
    end
    
    
    % Compute the Jacobian matrix; computed numerically
    if(nargout>1)
        dx = 0.00001;
        J = zeros(length(E), length(x));
        if nargin < 3
            for i = 1:length(x); % compute the derivatives for all model parameters
                xpert = x;
                xpert(i) = xpert(i)*(1+dx);
                if i == 2% radius -> change eccentricity to yield the same L
                    xpert(i+1) = xpert(i+1)./(1+dx);
                end
                Epert = WatsonFiniteCylinders(xpert, protocol);
                dEdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEdx;
            end
        else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
            
            for i = 1:length(x_deriv);
                if x_deriv(i) ~= 0
                    xpert = x;
                    xpert(i) = xpert(i)*(1+dx);
                    if i == 2% radius -> change eccentricity to yield the same L
                        xpert(i+1) = xpert(i+1)./(1+dx);
                    end                    
                    Epert = WatsonFiniteCylinders(xpert, protocol);
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
    kappa = x(4);
    
    fibre_dirs = ReadCaminoElecPS(sprintf('PointSets/Elec%03i.txt',Ndir))';  
    [az incl] = cart2sph(fibre_dirs(1,:),fibre_dirs(2,:),fibre_dirs(3,:)); % azimuth and inclinations

    maindir = GetFibreOrientation('WatsonFiniteCylinders', x); 

    WatsonWeights = hypergeom(1/2,3/2,kappa)^(-1).*exp(kappa.*(maindir'*fibre_dirs).^2);    
    WatsonWeights = WatsonWeights./sum(WatsonWeights);

    if isfield(protocol,'smalldel')
        M = size(protocol.smalldel,2);       
    elseif isfield(protocol,'smalldel1')
        M = size(protocol.smalldel1,2);   
    else
         M = size(protocol.G,1);       
    end

    
    if nargout ==1    
       
            E0 = zeros(M*2,Ndir);
            tmp = protocol.complex;
            protocol.complex = 'complex';
            for i = 1:Ndir
                params = [d R ecc pi/2-incl(i) az(i)];
                E0(:,i) = fun(params,protocol);
            end
            protocol.complex = tmp;
            if strcmp(protocol.complex,'complex')
                E = sum(E0.*repmat(WatsonWeights,2,1),2);
            elseif strcmp(protocol.complex,'real')
                E = sum(E0(1:M,:).*WatsonWeights,2);
            elseif strcmp(protocol.complex,'abs')
                E = abs(sum(E0(1:M,:).*WatsonWeights+1i*E0(M+1:end,:).*WatsonWeights,2));
            else error unknown protocol.complex
            end
     
    
    else % calculate jacobians      
     
        
        if isfield(protocol,'pert') 
            dx = protocol.pert;
        else
            dx = 1E-5;
        end

        WatsonWeights_pertk = hypergeom(1/2,3/2,kappa*(1+dx))^(-1).*exp(kappa.*(1+dx).*(maindir'*fibre_dirs).^2);
        WatsonWeights_pertk = WatsonWeights_pertk./sum(WatsonWeights_pertk);
        WatsonWeights_pertk = repmat(WatsonWeights_pertk,M,1);
        WatsonWeights_dk = (WatsonWeights_pertk-WatsonWeights)./(kappa.*dx); 

        theta_idx = GetParameterIndex('WatsonFiniteCylinders','theta');
        xpertth = x; xpertth(theta_idx) = x(theta_idx)+dx;    
        maindir_pertth = GetFibreOrientation('WatsonFiniteCylinders', xpertth); 
        WatsonWeights_pertth = hypergeom(1/2,3/2,kappa)^(-1).*exp(kappa.*(maindir_pertth'*fibre_dirs).^2);
        WatsonWeights_pertth = WatsonWeights_pertth./sum(WatsonWeights_pertth);
        WatsonWeights_pertth = repmat(WatsonWeights_pertth,M,1);
        WatsonWeights_dth = (WatsonWeights_pertth-WatsonWeights)./dx;  

        phi_idx = GetParameterIndex('WatsonFiniteCylinders','phi');
        xpertphi = x; xpertphi(phi_idx) = x(phi_idx)+dx;    
        maindir_pertphi = GetFibreOrientation('WatsonFiniteCylinders', xpertphi); 
        WatsonWeights_pertphi = hypergeom(1/2,3/2,kappa)^(-1).*exp(kappa.*(maindir_pertphi'*fibre_dirs).^2);
        WatsonWeights_pertphi = WatsonWeights_pertphi./sum(WatsonWeights_pertphi);
        WatsonWeights_pertphi = repmat(WatsonWeights_pertphi,M,1);
        WatsonWeights_dphi = (WatsonWeights_pertphi-WatsonWeights)./dx; 
    
        if nargin < 3            
      
                E0 = zeros(M*2,Ndir);
                J0 = zeros(size(E0,1),5,Ndir);
                tmp = protocol.complex;
                protocol.complex = 'complex';
                for i = 1:Ndir
                    params = [d R ecc pi/2-incl(i) az(i)];
                    [E0(:,i), J0(:,:,i)] = fun(params,protocol,[1 1 1 0 0]);
                end
                protocol.complex = tmp;
                if strcmp(protocol.complex,'complex')
                    E = sum(E0.*repmat(WatsonWeights,2,1),2);  
                    J = zeros(length(E),6);
                    J(:,1:3) = sum(J0(:,1:3,:).*permute(repmat(WatsonWeights,[2 1 3]),[1 3 2]),3); 
                    J(:,4) = sum(E0.*repmat(WatsonWeights_dk,2,1),2); 
                    J(:,5) = sum(E0.*repmat(WatsonWeights_dth,2,1),2); 
                    J(:,6) = sum(E0.*repmat(WatsonWeights_dphi,2,1),2);
                elseif strcmp(protocol.complex,'real')
                    E = sum(E0(1:M,:).*WatsonWeights,2);
                    J = zeros(length(E),6);
                    J(:,1:3) = sum(J0(1:M,1:3,:).*permute(repmat(WatsonWeights,[1 1 3]),[1 3 2]),3); 
                    J(:,4) = sum(E0(1:M,:).*WatsonWeights_dk,2); 
                    J(:,5) = sum(E0(1:M,:).*WatsonWeights_dth,2); 
                    J(:,6) = sum(E0(1:M,:).*WatsonWeights_dphi,2);
                elseif strcmp(protocol.complex,'abs')
                    E = abs(sum(E0(1:M,:).*WatsonWeights+1i*E0(M+1:end,:).*WatsonWeights,2));
                    E1 = E0(1:M,:); E2 = E0(M+1:end,:);
                    E1_mat = permute(repmat(E1,[1 1 3]),[1 3 2]); E2_mat = permute(repmat(E2,[1 1 3]),[1 3 2]);
                    WW_mat = permute(repmat(WatsonWeights,[1 1 3]),[1 3 2]);
                    J = zeros(length(E),6);
                    ind = find(E>0);
                    J(:,1:3) = (sum((E1_mat(ind,:,:).*J0(ind,1:3,:) + E2_mat(ind,:,:).*J0(M+ind,1:3,:)).*WW_mat(ind,:,:),3))./repmat(E(ind),[1 3]); 
                    
                    J(:,4) = (sum(E1(ind,:).*WatsonWeights(ind,:),2).*sum(E1(ind,:).*WatsonWeights_dk(ind,:),2) + ...
                       sum(E2(ind,:).*WatsonWeights(ind,:),2).*sum(E2(ind,:).*WatsonWeights_dk(ind,:),2))./E(ind);
                    J(:,5) = (sum(E1(ind,:).*WatsonWeights(ind,:),2).*sum(E1(ind,:).*WatsonWeights_dth(ind,:),2) + ...
                       sum(E2(ind,:).*WatsonWeights(ind,:),2).*sum(E2(ind,:).*WatsonWeights_dth(ind,:),2))./E(ind);
                    J(:,6) = (sum(E1(ind,:).*WatsonWeights(ind,:),2).*sum(E1(ind,:).*WatsonWeights_dphi(ind,:),2) + ...
                       sum(E2(ind,:).*WatsonWeights(ind,:),2).*sum(E2(ind,:).*WatsonWeights_dphi(ind,:),2))./E(ind);
                    
                else error unknown protocol.complex
                end  
        
        
        else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian      
        
                E0 = zeros(M*2,Ndir);
                J0 = zeros(size(E0,1),5,Ndir);
                tmp = protocol.complex;
                protocol.complex = 'complex';
                for i = 1:Ndir
                    params = [d R ecc pi/2-incl(i) az(i)];
                    [E0(:,i), J0(:,:,i)] = fun(params,protocol,[x_deriv(1:3) 0 0]);
                end
                protocol.complex = tmp;
                if strcmp(protocol.complex,'complex')
                    E = sum(E0.*repmat(WatsonWeights,2,1),2);  
                    J = zeros(length(E),6);
                    J(:,1:3) = sum(J0(:,1:3,:).*permute(repmat(WatsonWeights,[2 1 3]),[1 3 2]),3); 
                    if x_deriv(4) ~= 0
                    J(:,4) = sum(E0.*repmat(WatsonWeights_dk,2,1),2); 
                    end
                    if x_deriv(5) ~= 0
                    J(:,5) = sum(E0.*repmat(WatsonWeights_dth,2,1),2); 
                    end
                     if x_deriv(6) ~= 0
                    J(:,6) = sum(E0.*repmat(WatsonWeights_dphi,2,1),2);
                     end
                elseif strcmp(protocol.complex,'real')
                    E = sum(E0(1:M,:).*WatsonWeights,2);
                    J = zeros(length(E),6);
                    J(:,1:3) = sum(J0(1:M,1:3,:).*permute(repmat(WatsonWeights,[1 1 3]),[1 3 2]),3); 
                    if x_deriv(4) ~= 0
                    J(:,4) = sum(E0(1:M,:).*WatsonWeights_dk,2); 
                    end
                    if x_deriv(5) ~= 0
                    J(:,5) = sum(E0(1:M,:).*WatsonWeights_dth,2); 
                    end
                     if x_deriv(6) ~= 0
                    J(:,6) = sum(E0(1:M,:).*WatsonWeights_dphi,2);
                     end
                elseif strcmp(protocol.complex,'abs')
                    E = abs(sum(E0(1:M,:).*WatsonWeights+1i*E0(M+1:end,:).*WatsonWeights,2));
                    E1 = E0(1:M,:); E2 = E0(M+1:end,:);
                    E1_mat = permute(repmat(E1,[1 1 3]),[1 3 2]); E2_mat = permute(repmat(E2,[1 1 3]),[1 3 2]);
                    WW_mat = permute(repmat(WatsonWeights,[1 1 3]),[1 3 2]);
                    ind = find(E>0);
                    J = zeros(length(E),6);
                    J(:,1:3) = (sum((E1_mat(ind,:,:).*J0(ind,1:3,:) + E2_mat(ind,:,:).*J0(M+ind,1:3,:)).*WW_mat(ind,:,:),3))./repmat(E(ind),[1 3]); 
                    if x_deriv(4) ~= 0
                    J(:,4) =  (sum(E1(ind,:).*WatsonWeights(ind,:),2).*sum(E1(ind,:).*WatsonWeights_dk(ind,:),2) + ...
                       sum(E2(ind,:).*WatsonWeights(ind,:),2).*sum(E2(ind,:).*WatsonWeights_dk(ind,:),2))./E(ind);
                    end
                    if x_deriv(5) ~= 0
                    J(:,5) = (sum(E1(ind,:).*WatsonWeights(ind,:),2).*sum(E1(ind,:).*WatsonWeights_dth(ind,:),2) + ...
                       sum(E2(ind,:).*WatsonWeights(ind,:),2).*sum(E2(ind,:).*WatsonWeights_dth(ind,:),2))./E(ind);
                    end
                    if x_deriv(6) ~= 0
                    J(:,6) = (sum(E1(ind,:).*WatsonWeights(ind,:),2).*sum(E1(ind,:).*WatsonWeights_dphi(ind,:),2) + ...
                       sum(E2(ind,:).*WatsonWeights(ind,:),2).*sum(E2(ind,:).*WatsonWeights_dphi(ind,:),2))./E(ind);
                    
                    end
                else error unknown protocol.complex
                end
            
        end           
      
     end    
    
end

