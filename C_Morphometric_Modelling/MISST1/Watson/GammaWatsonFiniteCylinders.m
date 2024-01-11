function [E,J]=GammaWatsonFiniteCylinders(x, protocol,x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a GammaWatsonFiniteCylinders compartment.
%
% [E,J]=GammaWatsonFiniteCylinders(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate of finite cylinders with a Gamma
% distribution of sizes and a Watson distribution of orientations and a 
% diffusion protocol specified in the input 
%
% Substrate: Finite cylinders with Gamma distributed radii and Watson orientation
% Diffusion pulse sequence: wrapper for various sequences
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% x - size 7 vector of model parameters in SI units for AstroCylinders:
%       x(1) - free diffusivity of the material inside the cylinders.
%       x(2) - mean radius of the cylinders.
%       x(3) - shape of the Gamma distribution
%       x(4) - eccentricity (ratio between cylinder length and diameter - 
% the same for all sizes)
%       x(5) - concentration parameter of Watson distribution
%       x(6) is the angle the polar angle theta in spherical coordinates desbribing the fibre
% direction
%       x(7) is the azimuthal angle phi in spherical coordinates describing the
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
% Authors:%   
%   Gary Hui Zhang (gary.zhang@ucl.ac.uk)
%   Maira Tariq (maira.tariq.11@ucl.ac.uk)
%   Andrada Ianus (a.ianus.11@ucl.ac.uk)
%   

if length(x) ~= 7
    error('the model should have exactly 7 parameters');
end


if isfield(protocol,'approx')
    approx = protocol.approx;
else
    approx = 'GPD';
end

if strcmp(approx,'GPD')
    meanR = x(2);
    a = x(3);
    ecc = x(4);
    b = meanR./a;

    if ~isfield(protocol,'fixedR')
         if ~isfield(protocol,'GammaRdist') || strcmp(protocol.GammaRdist,'pdf') % weights come from the pdf of gamma function
            intR_steps = 10;
            [R, weight] = GammaRadList(a, b, intR_steps);
         elseif strcmp(protocol.GammaRdist,'vol') % diffusion signal from different radii is weighted by volume
            intR_steps = 10;
%             [R, weight] = GammaRadByVolList(a, b, intR_steps);
            [R, weight0] = GammaRadList(a, b, intR_steps);
            weight = weight0.*(R*1E6).^3;
            weight = weight./sum(weight);
         else error Unknown GammaRdist option

        end
    else
         R = protocol.GammaR;
         weight = protocol.GammaWeights;   
    end
    
    
    if  strcmp(protocol.pulseseq,'PGSE') || strcmp(protocol.pulseseq,'SWOGSE') || strcmp(protocol.pulseseq,'TWOGSE') || ...
        strcmp(protocol.pulseseq,'DSE') || strcmp(protocol.pulseseq,'STEAM')% sequences with 1 gradient direction
    
        fun = str2func('WatsonFiniteCylinders');
        

        E0 = zeros(length(protocol.smalldel),length(R));
        for i = 1:length(R)  
            params = [x(1) R(i) ecc x(5:end)];
            E0(:,i) = fun(params,protocol);
        end
         E = sum(E0.*repmat(weight,length(protocol.smalldel),1),2);
    else % sequences with more gradient orientations
        fun = str2func(['FiniteCylinder_' approx '_' protocol.pulseseq]);
        if isfield(protocol,'Ndir')
            Ndir = protocol.Ndir;
        else
            Ndir = 100;
        end
        
        protocol.Rweight = weight;       
        kappa = x(5);

        fibre_dirs = ReadCaminoElecPS(sprintf('PointSets/Elec%03i.txt',Ndir))';  
        [az incl] = cart2sph(fibre_dirs(1,:),fibre_dirs(2,:),fibre_dirs(3,:)); % azimuth and inclinations

        maindir = GetFibreOrientation('GammaWatsonFiniteCylinders', x);
        WatsonWeights = hypergeom(1/2,3/2,kappa)^(-1).*exp(kappa.*(maindir'*fibre_dirs).^2);
        WatsonWeights = WatsonWeights./sum(WatsonWeights);
        theta = pi/2-incl;
        phi = az;       
        
        xnew{1} = x(1);
        xnew{2} = R;
        xnew{3} = ecc;
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
                Epert = GammaWatsonFiniteCylinders(xpert, protocol);
                dEdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEdx;
            end
        else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian

            for i = 1:length(x_deriv);
                if x_deriv(i) ~= 0                  
                    xpert = x;
                    xpert(i) = xpert(i)*(1+dx);
                    Epert = GammaWatsonFiniteCylinders(xpert, protocol);
                    dEdx = (Epert - E)/(xpert(i)*dx);
                    J(:,i) = dEdx;
                end
            end
        end   
    end

else % signal computed with MM
    fun = str2func(['GammaFiniteCylinders_' approx '_' protocol.pulseseq]);
    if isfield(protocol,'Ndir')
        Ndir = protocol.Ndir;
    else
        Ndir = 100;  
    end
    fibre_dirs = ReadCaminoElecPS(sprintf('PointSets/Elec%03i.txt',Ndir))';  
    [az incl] = cart2sph(fibre_dirs(1,:),fibre_dirs(2,:),fibre_dirs(3,:)); % azimuth and inclinations
    kappa = x(5);

    maindir = GetFibreOrientation('GammaWatsonFiniteCylinders', x);
    WatsonWeights = hypergeom(1/2,3/2,kappa)^(-1).*exp(kappa.*(maindir'*fibre_dirs).^2);
    WatsonWeights = WatsonWeights./sum(WatsonWeights);
    
    if isfield(protocol,'smalldel')
        M = size(protocol.smalldel,2);       
    elseif isfield(protocol,'smalldel1')
        M = size(protocol.smalldel1,2);   
    else
         M = size(protocol.G,1);       
    end  

    E0 = zeros(2*M,Ndir);
    tmp = protocol.complex;
    protocol.complex = 'complex';
    
    if nargout == 1
             
            for i = 1:Ndir
                params = [x(1:4) pi/2-incl(i) az(i)];
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
    else        
        
        if isfield(protocol,'pert') 
            dx = protocol.pert;
        else
            dx = 1E-5;
        end

        WatsonWeights_pertk = hypergeom(1/2,3/2,kappa*(1+dx))^(-1).*exp(kappa.*(1+dx).*(maindir'*fibre_dirs).^2);
        WatsonWeights_pertk = WatsonWeights_pertk./sum(WatsonWeights_pertk);
        WatsonWeights_pertk = repmat(WatsonWeights_pertk,M,1);
        WatsonWeights_dk = (WatsonWeights_pertk-WatsonWeights)./(kappa.*dx); 

        theta_idx = GetParameterIndex('GammaWatsonFiniteCylinders','theta');
        xpertth = x; xpertth(theta_idx) = x(theta_idx)+dx;    
        maindir_pertth = GetFibreOrientation('GammaWatsonFiniteCylinders', xpertth); 
        WatsonWeights_pertth = hypergeom(1/2,3/2,kappa)^(-1).*exp(kappa.*(maindir_pertth'*fibre_dirs).^2);
        WatsonWeights_pertth = WatsonWeights_pertth./sum(WatsonWeights_pertth);
        WatsonWeights_pertth = repmat(WatsonWeights_pertth,M,1);
        WatsonWeights_dth = (WatsonWeights_pertth-WatsonWeights)./dx;  

        phi_idx = GetParameterIndex('GammaWatsonFiniteCylinders','phi');
        xpertphi = x; xpertphi(phi_idx) = x(phi_idx)+dx;    
        maindir_pertphi = GetFibreOrientation('GammaWatsonFiniteCylinders', xpertphi); 
        WatsonWeights_pertphi = hypergeom(1/2,3/2,kappa)^(-1).*exp(kappa.*(maindir_pertphi'*fibre_dirs).^2);
        WatsonWeights_pertphi = WatsonWeights_pertphi./sum(WatsonWeights_pertphi);
        WatsonWeights_pertphi = repmat(WatsonWeights_pertphi,M,1);
        WatsonWeights_dphi = (WatsonWeights_pertphi-WatsonWeights)./dx; 
         
        J0 = zeros(size(E0,1),6,Ndir);
          
        if nargin < 3 
          
          
            for i = 1:Ndir
                params = [x(1:4) pi/2-incl(i) az(i)];
                [E0(:,i), J0(:,:,i)] = fun(params,protocol,[1 1 1 1 0 0]);
            end
            protocol.complex = tmp;
            if strcmp(protocol.complex,'complex')
                E = sum(E0.*repmat(WatsonWeights,2,1),2);  
                J = zeros(length(E),7);
                J(:,1:4) = sum(J0(:,1:4,:).*permute(repmat(WatsonWeights,[2 1 4]),[1 3 2]),3); 
                J(:,5) = sum(E0.*repmat(WatsonWeights_dk,2,1),2); 
                J(:,6) = sum(E0.*repmat(WatsonWeights_dth,2,1),2); 
                J(:,7) = sum(E0.*repmat(WatsonWeights_dphi,2,1),2);
            elseif strcmp(protocol.complex,'real')
                E = sum(E0(1:M,:).*WatsonWeights,2);
                J = zeros(length(E),7);
                J(:,1:4) = sum(J0(1:M,1:4,:).*permute(repmat(WatsonWeights,[1 1 4]),[1 3 2]),3); 
                J(:,5) = sum(E0(1:M,:).*WatsonWeights_dk,2); 
                J(:,6) = sum(E0(1:M,:).*WatsonWeights_dth,2); 
                J(:,7) = sum(E0(1:M,:).*WatsonWeights_dphi,2);
            elseif strcmp(protocol.complex,'abs')
                E = abs(sum(E0(1:M,:).*WatsonWeights+1i*E0(M+1:end,:).*WatsonWeights,2));
                E1 = E0(1:M,:); E2 = E0(M+1:end,:);
                E1_mat = permute(repmat(E1,[1 1 4]),[1 3 2]); E2_mat = permute(repmat(E2,[1 1 4]),[1 3 2]);
                WW_mat = permute(repmat(WatsonWeights,[1 1 4]),[1 3 2]);
                J = zeros(length(E),7);
                ind = find(E>0);
                J(:,1:4) = (sum((E1_mat(ind,:,:).*J0(ind,1:4,:) + E2_mat(ind,:,:).*J0(M+ind,1:4,:)).*WW_mat(ind,:,:),3))./repmat(E(ind),[1 4]); 

                J(:,5) = (sum(E1(ind,:).*WatsonWeights(ind,:),2).*sum(E1(ind,:).*WatsonWeights_dk(ind,:),2) + ...
                   sum(E2(ind,:).*WatsonWeights(ind,:),2).*sum(E2(ind,:).*WatsonWeights_dk(ind,:),2))./E(ind);
                J(:,6) = (sum(E1(ind,:).*WatsonWeights(ind,:),2).*sum(E1(ind,:).*WatsonWeights_dth(ind,:),2) + ...
                   sum(E2(ind,:).*WatsonWeights(ind,:),2).*sum(E2(ind,:).*WatsonWeights_dth(ind,:),2))./E(ind);
                J(:,7) = (sum(E1(ind,:).*WatsonWeights(ind,:),2).*sum(E1(ind,:).*WatsonWeights_dphi(ind,:),2) + ...
                   sum(E2(ind,:).*WatsonWeights(ind,:),2).*sum(E2(ind,:).*WatsonWeights_dphi(ind,:),2))./E(ind);

            else error unknown protocol.complex
            end 
        else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian

         
           
            for i = 1:Ndir
                params = [x(1:4) pi/2-incl(i) az(i)];
                [E0(:,i), J0(:,:,i)] = fun(params,protocol,[x_deriv(1:4) 0 0]);
            end
            protocol.complex = tmp;
            
            if strcmp(protocol.complex,'complex')
                E = sum(E0.*repmat(WatsonWeights,2,1),2);  
                J = zeros(length(E),7);
                J(:,1:4) = sum(J0(:,1:4,:).*permute(repmat(WatsonWeights,[2 1 4]),[1 3 2]),3); 
                if x_deriv(5) ~= 0
                J(:,5) = sum(E0.*repmat(WatsonWeights_dk,2,1),2); 
                end
                if x_deriv(6) ~= 0
                J(:,6) = sum(E0.*repmat(WatsonWeights_dth,2,1),2); 
                end
                 if x_deriv(7) ~= 0
                J(:,7) = sum(E0.*repmat(WatsonWeights_dphi,2,1),2);
                 end
            elseif strcmp(protocol.complex,'real')
                E = sum(E0(1:M,:).*WatsonWeights,2);
                J = zeros(length(E),7);
                J(:,1:4) = sum(J0(1:M,1:4,:).*permute(repmat(WatsonWeights,[1 1 4]),[1 3 2]),3); 
                if x_deriv(5) ~= 0
                J(:,5) = sum(E0(1:M,:).*WatsonWeights_dk,2); 
                end
                if x_deriv(6) ~= 0
                J(:,6) = sum(E0(1:M,:).*WatsonWeights_dth,2); 
                end
                 if x_deriv(7) ~= 0
                J(:,7) = sum(E0(1:M,:).*WatsonWeights_dphi,2);
                 end
            elseif strcmp(protocol.complex,'abs')
                E = abs(sum(E0(1:M,:).*WatsonWeights+1i*E0(M+1:end,:).*WatsonWeights,2));
                E1 = E0(1:M,:); E2 = E0(M+1:end,:);
                E1_mat = permute(repmat(E1,[1 1 4]),[1 3 2]); E2_mat = permute(repmat(E2,[1 1 4]),[1 3 2]);
                WW_mat = permute(repmat(WatsonWeights,[1 1 4]),[1 3 2]);
                ind = find(E>0);
                J = zeros(length(E),7);
                J(:,1:4) = (sum((E1_mat(ind,:,:).*J0(ind,1:4,:) + E2_mat(ind,:,:).*J0(M+ind,1:4,:)).*WW_mat(ind,:,:),3))./repmat(E(ind),[1 4]); 
                if x_deriv(5) ~= 0
                J(:,5) =  (sum(E1(ind,:).*WatsonWeights(ind,:),2).*sum(E1(ind,:).*WatsonWeights_dk(ind,:),2) + ...
                   sum(E2(ind,:).*WatsonWeights(ind,:),2).*sum(E2(ind,:).*WatsonWeights_dk(ind,:),2))./E(ind);
                end
                if x_deriv(6) ~= 0
                J(:,6) = (sum(E1(ind,:).*WatsonWeights(ind,:),2).*sum(E1(ind,:).*WatsonWeights_dth(ind,:),2) + ...
                   sum(E2(ind,:).*WatsonWeights(ind,:),2).*sum(E2(ind,:).*WatsonWeights_dth(ind,:),2))./E(ind);
                end
                if x_deriv(7) ~= 0
                J(:,7) = (sum(E1(ind,:).*WatsonWeights(ind,:),2).*sum(E1(ind,:).*WatsonWeights_dphi(ind,:),2) + ...
                   sum(E2(ind,:).*WatsonWeights(ind,:),2).*sum(E2(ind,:).*WatsonWeights_dphi(ind,:),2))./E(ind);

                end
            else error unknown protocol.complex
            end         

        end       
       
    end

    
    
end

