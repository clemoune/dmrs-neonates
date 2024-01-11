function [E,J]=NormalFiniteAstroCylinders(x, protocol,x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating the NormalFiniteAstroCylinders compartment.
% 
% [E,J]=NormalFiniteAstroCylinders(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate consisting of isotropically 
% oriented finite cylinders with a normal distribution of radii and a diffusion 
% protocol specified in the input
% Substrate: Isotropically distributed, impermeable finite cylinders with a normal
% distribution of radii
% Diffusion pulse sequence: wrapper for various sequences
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% x - size 4 vector of model parameters in SI units for NormalFiniteAstroCylinders:
%       x(1) - free diffusivity of the material inside the cylinders.
%       x(2) - mean radius of the distribution
%       x(3) - variance of the normal distribution
%       x(4) - eccentricity (ratio between length and diameter)
% protocol - structure which includes all the information related to the 
%       diffusion protocol and required to generate the signal. 
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

if isfield(protocol,'approx')
    approx = protocol.approx;
else
    approx = 'GPD';
end
if strcmp(approx,'GPD')
    rad = x(2);
    var = x(3);
    ecc = x(4);


   if ~isfield(protocol,'fixedR')
         if ~isfield(protocol,'NormalRdist') || strcmp(protocol.NormalRdist,'pdf') % weights come from the pdf of gamma function
            intR_steps = 10;
            [R, weight] = NormalRadList(rad, var, intR_steps);
         elseif strcmp(protocol.NormalRdist,'vol') % diffusion signal from different radii is weighted by volume
            intR_steps = 10;
            [R, weight] = NormalRadByVolList(rad, var, intR_steps);
         else error Unknown GammaRdist option

        end
    else
         R = protocol.NormalR;
         weight = protocol.NormalWeights;   
    end
   
    if strcmp(protocol.pulseseq,'PGSE') || strcmp(protocol.pulseseq,'SWOGSE') || strcmp(protocol.pulseseq,'TWOGSE') || ...
            strcmp(protocol.pulseseq,'DSE') || strcmp(protocol.pulseseq,'STEAM')% sequences with 1 gradient direction


            fun = str2func('FiniteAstroCylinders'); % integrate over orientation and sum over radii

       

        E0 = zeros(length(protocol.smalldel),length(R));
        for i = 1:length(R)  
            params = [x(1) R(i) ecc];
            E0(:,i) = fun(params,protocol);
        end
         E = sum(E0.*repmat(weight,length(protocol.smalldel),1),2);

          % Compute the Jacobian matrix; computed numerically
        
    else % sequences with varying gradient orientation
        % need to loop over different directions
        % get Ndir directions on a spehere

         fun = str2func(['FiniteCylinder_' approx '_' protocol.pulseseq]);
        Ndir = 50;    
        fibre_dirs = ReadCaminoElecPS(sprintf('PointSets/Elec%03i.txt',Ndir))';  
        [az incl] = cart2sph(fibre_dirs(1,:),fibre_dirs(2,:),fibre_dirs(3,:)); % azimuth and inclinations
        theta = pi/2-incl;
        phi = az;

        
        protocol.Rweight = weight;
        xnew{1} = x(1);
        xnew{2} = R;
        xnew{3} = ecc;
        xnew{4} = theta;
        xnew{5} = phi;
        E = fun(xnew,protocol);
    end
    if(nargout>1)
            dx = 0.00001;
            J = zeros(length(E), length(x));
            if nargin < 3          
                for i = 1:length(x); % compute the derivatives for all model parameters
                    xpert = x;
                    xpert(i) = xpert(i)*(1+dx);
                    Epert = NormalFiniteAstroCylinders(xpert, protocol);
                    dEdx = (Epert - E)/(xpert(i)*dx);
                    J(:,i) = dEdx;
                end
            else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian

                for i = 1:length(x_deriv);
                    if x_deriv(i) ~= 0                  
                        xpert = x;
                        xpert(i) = xpert(i)*(1+dx);
                        Epert = NormalFiniteAstroCylinders(xpert, protocol);
                        dEdx = (Epert - E)/(xpert(i)*dx);
                        J(:,i) = dEdx;
                    end
                end
            end   
    end
   
else
    fun = str2func(['NormalFiniteCylinders_' approx '_' protocol.pulseseq]);
    Ndir = 20;    
    fibre_dirs = ReadCaminoElecPS(sprintf('PointSets/Elec%03i.txt',Ndir))';  
    [az incl] = cart2sph(fibre_dirs(1,:),fibre_dirs(2,:),fibre_dirs(3,:)); % azimuth and inclinations

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
            params = [x pi/2-incl(i) az(i)];
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
    else
        
        J0 = zeros(size(E0,1),length(x)+2,Ndir);        
        
        if nargin < 3 
            for i = 1:Ndir
                params = [x pi/2-incl(i) az(i)];
                [E0(:,i) J0(:,:,i)] = fun(params,protocol,[ones(size(x)) 0 0]);
            end
        else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
            for i = 1:Ndir
                params = [x pi/2-incl(i) az(i)];
                [E0(:,i) J0(:,:,i)] = fun(params,protocol,[x_deriv 0 0]);
            end         

        end 

       protocol.complex = tmp;
       if strcmp(protocol.complex,'complex')
            E = mean(E0,2);
            J = mean(J0(:,1:end-2,:),3);
        elseif strcmp(protocol.complex,'real')
            E = mean(E0(1:M,:),2);
            J = mean(J0(1:M,1:end-2,:),3);
        elseif strcmp(protocol.complex,'abs')
            E = abs(mean(E0(1:M,:)+1i*E0(M+1:end,:),2)); 
            J = zeros(length(E),length(x));
            for i = 1:length(x)
                ind = find(E>0);
                J(:,i) = (squeeze(mean(J0(ind,i,:),3)).*mean(E0(ind,:),2)+squeeze(mean(J0(ind+M,i,:),3)).*mean(E0(ind+M,:),2))./ E(ind);    
            end
        else error unknown protocol.complex
        end
       
       
       
   end
end