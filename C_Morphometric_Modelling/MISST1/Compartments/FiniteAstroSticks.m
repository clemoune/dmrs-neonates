function [E,J]=FiniteAstroSticks(x, protocol,x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating the FiniteAstroCylinders compartment.
% 
% [E,J]=FiniteAstroSticks(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate consisting of isotropically 
% oriented finite sticks and a diffusion protocol specified in the input
% Substrate: Isotropically distributed impermeable finite sticks
% Diffusion pulse sequence: wrapper for various sequences
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% x - size 2 vector of model parameters in SI units for FiniteAstroCylinders:
%       x(1) - free diffusivity of the material inside the cylinders.
%       L(2) - length of the sticks.
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
%   Andrada Ianus (a.ianus.11@ucl.ac.uk)
%   Daniel C. Alexander (d.alexander@ucl.ac.uk)
%       

if isfield(protocol,'approx')
    approx = protocol.approx;
else
    approx = 'GPD';
end

if strcmp(approx,'GPD')
    if strcmp(protocol.pulseseq,'PGSE') || strcmp(protocol.pulseseq,'SWOGSE') || strcmp(protocol.pulseseq,'TWOGSE') || ...
            strcmp(protocol.pulseseq,'DSE') || strcmp(protocol.pulseseq,'STEAM')% sequences with 1 gradient direction
        fun = str2func(['FiniteStick_' approx '_' protocol.pulseseq]);
        G = protocol.G';
        protocol1 = protocol; %  make protocol with gradient gradient oriented along x direction
        protocol1.grad_dirs = repmat([1 0 0],length(protocol.smalldel),1);

        Lperp = log(fun([x 0 0],protocol1)); % compute perpendicular component (cylinder along z, gradient along x)
        Lpar = log(fun([x pi/2 0],protocol1)); % compute parallel component (cylinder along x, gradient along x)

        Ldiff = zeros(size(G));     
        Ldiff(G>0) = Lperp(G>0)-Lpar(G>0);
        Ldiff(Ldiff<0) = 0;
        E_r = zeros(size(G)); 
        E_r(G == 0) = 1;
        E_r(Ldiff>0) = sqrt(pi)./(2.*sqrt(Ldiff(Ldiff>0))).*exp(Lperp(Ldiff>0)).*erf(sqrt(Ldiff(Ldiff>0)));    

        E=E_r;
    else
        fun = str2func(['FiniteStick_' approx '_' protocol.pulseseq]);      
        Ndir = 100;
       
        fibre_dirs = ReadCaminoElecPS(sprintf('PointSets/Elec%03i.txt',Ndir))';  
        [az incl] = cart2sph(fibre_dirs(1,:),fibre_dirs(2,:),fibre_dirs(3,:)); % azimuth and inclinations
        
         theta = pi/2-incl;
        phi = az;
        xnew{1} = x(1);
        xnew{2} = x(2);
        xnew{3} = theta;
        xnew{4} = phi;
        E = fun(xnew,protocol);
    end

    % Compute the Jacobian matrix; computed numerically
    if(nargout>1)
        dx = 0.00001;
        J = zeros(length(E), 3);
        if nargin < 3          
            for i = 1:length(x); % compute the derivatives for all model parameters
                xpert = x;
                xpert(i) = xpert(i)*(1+dx);
                Epert = FiniteAstroSticks(xpert, protocol);
                dEdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEdx;
            end
        else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian

            for i = 1:length(x_deriv);
                if x_deriv(i) ~= 0                  
                    xpert = x;
                    xpert(i) = xpert(i)*(1+dx);
                    Epert = FiniteAstroSticks(xpert, protocol);
                    dEdx = (Epert - E)/(xpert(i)*dx);
                    J(:,i) = dEdx;
                end
            end
        end   
    end
else % sequences that use the MM method for the signal calculation
    % need to loop over different directions
    % get Ndir directions on a spehere
    fun = str2func(['FiniteStick_' approx '_' protocol.pulseseq]);   
    Ndir = 64; 

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


    
    % Compute the Jacobian matrix
    if(nargout>1)
        J = zeros(length(E),length(x));
        dx = protocol.pert;
        if nargin < 3          
            for i = 1:length(x)
            xpert = x;
            protocol.diff=i;   
            xpert(i) = xpert(i)*(1+dx);    
            Epert = FiniteAstroSticks(xpert,protocol);
            dEtdx = (Epert - E)/(xpert(i)*dx);
            J(:,i) = dEtdx;
            end
        else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian      
             for i = 1:length(x)        
                if x_deriv(i) ~= 0  

                    xpert = x;
                    protocol.diff=i;   
                    xpert(i) = xpert(i)*(1+dx);    
                    Epert = FiniteAstroSticks(xpert,protocol,x_deriv);
                    dEtdx = (Epert - E)/(xpert(i)*dx);
                    J(:,i) = dEtdx;
                end     

             end             
        end   
        protocol.diff=0;   
    end   
    
end
end
