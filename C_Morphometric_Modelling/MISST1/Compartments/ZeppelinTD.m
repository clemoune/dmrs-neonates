function [E J]=ZeppelinTD(x,protocol,x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating the Zeppelin compartment with time dependent 
% perpendicular diffusivity according to Novikov model.
% 
% [E,J]=ZeppelinTD(x, protocol,x_deriv)
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
%       x(2) - perpendicular diffusivity (hindered) at t = intfy
%       x(4) - A factor
%       x(5) - polar angle theta in spherical coordinates desbribing the
%       principal direction
%       x(6) - azimuthal angle phi in spherical coordinates describing the
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
    A=cell2mat(x(3)); 
    theta = cell2mat(x(4));
    phi = cell2mat(x(5));
else
    dPar = x(1);
    dPerp = x(2);
    A = x(3);
    theta = x(4);
    phi = x(5);
end

% calculate fibre direction from the specified angles
fibredir = [sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];



GAMMA = 2.675987E8;
if strcmp(protocol.pulseseq,'PGSE') || ...
        strcmp(protocol.pulseseq,'DSE') || strcmp(protocol.pulseseq,'STEAM')% sequences with one gradient orientation
        grad_dirs = protocol.grad_dirs;
        % Angles between gradient directions and fibre direction.
        cosTheta = grad_dirs*fibredir;
        cosThetaSq = cosTheta.^2;

        Bval = GetBvalues(protocol)';
        Bval_par = Bval.*cosThetaSq;
        Bval_perp = Bval.*(1-cosThetaSq);
        ePar=exp(-Bval_par.*dPar); 
        dPerp_corrected = dPerp*ones(size(Bval_perp));
        dPerp_corrected(protocol.smalldel >0) = dPerp + A./(protocol.delta(protocol.smalldel >0) -...
            protocol.smalldel(protocol.smalldel >0)/3).*(log(protocol.delta(protocol.smalldel >0)./protocol.smalldel(protocol.smalldel >0))+3/2);
        ePerp = exp(-Bval_perp.*dPerp_corrected);
        E = ePar.*ePerp;
else
    error('Time dependendt diffusivity not implemented yet for other sequences')
  
end
     

% Compute the Jacobian matrix
if(nargout>1)
   dx = 0.00001;
   J = zeros(length(E), 5);
    if nargin < 3 
         
         dEddPar = -Bval_par.*ePar.*ePerp;
         J(:,1) = dEddPar; 
      
         for i = 2:5
            xpert = x;
            xpert(i) = xpert(i)*(1+dx);
            Epert = ZeppelinTD(xpert, protocol);
            dEdx = (Epert - E)/(xpert(i)*dx);
            J(:,i) = dEdx;
         end         
        
    else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
       
        for i = 1:length(x_deriv);
            if x_deriv(i) ~= 0  
                
                xpert = x;
                xpert(i) = xpert(i)*(1+dx);
                Epert = ZeppelinTD(xpert, protocol);
                dEdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEdx;
            end
        end
        
       
        
    end   
   
end