function [E,J]=GammaFiniteCylinders(x, protocol,x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating the GammaFiniteCylinders compartment.
% 
% [E,J]=GammaFiniteCylinders(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate consisting of parallel
% finite cylinders with a Gamma distribution of radii and a diffusion 
% protocol specified in the input
% Substrate: Parallel, impermeable finite cylinders with a Gamma
% distribution of radii
% Diffusion pulse sequence: wrapper for various sequences
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% x - size 6 vector of model parameters in SI units for GammaFiniteCylinders:
%       x(1) - free diffusivity of the material inside the cylinders.
%       x(2) - mean radius of the distribution.
%       x(3) - shape parameter of the Gamma distribution
%       x(4) - eccentricity (ratio between length and diameter, the same for all sizes)
%       x(5) - polar angle theta in spherical coordinates desbribing the fibre
% direction
%       x(6) - azimuthal angle phi in spherical coordinates describing the
% fibre direction
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
  
    % 
    Dres = x(1);
    meanR = x(2);
    a = x(3);
    ecc = x(4);
    theta = x(5);
    phi = x(6);

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

     fun_name = ['FiniteCylinder_' approx '_' protocol.pulseseq];
     fun = str2func(fun_name);
     protocol.Rweight = weight;
     xnew{1} = x(1);
     xnew{2} = R;
     xnew{3} = ecc;
     xnew{4} = theta;
     xnew{5} = phi;

     E = fun(xnew,protocol);

     if(nargout>1)
        dx = 0.00001;
        J = zeros(length(E), length(x));
        if nargin < 3 

            for i = 1:length(x); % compute the derivatives for all model parameters
                xpert = x;
                xpert(i) = xpert(i)*(1+dx);
                Epert = GammaFiniteCylinders(xpert, protocol);
                dEdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEdx;
            end


        else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian

            for i = 1:length(x_deriv);
                if x_deriv(i) ~= 0  

                    xpert = x;
                    xpert(i) = xpert(i)*(1+dx);
                    Epert = GammaFiniteCylinders(xpert, protocol);
                    dEdx = (Epert - E)/(xpert(i)*dx);
                    J(:,i) = dEdx;
                end
            end        

        end   
     end
else
        fun_name = ['GammaFiniteCylinders_' approx '_' protocol.pulseseq]; 
       fun = str2func(fun_name);
    if nargout == 1
           
       E = fun(x,protocol);
    else        
        if nargin < 3 
        [E J] = fun(x,protocol);
        else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
        [E J] = fun(x,protocol,x_deriv);           

        end   
   end
   
   
end    


% fun_name = ['GammaFiniteCylinders_' approx '_' protocol.pulseseq];
% fun = str2func(fun_name);
% 
% if nargin <3
%     if nargout == 1
%         E = fun(x,protocol);
%     else
%         [E J] = fun(x,protocol);
% 
%     end
% else
%     if nargout == 1
%         E = fun(x,protocol, x_deriv);
%     else
%         [E J] = fun(x,protocol, x_deriv);
%     end
% end