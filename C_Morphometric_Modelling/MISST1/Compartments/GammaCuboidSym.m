function [E,J]=GammaCuboidSym(x, protocol,x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating the GammaCuboidSym compartment.
% 
% [E,J]=GammaCuboid(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate consisting of a cuboid/ 
% an ensemble of cuboids with 2 equal sides and a Gamma distribution of sizes
% and a diffusion protocol specified in the input
% Substrate: cuboid / ensemble of cuboids with 2 equal sides and a Gamma
% distribution of sizes
% Diffusion pulse sequence: wrapper for various sequences
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% x - size 4 vector of model parameters in SI units for GammaCuboidSym:
%       x(1) - free diffusivity of the material inside the cylinders.
%       x(2) - length of the cuboid (lx = ly)
%       x(3) - shape parameter of the Gamma distribution
%       x(4) - eccentricity (ratio between height of the cuboid lz and 
%       length lx)
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
    Dres = x(1);
    meanL = x(2);
    a = x(3);
    ecc = x(4);

    b = meanL./a;


    if ~isfield(protocol,'fixedR')
         if ~isfield(protocol,'GammaRdist') || strcmp(protocol.GammaRdist,'pdf') % weights come from the pdf of gamma function
            intR_steps = 10;
            [L, weight] = GammaRadList(a, b, intR_steps);
         elseif strcmp(protocol.GammaRdist,'vol') % diffusion signal from different radii is weighted by volume
            intR_steps = 10;
%             [L, weight] = GammaRadByVolList(a, b, intR_steps);
            [L, weight0] = GammaRadList(a, b, intR_steps);
            weight = weight0.*(L*1E6).^3;
            weight = weight./sum(weight);
         else error Unknown GammaRdist option

        end
    else
         L = protocol.GammaR;
         weight = protocol.GammaWeights;   
    end

     fun_name = ['CuboidSym_' approx '_' protocol.pulseseq];
     fun = str2func(fun_name);
     if isfield(protocol,'G')
        E0 = zeros(length(protocol.G),intR_steps);
     else
         E0 = zeros(length(protocol.G1),intR_steps);
     end
     
     for i = 1:intR_steps
         params = [x(1) L(i) ecc];
         E0(:,i) = fun(params,protocol);
     end
         weight_mat = repmat(weight,size(E0,1),1);
         E = sum(E0.*weight_mat,2);
    

     if(nargout>1)
        dx = 0.00001;
        J = zeros(length(E), length(x));
        if nargin < 3 

            for i = 1:length(x); % compute the derivatives for all model parameters
                xpert = x;
                xpert(i) = xpert(i)*(1+dx);
                Epert = GammaCuboidSym(xpert, protocol);
                dEdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEdx;
            end


        else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian

            for i = 1:length(x_deriv);
                if x_deriv(i) ~= 0  

                    xpert = x;
                    xpert(i) = xpert(i)*(1+dx);
                    Epert = GammaCuboidSym(xpert, protocol);
                    dEdx = (Epert - E)/(xpert(i)*dx);
                    J(:,i) = dEdx;
                end
            end        

        end   
     end
else
    
          meanL = x(2);
        a = x(3);
        ecc = x(4);

        b = meanL./a;


    if ~isfield(protocol,'fixedR')
         if ~isfield(protocol,'GammaRdist') || strcmp(protocol.GammaRdist,'pdf') % weights come from the pdf of gamma function
            intR_steps = 10;
            [L, weight] = GammaRadList(a, b, intR_steps);
         elseif strcmp(protocol.GammaRdist,'vol') % diffusion signal from different radii is weighted by volume
            intR_steps = 10;
%             [L, weight] = GammaRadByVolList(a, b, intR_steps);
            [L, weight0] = GammaRadList(a, b, intR_steps);
            weight = weight0.*(L*1E6).^3;
            weight = weight./sum(weight);
         else error Unknown GammaRdist option

        end
    else
         L = protocol.GammaR;
         weight = protocol.GammaWeights;   
    end
   
      fun_name = ['CuboidSym_' approx '_' protocol.pulseseq]; 
       fun = str2func(fun_name); 
     if isfield(protocol,'G')
        E0 = zeros(length(protocol.G),intR_steps);
     else
         E0 = zeros(length(protocol.G1),intR_steps);
     end
     
       
       
 
    
    if nargout == 1
        for i = 1:length(L)
        modeltemp.name = 'CuboidSym';
        modeltemp.params = [x(1) L(i) ecc];
        protocol = MMConstants(modeltemp,protocol);        

        E0(:,i) = fun(modeltemp.params,protocol);
        end
         weight_mat = repmat(weight,size(E0,1),1);
         E = sum(E0.*weight_mat,2);
             
    else       
        error('Not yet implemented')
   end
   
   
end
% 
% fun_name = ['GammaCylinders_' approx '_' protocol.pulseseq];
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