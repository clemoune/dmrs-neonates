function f=MMnorm(roots,a,model) 
%
% camino.m--------------------------------------------------------------
% Normalization coefficients for the diffusion eigenfunctions
% 
% f=MMnorm(roots,a,model)
% 
% Description: Returns the normalization coefficients of the diffusion
% eigenfunction for a given restriction
%
% Parameters: 
% f - normalization coefficient 
% roots - stores information about the order of the eigenfunction and the
%       order of the root, needed for computing the normalization
%       coefficient
% a - size of restriction (radius for sphere and cylinder restriction,
% half the distance between plates for planar restriction)
% model - 0 = sphere, 1 = cylinder, 2 = parallel planes 
%
%------------------------------------------------------------------------
% This file is part of the camino.m toolbox.
% Copyright (c) 2015, UCL Microstructure Imaging Group (MIG), All rights reserved.
% Distributed under the Modified BSD Licence (see: LICENSE.pdf).
%
% Authors:
%   Ivana Drobnjak (i.drobnjak@ucl.ac.uk)
%   Andrada Ianus (a.ianus.11@ucl.ac.uk)


% 0 is for sphere
% 1 is for cylinder
% 2 is for planes
SMALL=1E-5;
if model==0
    n = roots(1);
    m = roots(2); % upper index of the associated legendre polynomials P_n^m
    alpha = roots(4);
    if n ==0 && alpha <SMALL
       f = (a^3*4*pi/3)^(-1/2);        
    else
    f1=(a^3*pi*((spherical_besselj(n,alpha)).^2-spherical_besselj(n-1,alpha).*...
    spherical_besselj(n+1,alpha))).^(-1/2);
    f = f1*(factorial(n+m)*2/(2*n+1)/(factorial(n-m)))^(-1/2);
    end
elseif model==1
    n = roots(1);
    alpha = roots(3);
    if n==0
        k=1;
    else
        k=2;
    end
    tmp1=k/(pi*a^2*(besselj(n,alpha)).^2);
    if n==0 && alpha<SMALL
        tmp2=1;
    else
        tmp2=alpha^2/(alpha^2-n^2);
    end
    f=(tmp1*tmp2)^(1/2);
elseif model==2
    n = roots(1);
    alpha = roots(2);
    if alpha == 0    
    % add the case for eigenval 0
    f = 1./sqrt(a*2);
    else
        f=sqrt(1./(a*(1+(-1)^n*sin(2*alpha)/(2*alpha))));
    end
else
    error('Model index not specified. 0 for sphere, 1 for cylinder, 2 for planar.')
end
    