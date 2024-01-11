function f=spherical_besselj(n,x)
%
% camino.m--------------------------------------------------------------
% Spherical bessel J function
% 
% f=spherical_besselj(n,x)
% 
% Description: Fast implementation of the spherical bessel J function
%
% Parameters:   
% j - evaluated Jn at x
% n - degree, n>=0
% x - variable
%
%------------------------------------------------------------------------
% This file is part of the camino.m toolbox.
% Copyright (c) 2015, UCL Microstructure Imaging Group (MIG), All rights reserved.
% Distributed under the Modified BSD Licence (see: LICENSE.pdf).
%
% Authors:
%   Andrada Ianus (a.ianus.11@ucl.ac.uk)
%
SMALL=1E-20;
f = zeros(size(x));
for k=1:numel(x)
  if x(k)<SMALL
    if n==0
      f(k)=1;
    else
      f(k)=0;
    end
  else
    f(k)=(pi./(2.*x(k))).^(1/2).*besselj(n+1/2,x(k));
  end
end