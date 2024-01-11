function val = associated_legendre(n,m,x)
%
% camino.m--------------------------------------------------------------
% Associated legendre polynomials of degree n and order m
% 
% val = associated_legendre(n,m,x)
% 
% Description: Returns the value of the associated legendre polynomial of
% degree n and order m, evaluate at x. Faster implementation than the
% standard Matlab function which returns the values for all the orders m<=n
%
% Parameters:   
% val - evaluated P_n,m at x
% n - degree, n>=0
% m - order, -n <= m <= n
% x - variable, -1 <= x <= 1
%
%------------------------------------------------------------------------
% This file is part of the camino.m toolbox.
% Copyright (c) 2015, UCL Microstructure Imaging Group (MIG), All rights reserved.
% Distributed under the Modified BSD Licence (see: LICENSE.pdf).
%
% Authors:
%   Andrada Ianus (a.ianus.11@ucl.ac.uk)
%


if ~isscalar(n) || n < 0 || n ~= round(n)
    error('invalid n');
end

if ~isreal(x) || max(abs(x(:))) > 1
    error('x should be between -1 and 1')
end
if ~isscalar(m) || abs(m) >n || m ~= round(m)
    error('invalid m');
end

% recursively calculate the terms
if m>= 0
if n == m
val =  (-1).^m*factorial(2.*m)./factorial(m)./2.^m.*(1-x.^2).^(m./2); % P_m^m
return
end
if n == m+1
term1 = (-1).^m*factorial(2.*m)./factorial(m)./2.^m.*(1-x.^2).^(m./2); % P_m^m
val = x.*(2.*m+1).*term1; % P_m+1^m
return
end

val = 1./(n-m).*(associated_legendre(n-1,m,x).*(2.*n-1).*x - associated_legendre(n-2,m,x).*(n-1+m));
else
val = (-1).^abs(m).*factorial(n-abs(m))./factorial(n+abs(m))*associated_legendre(n,abs(m),x);

end