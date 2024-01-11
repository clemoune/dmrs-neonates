% Finds the roots of the equation J'_1(x) = 0 where J_1 is the first order
% Bessel function of the first kind.
% The argument starts is the number of starting points in the search for
% roots.  The larger it is, the more roots are returned in the root array.
% The default is 20, which returns 6 roots.  starts = 100 returns 32 roots.
% Generally the number of roots is just under 1/3 the number of starting
% points.
function root = BesselJ_RootsSphere(starts)

if(nargin==0)
    starts = 20;
end

f = @(x)(sqrt(2./(pi*(x.^3))).*(2*x.*cos(x) + (x.^2 - 2).*sin(x)));

% Get a good list of starting points
y = 0.00001:0.1:starts;
dj1 = f(y);
starts = (find(dj1(1:end-1).*dj1(2:end)<0)+1)*0.1;

for s=1:length(starts)
    root(s) = fzero(f,starts(s));
end




