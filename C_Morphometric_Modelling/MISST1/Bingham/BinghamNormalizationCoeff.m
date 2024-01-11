function c = BinghamNormalizationCoeff(kappa, beta, theta)
% function c = BinghamNormalizationCoeff(kappa, beta, theta)
%
% Compute the normalization coefficient of Bingham's distribution
%

alpha = 2;

if nargin < 2
    M = mhgOptimalM(kappa);
    c = mhg(M, alpha, 0.5, 1.5, [kappa]);
    return;
elseif nargin < 3
    theta = 0;
end

Y = sort(real([kappa beta theta]));
Baseline = Y(1);
kappa = Y(3) - Baseline;
beta = Y(2) - Baseline;
M = mhgOptimalM(kappa);
c = mhg([M], alpha, 0.5, 1.5, [kappa beta]);
c = c*exp(Baseline);

