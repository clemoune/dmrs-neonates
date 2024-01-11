function E=WatsonCylNeumanComponent(z, phi, Lpmp, kappa, watsonNC, dp)
%
% For Internal use ONLY!!!
%
% E=WatsonCylNeumanComponent_PGSE(z, phi, Lpmp, kappa, watsonNC, nu)
%
% $Id: WatsonCylNeumanComponent_PGSE.m,v 1.7 2010/02/02 13:38:35 gzhang Exp $

cosTheta = z;
sinTheta = sqrt(1 - z.*z);
cosPhi = cos(phi);
sinPhi = sin(phi);

% the gradient and fibre direction in the local reference frame of the symmetry
% axis of the Watson's distribution
g = [sqrt(1-dp*dp), 0, dp];
n = [sinTheta*cosPhi; sinTheta*sinPhi; cosTheta];

% compute the angle between the gradient and the fibre direction
w = g*n;

watsonWeighting = watsonNC*WatsonUnNormalized(kappa,cosTheta);
signal = exp(Lpmp*w.^2);

E = signal.*watsonWeighting;
