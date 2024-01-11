function fanningdir = computeFanningOrientation(theta, phi, psi)
% compute fanning orientation from euler angles

cosTheta = cos(theta);
sinTheta = sin(theta);
cosPhi = cos(phi);
sinPhi = sin(phi);
cosPsi = cos(psi);
sinPsi = sin(psi);
fanningdir = [ -cosPsi.*sinPhi - cosTheta.*cosPhi.*sinPsi cosPhi.*cosPsi - cosTheta.*sinPhi.*sinPsi sinTheta.*sinPsi]';

