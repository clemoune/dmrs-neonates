function c = WatsonNormalizationCoeff(kappa)
%
% camino.m--------------------------------------------------------------
% function c = WatsonNormalizationCoeff(kappa)
%
% Compute the normalization coefficient of Watson's distribution
%
% The accuracy for large kappa is limited by the accuracy of erfi.
%
% $Id: WatsonNormalizationCoeff.m,v 1.2 2009/12/12 23:42:54 gzhang Exp $
%------------------------------------------------------------------------
% This file is part of the camino.m toolbox.
% Copyright (c) 2015, UCL Microstructure Imaging Group (MIG), All rights reserved.
% Distributed under the Modified BSD Licence (see: LICENSE.pdf).
%
% Authors:   
%   Gary Hui Zhang (gary.zhang@ucl.ac.uk)

ksqrt = sqrt(kappa);

if ksqrt < 1E-10;
	c = 1;
else
	c = 2*ksqrt/(sqrt(pi)*erfi(ksqrt));
end

