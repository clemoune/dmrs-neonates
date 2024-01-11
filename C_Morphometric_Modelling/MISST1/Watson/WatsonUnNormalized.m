function f = WatsonUnNormalized(kappa,x)
% function f = WatsonUnNormalized(kappa,x)
%
% Compute the Watson's distribution without the normalization coefficient
%
% $Id: WatsonUnNormalized.m,v 1.2 2009/12/12 23:42:54 gzhang Exp $

f = exp(kappa*x.^2);

