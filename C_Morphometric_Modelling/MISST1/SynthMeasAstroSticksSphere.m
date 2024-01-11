function [E, J] = SynthMeasAstroSticksSphere(params,protocol,params_deriv)
% function [E J] = SynthMeasZeppelinCylinder(model, protocol, params_deriv)
% returns the diffusion signal and signal jacobian for the given protocol
% and ZeppelinCylinder tissue model
%
% author: Marco Palombop (mrc.palombo@gmail.com)

model1.name = 'AstroSticks';
model2.name = 'Sphere';

index1 = 2; % parameters for astrosticks
index2 = [2 4]; % parameters for sphere
f1 = params(1); % volume fraction of astrostick compartment
f2 = params(3); % volume fraction of sphere compartment

model1.params = params(index1);
model2.params = params(index2);

protocolGEN1 = MMConstants(model1,protocol);
protocolGEN1.pulseseq = 'GEN';
protocolGEN1.G = wave_form(protocol);

protocolGEN2 = MMConstants(model2,protocol);
protocolGEN2.pulseseq = 'GEN';
protocolGEN2.G = wave_form(protocol);

if nargout == 1
        E1 = SynthMeas(model1,protocolGEN1);
        E2 = SynthMeas(model2,protocolGEN2);
        E = f1*E1+f2*E2;
else 
    if nargin <3
        params_deriv = ones(size(params));
    end
    [E1, J1] = SynthMeas(model1,protocolGEN1);
    [E2, J2] = SynthMeas(model2,protocolGEN2);
    E = f1*E1+f2*E2;
    J = zeros(length(E),nnz(params_deriv));
    it = 0;
    if params_deriv(1), it = it+1; J(:,it) = E1; end % astrosticks volume fraction
    if params_deriv(2), it = it+1; J(:,it) = f1*J1(:,1) + f2*J2(:,1); end % Dintra
    if params_deriv(3), it = it+1; J(:,it) = E2; end % sphere volume fraction
    if params_deriv(4), it = it+1; J(:,it) = f2*J2(:,2); end % radius    
   
end
end