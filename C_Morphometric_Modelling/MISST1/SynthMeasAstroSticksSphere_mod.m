function E = SynthMeasAstroSticksSphere_mod(params,protocol)
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

E1 = SynthMeas_mod(model1,protocolGEN1);
E2 = SynthMeas_mod(model2,protocolGEN2);
E = f1*E1+f2*E2;

end