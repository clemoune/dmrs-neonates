% This script generates the example waveforms from  from Drobnjak et al 2011

% add the path to GENSym directory
parentpath = cd(cd('..'));
addpath(genpath(parentpath));
 

%% dPFG sequences - Fig 1b

clear all

Npoints = 25;
protocoldPFG1.pulseseq = 'dPFG';
protocoldPFG1.smalldel = [repmat(0.0015,1,Npoints) repmat(0.0045,1,Npoints) repmat(0.0075,1,Npoints)];
protocoldPFG1.G = [repmat(0.3757,1,Npoints) repmat(0.1252,1,Npoints) repmat(0.0751,1,Npoints)];
protocoldPFG1.delta = repmat(0.04,1,Npoints*3);
protocoldPFG1.theta = repmat(pi/2,1,Npoints*3);
protocoldPFG1.phi = repmat(0,1,Npoints*3); % first gradient along x direction
protocoldPFG1.theta1 = repmat(pi/2,1,Npoints*3);
protocoldPFG1.phi1 = repmat(linspace(0,2*pi,Npoints),1,3); % second gradient rotating in x-y plane
protocoldPFG1.tm = repmat(0,1,Npoints*3);
tau = 1E-4;
protocoldPFG1.tau = tau;

% create the GEN protocol:
protocolGEN.pulseseq = 'GEN';
protocolGEN.G = wave_form(protocoldPFG1);
protocolGEN.tau = tau;

% display the gradient waveform for the first puslse sequence:
figure();
G_plot = protocolGEN.G(1,1:3:end);
plot((0:tau:(length(G_plot)-1)*tau)*1E3,G_plot*1000,'LineWidth',2)
 xlim([0,(length(G_plot)+5)*tau*1E3])
 ylim([min(G_plot)*1200 max(G_plot)*1200])
 set(gca,'FontSize',16);
 xlabel('time (ms)','FontSize',16);
 ylabel('Gx (mT/m)','FontSize',16);
 title('Gradient wavefor for \phi = 0')
 
% make the model 
model.name = 'Cylinder';
di = 2E-9; % intrinsic diffusivity
rad = 5.2E-6; % cylinder radius
% angles in spherical coordinates describing the cylinder orientation; 
theta = 0; % angle from z axis
phi = 0; % azimuthal angle
model.params = [ di rad theta phi];

% add the matrices and other constants necessary for the matrix method
protocolGEN = MMConstants(model,protocolGEN);

% compute diffusion signal for the given protocol and model
signal = SynthMeas(model,protocolGEN);
%
% plot the signal as a function phi1
signal_plot = reshape(signal,Npoints,3);
figure();
hold on;  
plot(linspace(0,2*pi,Npoints),signal_plot(:,1),'xr',...
    linspace(0,2*pi,Npoints),signal_plot(:,2),'sb',...
    linspace(0,2*pi,Npoints),signal_plot(:,3),'<g',...
    'LineWidth',2,'MarkerSize',12);
xlabel('\phi','FontSize',16);
ylabel('Diffusion Signal','FontSize',16);
legend('\delta = 1.5 ms','\delta = 4.5 ms','\delta = 7.5 ms')
set(gca,'FontSize',16);
title('Diffusion signal for dPFG (Fig. 1b)')

%% dPFG sequences - Figure 1c

tau = 1E-4;
protocoldPFG2.pulseseq = 'dPFG';
protocoldPFG2.smalldel = repmat(0.0015,1,Npoints*3);
protocoldPFG2.G = [repmat(0.1377,1,Npoints) repmat(0.2567,1,Npoints) repmat(0.3757,1,Npoints)];
protocoldPFG2.delta = repmat(0.12,1,Npoints*3);
protocoldPFG2.theta = repmat(pi/2,1,Npoints*3);
protocoldPFG2.phi = repmat(0,1,Npoints*3); % first gradient along x direction
protocoldPFG2.theta1 = repmat(pi/2,1,Npoints*3);
protocoldPFG2.phi1 = repmat(linspace(0,2*pi,Npoints),1,3); % second gradient rotating in x-y plane
protocoldPFG2.tm = repmat(0,1,Npoints*3);
protocoldPFG2.tau = tau;

% create the GEN protocol:
protocolGEN.pulseseq = 'GEN';
protocolGEN.G = wave_form(protocoldPFG2);
protocolGEN.tau = tau;

% establish which tissue model will be used; all the available options are
% listed in the documentation The name of the models are consistent with
% the ones in Panagiotaki et al, 2012.

% an example for a white matter model:

model.name = 'Cylinder';
di = 2E-9; % intrinsic diffusivity
rad = 9.7E-6; % cylinder radius
% angles in spherical coordinates describing the cylinder orientation; 
theta = 0; % angle from z axis
phi = 0; % azimuthal angle
model.params = [ di rad theta phi];

% add the matrices and other constants necessary for the matrix method
protocolGEN = MMConstants(model,protocolGEN);

% compute diffusion signal for the given protocol and model
signal = SynthMeas(model,protocolGEN);
%
% plot the signal as a function phi1
signal_plot = reshape(signal,Npoints,3);
figure();

hold on;  
plot(linspace(0,2*pi,Npoints),signal_plot(:,1),'xr',...
    linspace(0,2*pi,Npoints),signal_plot(:,2),'sb',...
    linspace(0,2*pi,Npoints),signal_plot(:,3),'<g',...
    'LineWidth',2);

xlabel('\phi','FontSize',16);
ylabel('Diffusion Signal','FontSize',16);
legend('G = 137.7 mT/m','G = 256.7 mT/m','G = 375.7 mT/m')
set(gca,'FontSize',16);
title('Diffusion signal for dPFG (Fig. 1c)')


% display the gradient waveform for the first puslse sequence:
figure();
G_plot = protocolGEN.G(1,1:3:end);
plot((0:tau:(length(G_plot)-1)*tau)*1E3,G_plot*1000,'LineWidth',2)
 xlim([0,(length(G_plot)+5)*tau*1E3])
 ylim([min(G_plot)*1200 max(G_plot)*1200])
 xlabel('time (ms)');
 ylabel('Gx (mT/m)'); 
 
 
%% dPFG sequences - Figure 1d

protocoldPFG3.pulseseq = 'dPFG';
protocoldPFG3.smalldel = repmat(0.0015,1,Npoints*4) ;
protocoldPFG3.G = repmat(0.3757,1,Npoints*4) ;
protocoldPFG3.delta = repmat(0.04,1,Npoints*4);
protocoldPFG3.theta = repmat(pi/2,1,Npoints*4);
protocoldPFG3.phi = repmat(0,1,Npoints*4); % first gradient along x direction
protocoldPFG3.theta1 = repmat(pi/2,1,Npoints*4);
protocoldPFG3.phi1 = repmat(linspace(0,2*pi,Npoints),1,4); % second gradient rotating in x-y plane
protocoldPFG3.tm = [repmat(0,1,Npoints) repmat(0.005,1,Npoints) repmat(0.02,1,Npoints) repmat(0.1,1,Npoints)] ;
protocoldPFG3.tau = tau;

% create the GEN protocol:
protocolGEN.pulseseq = 'GEN';
protocolGEN.G = wave_form(protocoldPFG3);
protocolGEN.tau = tau;

% establish which tissue model will be used; all the available options are
% listed in the documentation The name of the models are consistent with
% the ones in Panagiotaki et al, 2012.

% an example for a white matter model:

model.name = 'Cylinder';
di = 2E-9; % intrinsic diffusivity
rad = 5.2E-6; % cylinder radius
% angles in spherical coordinates describing the cylinder orientation; 
theta = 0; % angle from z axis
phi = 0; % azimuthal angle
model.params = [ di rad theta phi];

% add the matrices and other constants necessary for the matrix method
protocolGEN = MMConstants(model,protocolGEN);

% compute diffusion signal for the given protocol and model
signal = SynthMeas(model,protocolGEN);
%
% plot the signal as a function phi1
signal_plot = reshape(signal,Npoints,4);
figure();

hold on;  
plot(linspace(0,2*pi,Npoints),signal_plot(:,1),'xr',...
    linspace(0,2*pi,Npoints),signal_plot(:,2),'sb',...
    linspace(0,2*pi,Npoints),signal_plot(:,3),'<g',...
    linspace(0,2*pi,Npoints),signal_plot(:,3),'^k',...
    'LineWidth',2);

xlabel('\phi','FontSize',16);
ylabel('Diffusion Signal','FontSize',16);
legend('\tau_m = 0 ms','\tau_m = 5 ms','\tau_m = 20 ms','\tau_m = 1000 ms')
set(gca,'FontSize',16);
title('Diffusion signal for dPFG (Fig. 1b)')


% display the gradient waveform for the first puslse sequence:
figure();
G_plot = protocolGEN.G(1,1:3:end);
plot((0:tau:(length(G_plot)-1)*tau)*1E3,G_plot*1000,'LineWidth',2)
 xlim([0,(length(G_plot)+5)*tau*1E3])
 ylim([min(G_plot)*1200 max(G_plot)*1200])
 xlabel('time (ms)');
 ylabel('Gx (mT/m)'); 

%% Helical waveforms - Figure 2

clear all
fosc = 0.2:0.2:3; % normalized frequency
smalldel = 0.01; % duration
delta = 0.0127;
protocolHelical.pulseseq = 'Helical';
protocolHelical.smalldel = repmat(smalldel,1,length(fosc)*3) ;
protocolHelical.G = [repmat(0.05,1,length(fosc)) repmat(0.08,1,length(fosc)) repmat(0.15,1,length(fosc))];
protocolHelical.delta = repmat(delta,1,length(fosc)*3);
protocolHelical.omega = repmat(2*pi*fosc/smalldel,1,3);
protocolHelical.slopez = repmat(1/5/smalldel,1,length(fosc)*3);
tau = 1E-4;
protocolHelical.tau = tau;

% create the GEN protocol:
protocolGEN.pulseseq = 'GEN';
protocolGEN.G = wave_form(protocolHelical);
protocolGEN.tau = tau;

% display the gradient waveform for the first puslse sequence:
figure();
G_plot = protocolGEN.G(end,1:3:end);
plot((0:tau:(length(G_plot)-1)*tau)*1E3,G_plot*1000,'LineWidth',2)
 xlim([0,(length(G_plot)+5)*tau*1E3])
 ylim([min(G_plot)*1200 max(G_plot)*1200])
 xlabel('time (ms)','FontSize',16);
 ylabel('Gx (mT/m)','FontSize',16);
 set(gca,'FontSize',16);
 title('Gradient waveform in x direction with the highest frequency','FontSize',16);

model.name = 'Cylinder';
di = 2E-9; % intrinsic diffusivity
rad = 5E-6; % cylinder radius
% angles in spherical coordinates describing the cylinder orientation; 
theta = 0; % angle from z axis
phi = 0; % azimuthal angle
model.params = [ di rad theta phi];

% add the matrices and other constants necessary for the matrix method
protocolGEN = MMConstants(model,protocolGEN);

% compute diffusion signal for the given protocol and model
signal = SynthMeas(model,protocolGEN);
%
% plot the signal as a function phi1
signal_plot = reshape(signal,length(fosc),3);
figure();

hold on;  
plot(fosc,signal_plot(:,1),'xr',...
    fosc,signal_plot(:,2),'sb',...
    fosc,signal_plot(:,3),'<g',...   
    'LineWidth',2,'MarkerSize',12);
xlabel('f','FontSize',16);
ylabel('Diffusion Signal','FontSize',16);
legend('G = 50 mT/m','G = 80 mT/m','G = 150 mT/m')
set(gca,'FontSize',16);
title('Diffusion signal for Helical waveform (Fig. 2)')




%% STEAM sequences - Fig 3 

clear all

protocolSTEAM.pulseseq = 'STEAM';
protocolSTEAM.smalldel = repmat(0.0079,1,50);
protocolSTEAM.sdels = repmat(0.001,1,50);
protocolSTEAM.sdelc = repmat(0.0015,1,50);
protocolSTEAM.G = repmat(0.14,1,50);
protocolSTEAM.Gs = repmat(0.139,1,50);
protocolSTEAM.Gc = repmat(0.04,1,50);
protocolSTEAM.grad_dirs = abs(ReadCaminoElecPS(sprintf('PointSets/Elec%03i.txt',50)));
protocolSTEAM.tau1 = repmat(0,1,50);
protocolSTEAM.tau2 = repmat(0,1,50);
protocolSTEAM.tm = repmat(0.01,1,50);
tau = 1E-4;
protocolSTEAM.tau = tau;

% create the GEN protocol:
protocolGEN.pulseseq = 'GEN';
protocolGEN.G = wave_form(protocolSTEAM);
protocolGEN.tau = tau;

% display the gradient waveform for the first puslse sequence:
figure();
Gz = protocolGEN.G(1,3:3:end);
plot((0:tau:(length(Gz)-1)*tau)*1E3,Gz*1000,'LineWidth',2)
 xlim([0,(length(Gz)+5)*tau*1E3])
 ylim([min(Gz)*1200 max(Gz)*1200])
 xlabel('time (ms)','FontSize',16);
 ylabel('Gz (mT/m)','FontSize',16);
 set(gca,'FontSize',16);

% make te tissue model
model.name = 'Cylinder';
di = 0.7E-9; % intrinsic diffusivity
rad = 5E-6; % cylinder radius
% angles in spherical coordinates describing the cylinder orientation; 
theta = 0; % angle from z axis
phi = 0; % azimuthal angle
model.params = [ di rad theta phi];

% add the matrices and other constants necessary for the matrix method
protocolGEN = MMConstants(model,protocolGEN);

% compute diffusion signal for the given protocol and model
signal = SynthMeas(model,protocolGEN);
%
% plot the signal as a function |Gz|/G
signal_plot = signal;
Gplot = abs(protocolSTEAM.grad_dirs(:,3));
figure();
plot(Gplot,signal_plot,'x','LineWidth',2,'MarkerSize',12);
xlabel('|G_z|/G','FontSize',16);
ylabel('Diffusion Signal','FontSize',16);
set(gca,'FontSize',16);
title('Diffusion signal for STEAM')



