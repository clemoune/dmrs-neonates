% This script sets up two example GEN protocols and computes the 
% corresponding diffusion MRI signal

% the gradient waveform is a M x3K matrix, where M is the number of measurements and
% K is the number of gradient points in one measurement; 
% the gradient waveform has the following structure:
% G1x(1) G1y(1) G1z(1) G1x(2) G1y(2) G1z(2) .... G1x(K) G1y(K) G1z(K)
%  ........
% GMx(1) GMy(1) GMz(1) GMx(2) GMy(2) GMz(2) .... GMx(K) GMy(K) GMz(K)
% 
% The signal computation is very flexible, however, for a diffusion MRI
% experiment, the integral of the gradient at the echo time must be 0
% 
% We provide a set of example gradient waveforms including pulse gradient 
% spin echo sequence (PGSE), sinusoidal oscillating gradients (OGSE), 
% square oscillating gradients (SWOGSE), trapezoidal oscillating gradients (TWOGSE),
% square oscillating gradients with different parameters in x,y and z directions (SWOGSE_3D), 
% double pulsed field gradients (dPFG), stimulated echo sequences (STEAM), 
% helical gradients (Helical). See documentation for details related to the
% parameters of each sequence.
% 
% The discretized gradients for the example waveforms are created by the
% wave_form(protocol) function, where protocol contains information about the desired sequence 
% 

% add the path to MISST directory
parentpath = cd(cd('..'));
addpath(genpath(parentpath));
%% Example 1 - square wave oscillating gradients (SWOGSE)

% setup an initial protocol which includes the necessary information 
% to compute the discretized waveform for a SWOGSE sequence
clear all
% add some combinations of parameters
delta = 0.015:0.005:0.04; % duration in s
smalldel = delta - 0.005;  % duration in s
Nosc = 1:5;  % number of lobes in the oscillating gradient. A full period has 2 lobes
Gmax = 0.08; % gradient strength in T/m;
tau = 1E-4; % time interval for waveform discretization
it = 0;
protocol_init.pulseseq = 'SWOGSE';
for i = 1:length(delta)
    for j = 1:length(Nosc)
        it = it+1;
        protocol_init.smalldel(it) = smalldel(i);
        protocol_init.delta(it) = delta(i);
        protocol_init.omega(it) = Nosc(j).*pi./smalldel(i);
        protocol_init.G(it) = Gmax;
        protocol_init.grad_dirs(it,:) = [1 0 0]; % gradient in x direction
    end
end
protocol_init.tau = tau;

% create the GEN protocol:
protocolGEN.pulseseq = 'GEN';
protocolGEN.G = wave_form(protocol_init);
protocolGEN.tau = tau;
% include smalldel and delta as they make computation slightly faster
protocolGEN.delta = protocol_init.delta;
protocolGEN.smalldel = protocol_init.smalldel;


% Add an example for a white matter tissue model:

model.name = 'FiniteAstroCylinders';
di = 1.7E-9; % intrinsic diffusivity
rad = 4E-6; % cylinder radius
lcyl = 16E-6;
model.params = [ di rad lcyl];

% add the matrices and other constants necessary for the matrix method
protocolGEN = MMConstants(model,protocolGEN);

% display the gradient waveform for the first puslse sequence:
figure();
Gx = protocolGEN.G(1,1:3:end);
plot((0:tau:(length(Gx)-1)*tau)*1E3,Gx*1000,'LineWidth',2)
xlim([0,(protocolGEN.smalldel(1)+protocolGEN.delta(1)+10*tau)*1E3])
ylim([min(Gx)*1200 max(Gx)*1200])
xlabel('time (ms)','FontSize',16);
ylabel('Gx (mT/m)','FontSize',16);
set(gca,'FontSize',16);

% compute diffusion signal for the given protocol and model
signal = SynthMeas(model,protocolGEN);

% plot the signal as a function of delta for Nosc = 1
signal_matrix = reshape(signal,length(Nosc),length(delta));
figure();
colors = jet(length(Nosc));
hold on;
for i  = 1:length(Nosc)  
    legend_name{i} = ['Nosc = ' num2str(Nosc(i))];
    plot(delta*1000,signal_matrix(i,:),'LineWidth',2,'Color',colors(i,:));
end

xlabel('\Delta (ms)','FontSize',16);
ylabel('Diffusion Signal','FontSize',16);
set(gca,'FontSize',16);
title('Diffusion signal as a function of \Delta for various number of oscillations, G = 80mT/m')
legend(legend_name,'FontSize',16)




%% Example 2: random gradient waveform. 
clear all
smalldel = 0.02; % duration of the waveform;
Npoints = 100; % 100 random points
tau = smalldel / Npoints; % sampling interval
p180 = 0.005; % time interval between the first and second waveforms (duration of rf pulse)
Gmax = 0.08; % % maximum gradient in T/m; in T/m;

% the gradient waveforms will be repeated after the 180 pulse
wf = [rand(1,Npoints)*Gmax; rand(1,Npoints)*Gmax; rand(1,Npoints)*Gmax]; 
% gradient waveform consists of random points along x, y and z  
% create protocolGEN in the necessary format
protocolGEN.G = [zeros(1,3) wf(:)' zeros(1,3* round(p180/tau)) -wf(:)' zeros(1,3)]; 
protocolGEN.tau=tau;
protocolGEN.pulseseq = 'GEN';

% display the one vector of random numbers and the resulting gradient waveform
figure();
subplot(1,2,1)
plot(wf(1,:),'LineWidth',2)
title('Initial vector with random points (for x direction)','FontSize',16);
set(gca,'FontSize',16);
subplot(1,2,2)
Gx = protocolGEN.G(1,1:3:end);
plot((0:protocolGEN.tau:(length(Gx)-1)*protocolGEN.tau)*1E3,Gx*1000,'LineWidth',2)
xlabel('time (ms)','FontSize',16);
ylabel('Gx (mT/m)','FontSize',16);
title('Generated diffusion gradient along x direction','FontSize',16);
set(gca,'FontSize',16);
xlim([0,(length(Gx)+1)*tau*1E3])

% Add an example for a white matter tissue model:

model.name = 'ZeppelinCylinder';
di = 1.7E-9; % intrinsic diffusivity
dh = 1.2E-9; % hindered diffusivity
rad = 4E-6; % cylinder radius
% angles in spherical coordinates describing the cylinder orientation; 
theta = 0; % angle from z axis
phi = 0; % azimuthal angle
ficvf = 0.7; % intracellular volume fraction
model.params = [ficvf di dh rad theta phi];

% add the matrices and other constants necessary for the matrix method
protocolGEN = MMConstants(model,protocolGEN);

% compute diffusion signal for the given protocol and model
signal = SynthMeas(model,protocolGEN);

