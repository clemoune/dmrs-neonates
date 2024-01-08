%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script can be used to co-fit the high b values data and the long
% diffusion times data. It fits together the morphometric moedel and the
% astrosticks & spheres model. However there is some redundancy in the
% model and it does not fit well data together. Therefore we chose to model
% the morphometric model alone on long diffusion times. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Marco Palombo
% Small modifications: Clémence Ligneul 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

clear all
close all
clc

%%% The MISST toolbox is used to 
addpath(genpath('~/MATLAB/code_231005/MachineLearningPipeline/MISST_1'))


%% Filter training set

Dintra_bounds = [0.01 1];
N_bounds = [2 20];
L_bounds = [20 500];
SDN_bounds = [2 3];
SDL_bounds = [5 10];

MachineLearning_ConstrainedDatabase_ClemenceCerebellum(Dintra_bounds, N_bounds, L_bounds, SDN_bounds, SDL_bounds);

%% Train Random Forest

disp('Training the Random Forest Regressor')

load RF_features_extended_ClemenceCerebellum_constrained.mat
load param_to_fit_extended_ClemenceCerebellum_constrained.mat

%% Expand the dictionary to include spherical compartment

Nset = 100;
fsphere = rand(Nset,1);
Rsphere = rand(Nset,1).*(15-3) + 3;
%Rcyl = rand(Nset,1).*(3-0.25) + 0.25;

%new_params = [fsphere, Rsphere, Rcyl];
new_params = [fsphere, Rsphere];

% For S(b) 
bvals = [0.035, 3.035, 6, 10, 20, 30];
smalldelta = 5;
delta = 100 + smalldelta/3;

%fmodel_b = @(p,x) p(1).*SynthMeasSphere([0.5.*1e-9, p(2).*1e-6], x) + (1 - p(1)).*SynthMeasAstroCylinders([p(4).*1e-9, p(3).*1e-6], x); 
fmodel_b = @(p,x) p(1).*SynthMeasSphere([0.5.*1e-9, p(2).*1e-6], x) + (1 - p(1)).*SynthMeasAstroSticks(p(3).*1e-9, x); 

protocol = make_protocol(bvals, delta, smalldelta);
protocol.roots_sphere = BesselJ_RootsSphere(100);
protocol.roots_cyl = BesselJ_RootsCyl(100);
protocol.pulseseq='GPD_PGSE';

new_training_set = zeros(size(new_params,1)*size(RF_features_reduced,1), numel(bvals));
new_training_params = zeros(size(new_params,1)*size(RF_features_reduced,1), size(param_to_fit, 2) + size(new_params, 2));

tic
h = waitbar(0,'Expanding Training Set');
k = 0;
for i = 1:size(new_params,1)
    waitbar(i/size(new_params,1))

parfor j = 1:size(RF_features_reduced,1)
    
    new_training_params(k+j,:) = [new_params(i,:), param_to_fit(j,:)];
    new_training_set(k+j,:) = fmodel_b([new_params(i,:), param_to_fit(j,1)], protocol);
    
end
k = k+size(RF_features_reduced,1);
end
close(h);

toc

%% For ADC(td)
td = [100, 100, 500, 500, 750, 750, 1000, 1000];
delta_td = td + smalldelta/3;
b0s =[0.035, 0.170, 0.260, 0.350];
bs = b0s + 3;
bvals_td = [b0s', bs']';

B = (bvals_td(:))';

GAMMA = 2.675987E8;
protocol = struct;
protocol.pulseseq = 'PGSE';
protocol.testrategy = 'fixed';
protocol.roots_sphere = BesselJ_RootsSphere(100);
protocol.roots_cyl = BesselJ_RootsCyl(100);

protocol.delta = delta_td.*1e-3;
protocol.smalldel = smalldelta.*1e-3.*ones(size(B));
protocol.G = sqrt(B.*1e9./(protocol.delta-protocol.smalldel/3))./(GAMMA.*protocol.smalldel);
protocol.grad_dirs = repmat([1, 1, 1]./sqrt(3), [numel(B), 1]);

new_training_set2 = zeros(size(new_params,1)*size(RF_features_reduced,1), numel(td)/2);

tic
h = waitbar(0,'Expanding Training Set');
k = 0;
for i = 1:size(new_params,1)
    waitbar(i/size(new_params,1))

parfor j = 1:size(RF_features_reduced,1)
    
    new_training_set2(k+j,:) = new_params(i,1).*compute_ADC_from_Sphere([0.5*1E-9, new_params(i,2)*1e-6], protocol) + (1-new_params(i,1)).*RF_features_reduced(j,:);
    
end
k = k+size(RF_features_reduced,1);
end
close(h);
toc

%%

training_dictionary = [new_training_set, new_training_set2];
training_params = new_training_params;
%features = {'fsph', 'Rsph,','Rcyl','Dintra','Nbranch', 'Lbranch', 'SDNbranch', 'SDLbranch'};
features = {'fsph', 'Rsph,','Dintra','Nbranch', 'Lbranch', 'SDNbranch', 'SDLbranch'};

%% Create the dictionary for sphere + randomly oriented cylinders

rng(1); % For reproducibility
Mdl = cell(size(training_params, 2),1);

tic

%RF
Ntree = 100;

h = waitbar(0,'Training RF');
for i =1:size(training_params, 2)
    waitbar(i/size(training_params, 2))
    Mdl{i} = TreeBagger(Ntree,training_dictionary,training_params(:,i), 'Method','regression','OOBPredictorImportance','on');
end
close(h)

% MLP
% training_performances = cell(size(training_params, 2),1);
% net_structure = 10;
% 
% % h = waitbar(0,'Training MLP');
% for i = 1:size(training_params,2)
%     Mdl{i} = feedforwardnet(net_structure);
%     Mdl{i}.trainParam.showWindow = false;
%     Mdl{i}.trainParam.showCommandLine = false;
% end
% 
% parfor i=1:size(training_params,2)
% %     waitbar(i/size(training_params, 2))
%     [Mdl{i}, training_performances{i}] = train(Mdl{i}, training_dictionary', training_params(:,i)');
% 
% end
% % close(h);

ttt = toc;

disp(['[DONE] - ' num2str(round(ttt)) ' sec.'])
%% Mouse Data 

time_points = {'P5', 'P10', 'P15', 'P20', 'P30'};
datafolder = '~/MATLAB/code_231005/Data';
plot_flag = 0;

for time = 1:numel(time_points)

disp(['Processing Cerebellum and Thalamus data for time point ' time_points{time}])

tic

eval(['load ' datafolder '/' time_points{time} '_211101']);

exp_Td = [100, 500, 750, 1000];

metab = {'tNAA', 'tCr', 'tCho', 'Glu', 'Ins', 'Tau'};
morphometry = struct;

for i=1:numel(metab)
    
    disp([' - Cerebellum: processing metabolite: ' metab{i} ' [ ' num2str(i) '/' num2str(numel(metab)) ' ]'])

    % Cerebellum
    
Y = [exp(Age.Cereb(i).mean)' Age.Cereb(i).ADCmean];

quartiles = zeros(numel(features), 3);

for j=1:numel(Mdl)
    
    mpgMean = predict(Mdl{j},Y);
    %net = Mdl{j};
    %mpgMean = net(Y');
    morphometry.Cereb(i).RFregression_mean(j) = mpgMean;
    
    %mpgQuartiles = quantilePredict(Mdl{j},Y,'Quantile',[0.25,0.5,0.75]);
    %quartiles(j,:) = mpgQuartiles;
    
    morphometry.Cereb(i).RFregression_median(j) = quartiles(j,2);
    morphometry.Cereb(i).RFregression_25thPercentile(j) = quartiles(j,1);
    morphometry.Cereb(i).RFregression_75thPercentile(j) = quartiles(j,3);
    morphometry.Cereb(i).feature{j} = features{j};
end

if plot_flag==1, figure('Name', [metab{i} ' in Cerebellum at ' time_points{time}]), end

morphometry.Cereb(i).exp_signal = Y;
%morphometry.Cereb(i).exp_Td = exp_Td;

%morphometry.Cereb(i).pred_signal_mean = plot_fit_ClemenceCerebellum(Y, exp_Td, morphometry.Cereb(i).RFregression_mean, dictionary,plot_flag, 'r-');
%morphometry.Cereb(i).pred_signal_25thPercentile = plot_fit_ClemenceCerebellum(Y, exp_Td, quartiles(:,1)', dictionary,plot_flag, 'k--');
%morphometry.Cereb(i).pred_signal_median = plot_fit_ClemenceCerebellum(Y, exp_Td, quartiles(:,2)', dictionary,plot_flag, 'k-');
%morphometry.Cereb(i).pred_signal_75thPercentile = plot_fit_ClemenceCerebellum(Y, exp_Td, quartiles(:,3)', dictionary,plot_flag, 'k--');

 % Thalamus
 
     disp([' - Thalamus: processing metabolite: ' metab{i} ' [ ' num2str(i) '/' num2str(numel(metab)) ' ]'])
 
Y = [exp(Age.Thalam(i).mean)', Age.Thalam(i).ADCmean];

quartiles = zeros(numel(features), 3);

for j=1:numel(features)
    
    mpgMean = predict(Mdl{j},Y);
    %net = Mdl{j};
    %mpgMean = net(Y');
    morphometry.Thalam(i).RFregression_mean(j) = mpgMean;
    
    %mpgQuartiles = quantilePredict(Mdl{j},Y,'Quantile',[0.25,0.5,0.75]);
    %quartiles(j,:) = mpgQuartiles;

    morphometry.Thalam(i).RFregression_median(j) = quartiles(j,2);
    morphometry.Thalam(i).RFregression_25thPercentile(j) = quartiles(j,1);
    morphometry.Thalam(i).RFregression_75thPercentile(j) = quartiles(j,3);
    morphometry.Thalam(i).feature{j} = features{j};
end

if plot_flag==1, figure('Name', [metab{i} ' in Thalamus at ' time_points{time}]), end

morphometry.Thalam(i).exp_signal = Y;
%morphometry.Thalam(i).exp_Td = exp_Td;

%morphometry.Thalam(i).pred_signal_mean = plot_fit_ClemenceCerebellum(Y, exp_Td, morphometry.Thalam(i).RFregression_mean, dictionary,plot_flag, 'r-');
%morphometry.Thalam(i).pred_signal_25thPercentile = plot_fit_ClemenceCerebellum(Y, exp_Td, quartiles(:,1)', dictionary,plot_flag, 'k--');
%morphometry.Thalam(i).pred_signal_median = plot_fit_ClemenceCerebellum(Y, exp_Td, quartiles(:,2)', dictionary,plot_flag, 'k-');
%morphometry.Thalam(i).pred_signal_75thPercentile = plot_fit_ClemenceCerebellum(Y, exp_Td, quartiles(:,3)', dictionary,plot_flag, 'k--');

end

save([datafolder '/RF_features_Raw_reduced_ClemenceCerebellum_' time_points{time}], 'morphometry');

ttt = toc;

disp(['[DONE] - ' num2str(round(ttt)) ' sec.'])

end

