%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script fit data with the morphometric model at long diffusion times. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Marco Palombo
% Small modifications: Clémence Ligneul 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 


clear all
close all
clc

%% Filter training set

Dintra_bounds = [0.1 1];
% Dintra_bounds = [0.4 0.5];
N_bounds = [2 30];
L_bounds = [5 100];
SDN_bounds = [2 3];
SDL_bounds = [5 10];

MachineLearning_ConstrainedDatabase_ClemenceCerebellum(Dintra_bounds, N_bounds, L_bounds,SDN_bounds,SDL_bounds);

%% Load main dictionary for plotting purposes

load Dictionary.mat

Y_full{1}=Y_dictionary;
Dintra_full{1}=Dintra;
N_full{1}=N;
L_full{1}=L;
SDN_full{1}=SDN;
SDL_full{1}=SDL;

load Dictionary_2_partial.mat

Y_full{2}=Y_dictionary;
Dintra_full{2}=Dintra;
N_full{2}=N;
L_full{2}=L;
SDN_full{2}=SDN;
SDL_full{2}=SDL;

load Dictionary_3.mat

Y_full{3}=Y_dictionary;
Dintra_full{3}=Dintra;
N_full{3}=N;
L_full{3}=L;
SDN_full{3}=SDN;
SDL_full{3}=SDL;

load Dictionary_Part_II.mat

Y_full{4}=Y_dictionary;
Dintra_full{4}=Dintra;
N_full{4}=N;
L_full{4}=L;
SDN_full{4}=SDN;
SDL_full{4}=SDL;

load Dictionary_Part_JV.mat

Y_full{5}=Y_dictionary;
Dintra_full{5}=Dintra;
N_full{5}=N;
L_full{5}=L;
SDN_full{5}=SDN;
SDL_full{5}=SDL;

close all

dictionary = struct;
dictionary.Y_full = Y_full;
dictionary.Dintra_full = Dintra_full;
dictionary.N_full = N_full;
dictionary.L_full = L_full;
dictionary.SDN_full = SDN_full;
dictionary.SDL_full = SDL_full;

%% Train Random Forest

disp('Training the Random Forest Regressor')

tic

load RF_features_extended_ClemenceCerebellum_constrained.mat
load param_to_fit_extended_ClemenceCerebellum_constrained.mat

features = {'Dintra','Nbranch', 'Lbranch', 'SDNbranch', 'SDLbranch'};

% Expand the dictionary to include spherical compartment

fsphere = linspace(0.01,0.95, 6);
Rsphere = linspace(3, 15, 6);
Rcyl = linspace(0.25, 3, 6);

% Create the dictionary for sphere + randomly oriented cylinders



rng(1); % For reproducibility
Mdl = cell(numel(features),1);
Ntree = 100;
tree_depth = [];

for i =1:numel(features)
    Mdl{i} = TreeBagger(Ntree,RF_features_reduced,param_to_fit(:,i), 'Method','regression', 'MaxNumSplits', tree_depth, 'Surrogate','on','PredictorSelection','curvature','OOBPredictorImportance','on');
end

ttt = toc;

disp(['[DONE] - ' num2str(round(ttt)) ' sec.'])
%% Mouse Data 

time_points = {'P5', 'P10', 'P15', 'P20', 'P30'};
datafolder = '~/MATLAB/code_231005/Data';
plot_flag = 0;

for time = 1:numel(time_points)

disp(['Processing Cerebellum and Thalamus data for time point ' time_points{time}])

tic

eval(['load ' datafolder '/' time_points{time} '_230825']);

exp_Td = [100, 500, 750, 1000];

metab = {'tNAA', 'tCr', 'tCho', 'Glu', 'Ins', 'Tau'};
morphometry = struct;

for i=1:numel(metab)
    
    disp([' - Cerebellum: processing metabolite: ' metab{i} ' [ ' num2str(i) '/' num2str(numel(metab)) ' ]'])

    % Cerebellum
    
Y = Age.Cereb(i).ADCCorrectedmean;

quartiles = zeros(numel(features), 3);

for j=1:numel(features)
    mpgMean = predict(Mdl{j},Y);
    morphometry.Cereb(i).RFregression_mean(j) = mpgMean;
    mpgQuartiles = quantilePredict(Mdl{j},Y,'Quantile',[0.25,0.5,0.75]);
    quartiles(j,:) = mpgQuartiles;
    
    morphometry.Cereb(i).RFregression_median(j) = quartiles(j,2);
    morphometry.Cereb(i).RFregression_25thPercentile(j) = quartiles(j,1);
    morphometry.Cereb(i).RFregression_75thPercentile(j) = quartiles(j,3);
    morphometry.Cereb(i).feature{j} = features{j};
end

if plot_flag==1, figure('Name', [metab{i} ' in Cerebellum at ' time_points{time}]), end

morphometry.Cereb(i).exp_signal = Y;
morphometry.Cereb(i).exp_Td = exp_Td;

morphometry.Cereb(i).pred_signal_mean = plot_fit_ClemenceCerebellum(Y, exp_Td, morphometry.Cereb(i).RFregression_mean, dictionary,plot_flag, 'r-');
morphometry.Cereb(i).pred_signal_25thPercentile = plot_fit_ClemenceCerebellum(Y, exp_Td, quartiles(:,1)', dictionary,plot_flag, 'k--');
morphometry.Cereb(i).pred_signal_median = plot_fit_ClemenceCerebellum(Y, exp_Td, quartiles(:,2)', dictionary,plot_flag, 'k-');
morphometry.Cereb(i).pred_signal_75thPercentile = plot_fit_ClemenceCerebellum(Y, exp_Td, quartiles(:,3)', dictionary,plot_flag, 'k--');

 % Thalamus
 
     disp([' - Thalamus: processing metabolite: ' metab{i} ' [ ' num2str(i) '/' num2str(numel(metab)) ' ]'])
 
Y = Age.Thalam(i).ADCCorrectedmean;

quartiles = zeros(numel(features), 3);

for j=1:numel(features)
    
    mpgMean = predict(Mdl{j},Y);
    morphometry.Thalam(i).RFregression_mean(j) = mpgMean;
    mpgQuartiles = quantilePredict(Mdl{j},Y,'Quantile',[0.25,0.5,0.75]);
    quartiles(j,:) = mpgQuartiles;

    morphometry.Thalam(i).RFregression_median(j) = quartiles(j,2);
    morphometry.Thalam(i).RFregression_25thPercentile(j) = quartiles(j,1);
    morphometry.Thalam(i).RFregression_75thPercentile(j) = quartiles(j,3);
    morphometry.Thalam(i).feature{j} = features{j};
end

if plot_flag==1, figure('Name', [metab{i} ' in Thalamus at ' time_points{time}]), end

morphometry.Thalam(i).exp_signal = Y;
morphometry.Thalam(i).exp_Td = exp_Td;

morphometry.Thalam(i).pred_signal_mean = plot_fit_ClemenceCerebellum(Y, exp_Td, morphometry.Thalam(i).RFregression_mean, dictionary,plot_flag, 'r-');
morphometry.Thalam(i).pred_signal_25thPercentile = plot_fit_ClemenceCerebellum(Y, exp_Td, quartiles(:,1)', dictionary,plot_flag, 'k--');
morphometry.Thalam(i).pred_signal_median = plot_fit_ClemenceCerebellum(Y, exp_Td, quartiles(:,2)', dictionary,plot_flag, 'k-');
morphometry.Thalam(i).pred_signal_75thPercentile = plot_fit_ClemenceCerebellum(Y, exp_Td, quartiles(:,3)', dictionary,plot_flag, 'k--');

end

save([datafolder '/RF_features_Raw_reduced_ClemenceCerebellum_OtherNandL_' time_points{time}], 'morphometry');

ttt = toc;

disp(['[DONE] - ' num2str(round(ttt)) ' sec.'])

end

