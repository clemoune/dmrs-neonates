%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script allows to visualise properties coming from the morphometric
% model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Marco Palombo
% Small modifications: Clémence Ligneul 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all
%%

datafolder = '~/MATLAB/code_231005/Data';
% Time points = Age
time_points = {'P5', 'P10', 'P15', 'P20', 'P30'};
T = [5 10 15 20 30];

% Metabolites of interest
metab = {'tNAA', 'tCr', 'tCho', 'Glu', 'Ins', 'Tau'};

% Features to be plotted
features = {'D_{intra} (\mum^2/ms)','N_{branch}', 'L_{branch} (\mum)'};
scales = [ 0.2 0.6; 2 8; 0 100];


% Colour coding by age and region
cerebcor=[0 0 0;0 0.12 0.25;0 0.25 0.5;0 0.4 0.75;0 .5 1.0];
thalamcor=[0 0 0;0.25 0 0;0.5 0 0;0.75 0 0;1.0 0 0];

mse_Cereb = zeros(numel(time_points), numel(metab));
mse_Thalam = zeros(numel(time_points), numel(metab));

for m = 1:numel(metab)
    
    h = figure('Name', ['Cell morphology time evolution for ' metab{m}])
    
    
    for f = 1:numel(features)
        
        feature_mean = zeros(numel(time_points), 1);
        feature_IQR = zeros(numel(time_points), 1);
        cumulative_L = zeros(numel(time_points), 1);
        for time = 1:numel(time_points)
                        
            % Cerebellum
            
            eval(['load ' datafolder '/RF_features_Raw_reduced_ClemenceCerebellum_OtherNandL_' time_points{time} ' morphometry']);
            

            M = table2cell(table(morphometry.Cereb.RFregression_mean));
            feature_mean(time) = M{m}(f);
            
            %M = table2cell(table(morphometry.Cereb.RFregression_25thPercentile));
            %feature_25th = M{m}(f);
            
            %M = table2cell(table(morphometry.Cereb.RFregression_75thPercentile));
            %feature_75th = M{m}(f);
            
            %feature_IQR(time) = (feature_75th - feature_25th)./2;
            
        M = table2cell(table(morphometry.Cereb.RFregression_mean));
%         mean_Nbranch = M{m}(5);
%         mean_L = M{m}(6);
        mean_Nbranch = M{m}(2);
        mean_L = M{m}(3);
        cumulative_L(time) = mean_Nbranch*mean_L;
        end
        
        subplot(2,3,f), hold on
        %errorbar(T, feature_mean, feature_IQR, 'ko-', 'LineWidth', 2.0, 'MarkerSize', 10), title([features{f} '- Cerebellum'])
        plot(T, feature_mean, 'k--', 'LineWidth', 2.0), xlabel('Days from birth'), ylabel(features{f}), axis([0 30 scales(f,:)])
        scatter(T, feature_mean,500,cerebcor,'filled','p')  %         subplot(2,4,numel(features)+1), hold on
%         plot(T, cumulative_L, 'bp--', 'LineWidth', 2.0, 'MarkerSize', 10), xlabel('Days from birth'), ylabel('Average dendritic length [\mum]'), title('Thalamus'), axis([0 30 0 800])
        ax = gca; % current axes
        ax.FontSize = 12;
        ax.TickDir = 'out';
        ax.TickLength = [0.05 0.05];
        ax.XColor ='k';
        ax.YColor ='k';
        ax.FontSmoothing = 'on';
    end
    
    for f = 1:numel(features)
        
        feature_mean = zeros(numel(time_points), 1);
        feature_IQR = zeros(numel(time_points), 1);
        cumulative_L = zeros(numel(time_points), 1);
        for time = 1:numel(time_points)
                        
            % Thalamus
            
            eval(['load ' datafolder '/RF_features_Raw_reduced_ClemenceCerebellum_OtherNandL_' time_points{time} ' morphometry']);
            

            M = table2cell(table(morphometry.Thalam.RFregression_mean));
            feature_mean(time) = M{m}(f);
            
            %M = table2cell(table(morphometry.Thalam.RFregression_25thPercentile));
            %feature_25th = M{m}(f);
            
            %M = table2cell(table(morphometry.Thalam.RFregression_75thPercentile));
            %feature_75th = M{m}(f);
            
            %feature_IQR(time) = (feature_75th - feature_25th)./2;
            
        M = table2cell(table(morphometry.Thalam.RFregression_mean));
%         mean_Nbranch = M{m}(5);
%         mean_L = M{m}(6);
        mean_Nbranch = M{m}(2);
        mean_L = M{m}(3);
        cumulative_L(time) = mean_Nbranch*mean_L;
        
        end
        
        subplot(2,3,f+3), hold on
        %errorbar(T, feature_mean, feature_IQR, 'ro-', 'LineWidth', 2.0, 'MarkerSize', 10), title([features{f} '- Thalamus'])
        plot(T, feature_mean, 'k--', 'LineWidth', 2.0), xlabel('Days from birth'), ylabel(features{f}), axis([0 30 scales(f,:)])
        scatter(T, feature_mean,500,thalamcor,'filled','p')        
%         subplot(2,4,numel(features)+1+4), hold on
%         plot(T, cumulative_L, 'rp-', 'LineWidth', 2.0, 'MarkerSize', 10), xlabel('Days from birth'), ylabel('Average dendritic length [\mum]'), title('Thalamus'), axis([0 30 0 800])
        ax = gca; % current axes
        ax.FontSize = 12;
        ax.TickDir = 'out';
        ax.TickLength = [0.05 0.05];
        ax.XColor ='k';
        ax.YColor ='k';
        ax.FontSmoothing = 'on';

    end
    
    %%% Save figures (1 figure / metabolite)
%     fig = gcf;
%     fig.PaperUnits = 'centimeters';
%     fig.Position = [0 0 20 10];
%     fig.PaperPosition = [0 0 20 11];
    
    %print(['Results-Figures/230920_Morphometric_features_' metab{m} '.png'],'-dpng','-r300')
     
end
