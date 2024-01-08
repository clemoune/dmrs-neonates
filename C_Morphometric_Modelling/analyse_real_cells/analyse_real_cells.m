clear all
close all
clc

%% Read data

mainFolder = fullfile(pwd, 'Cerebellum');
datasets = {'dusart'; 'kengaku'; 'watt'};

listFiles = cell(numel(datasets), 1);

ages = [6, 8, 9, 10, 11, 12, 13, 35, 37, 43];

tree = [];
tree_age = [];
Lstats = [];

c = 1;
for i=1:numel(datasets)
    
    listFiles{i} = dir(fullfile(mainFolder, datasets{i}, 'CNG version', '*.swc'));
    
    for j=1:numel(listFiles{i})
        
        fileName = fullfile(listFiles{i}(j).folder, listFiles{i}(j).name);
        tree{c} = load_tree(fileName);
        for k=1:numel(ages)
            if contains(fileName,['P' num2str(ages(k))])==1
                tree_age(c) = ages(k);
            end
        end
        
        c = c+1;
        
    end

end

tree_stats = [];
for i=1:numel(ages)
            tt = 1;
    for j=1:numel(tree_age)
        if tree_age(j) == ages(i)
            tree_stats{i,tt} = tree(j);
            tt = tt+1;
        end
    end
end

%% 
figure, hold on, 
stats = stats_tree(tree_stats);

%%
branch_len_dist = [];
BO_dist = [];
tot_len_dist = [];

branch_len_qnt = [];
BO_qnt = [];
tot_len_qnt = [];

branch_len_mean = [];
BO_mean = [];
tot_len_mean = [];

branch_len_std = [];
BO_std = [];
tot_len_std = [];

for i=1:numel(ages)
    tmp = stats.dstats(i).blen{1};
    branch_len_dist(i,:) = ksdensity(tmp, 0:1:25, 'support','positive');
    branch_len_qnt(i,:) = quantile(tmp, [0.25 0.50 0.75]);
    branch_len_mean(i,:) = nanmean(tmp(tmp>0));
    branch_len_std(i,:) = nanstd(tmp(tmp>0));

    
    tmp = stats.dstats(i).BO{1};
    tmp(tmp<=0) = 0.1;
    BO_dist(i,:) = ksdensity(tmp, 0:1:25, 'support','positive');
    BO_qnt(i,:) = quantile(tmp, [0.25 0.50 0.75]);
    BO_mean(i,:) = nanmean(tmp(tmp>0));
    BO_std(i,:) = nanstd(tmp(tmp>0));

    
    tmp = stats.dstats(i).Plen{1};
     tmp(tmp<=0) = 0.1;
    tot_len_dist(i,:) = ksdensity(tmp, 0:1:250, 'support','positive');
    tot_len_qnt(i,:) = quantile(tmp, [0.25 0.50 0.75]);
    tot_len_mean(i,:) = nanmean(tmp(tmp>0));
    tot_len_std(i,:) = nanstd(tmp(tmp>0));

end

figure, 
hold on
subplot(1,3,1), errorbar(ages, branch_len_mean, branch_len_std, 'ko--'), ylabel('Branch Length [\mum]'), xlabel('age')
subplot(1,3,2), errorbar(ages, BO_mean, BO_std, 'ko--'), ylabel('Branching Order'), xlabel('age')
subplot(1,3,3), errorbar(ages, tot_len_mean, tot_len_std, 'ko--'), ylabel('Total Process Length [\mum]'), xlabel('age')


%% Thalamus

mainFolder = fullfile(pwd, 'Thalamus');
datasets = {'vonengelhardt'};

listFiles = cell(numel(datasets), 1);

ages = [9, 10, 11, 26, 30];

tree = [];
tree_age = [];
Lstats = [];

c = 1;
for i=1:numel(datasets)
    
    listFiles{i} = dir(fullfile(mainFolder, datasets{i}, 'CNG version', '*.swc'));
    
    for j=1:numel(listFiles{i})
        
        fileName = fullfile(listFiles{i}(j).folder, listFiles{i}(j).name);
        tree{c} = load_tree(fileName);
        for k=1:numel(ages)
            if contains(fileName,['P' num2str(ages(k))])==1
                tree_age(c) = ages(k);
            end
        end
        
        c = c+1;
        
    end

end

tree_stats = [];
for i=1:numel(ages)
            tt = 1;
    for j=1:numel(tree_age)
        if tree_age(j) == ages(i)
            tree_stats{i,tt} = tree(j);
            tt = tt+1;
        end
    end
end

%% 
figure, hold on, 

stats = stats_tree(tree_stats);

%%
branch_len_dist = [];
BO_dist = [];
tot_len_dist = [];

branch_len_qnt = [];
BO_qnt = [];
tot_len_qnt = [];

branch_len_mean = [];
BO_mean = [];
tot_len_mean = [];

branch_len_std = [];
BO_std = [];
tot_len_std = [];

for i=1:numel(ages)
    tmp = stats.dstats(i).blen{1};
    branch_len_dist(i,:) = ksdensity(tmp, 0:1:25, 'support','positive');
    branch_len_qnt(i,:) = quantile(tmp, [0.25 0.50 0.75]);
    branch_len_mean(i,:) = nanmean(tmp(tmp>0));
    branch_len_std(i,:) = nanstd(tmp(tmp>0));

    
    tmp = stats.dstats(i).BO{1};
    tmp(tmp<=0) = 0.1;
    BO_dist(i,:) = ksdensity(tmp, 0:1:25, 'support','positive');
    BO_qnt(i,:) = quantile(tmp, [0.25 0.50 0.75]);
    BO_mean(i,:) = nanmean(tmp(tmp>0));
    BO_std(i,:) = nanstd(tmp(tmp>0));

    
    tmp = stats.dstats(i).Plen{1};
     tmp(tmp<=0) = 0.1;
    tot_len_dist(i,:) = ksdensity(tmp, 0:1:250, 'support','positive');
    tot_len_qnt(i,:) = quantile(tmp, [0.25 0.50 0.75]);
    tot_len_mean(i,:) = nanmean(tmp(tmp>0));
    tot_len_std(i,:) = nanstd(tmp(tmp>0));

end

figure, 
hold on
subplot(1,3,1), errorbar(ages, branch_len_mean, branch_len_std, 'ko--'), ylabel('Branch Length [\mum]'), xlabel('age')
subplot(1,3,2), errorbar(ages, BO_mean, BO_std, 'ko--'), ylabel('Branching Order'), xlabel('age')
subplot(1,3,3), errorbar(ages, tot_len_mean, tot_len_std, 'ko--'), ylabel('Total Process Length [\mum]'), xlabel('age')


