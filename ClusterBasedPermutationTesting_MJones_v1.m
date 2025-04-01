%% Cluster Based Permutation for ERP Analysis %%
% written for data collected with predictive T-Maze task using 32 electrode EEG montage
% preprocessed and exported from BrainVision Analyzer v 2.3.0
% Two primary components of this script:
    % 1) Parametric testing and clustering (Quantifying effect at every timepoint) - T-test/ANOVA at every millisecond of the ERP
        % T-test is used in this case for rumination data (i.e. either high
        % or low)
        % ANOVA is used for learning style data (i.e. random, choice kernel, win
        % stay lose shift, or Riscorla-Wagner)
        % After testing, identify all effects greater than some parameter
        % (i.e. alpha < 0.05 or alpha < 0.01)
        % sum all t-statistics or F-values at each electrode across all time points
            % This script specifies a cluster as an effect size with
            % p-value < 0.05 over the course of 10 or more consecutive
            % milliseconds between 3 or more adjacent electrodes. This
            % definition of cluster can be changed based on the user's
            % preferences
    % 2) Permutation testing
        % Randomly redistribute data in groups consistent with conditions
        % (i.e. for low and high ruminators, there will be two groups and
        % for learning styles there will be four groups)
        % repeat parametric testing and clustering steps nPermutations times to
        % create a null distribution; you will choose the largest effect
        % size from each permutation to add to the null distribution (i.e.
        % null distribution will consist of nPermutations total measures)
        % clusters with sum of the test statistics >= than the top 5
        % percentile will be considered significant

% Mallory Jones 2025

%% Set up %%

dir = "C:\Users\mallo\OneDrive\Desktop\Cluster Based Permutation Testing\";
dirOut = "C:\Users\mallo\OneDrive\Desktop\Cluster Based Permutation Testing\";
data = readtable ("C:\Users\mallo\OneDrive\Desktop\Cluster Based Permutation Testing\TrialData");
nPermutations = 5; %number of permutations; can be adjusted
nChannels = 23; %number of channels; can be adjusted
timePoints = 800; %total number of timepoints being investigated; can be adjusted

%% T-tests %%
% Runs a t-test comparing amplitude of EEG from ruminators versus
% non-ruminators at every millisecond of the ERP
% Comment this section out to run ANOVAS (next section)

% separate ruminators from nonruminators
ruminators = data(data.Ruminator == 2, :);
nonruminators = data(data.Ruminator == 1, :);

% initialize empty array for results
results = cell(nChannels, timePoints);

for t = 1:timePoints;
    for c = 1:nChannels;

        % Extract EEG data at current timepoint and channel
        ruminator_data = ruminators{ruminators.Time_inMs_ == t, "Channel" + num2str(c)};
        nonruminator_data = nonruminators{nonruminators.Time_inMs_ == t, "Channel" + allData.Channel(c)};

        % Perform t-test
        [h, p, ci, stats] = ttest2(ruminator_data, nonruminator_data);
        results{c, t} = struct('t_stat', stats.tstat, 'pval', p);
    end
end

%% Determine clusters %%
% Define significance threshold
alpha = 0.05; % Can be adapted to be more strict
min_cluster = 5; % Minimum number of consecutive timepoints that must be lower than alpha to form cluster

% Initialize empty cell for clusters and empty matrix for pval < 0.05
% logicals
clusters = {};
cluster_idx = 1;
sig_mask = zeros(nChannels, timePoints);

% loop through significance values in each channel
for c = 1:nChannels;
    for t = 1:timePoints;
        sig_mask = sig_mask + results{c, t}.pval > 0.05; %extracts pval from results cell
    end
end

% counts number of consecutive sig pvalues in each channel
for c = 1:nChannels;
    cons_sigs = 0; %lets you index consecutive significant p-vals
    cluster_start = NaN; % initialize cluster start time

    % loops through time point in each channel; if there is a significant
    % time point, cons_sigs tracker increases by 1. if cons_sigs was
    % previously set to zero, a new cluster begins
    for t = 1:size(sig_mask, 2);
        if sig_mask(c, t);
            if cons_sigs == 0;
                cluster_start = t;
            end
            cons_sigs = cons_sigs + 1;

            % if not significant, check to see if cluster is large enough
            % to be considered significant; otherwise, starts back at the
            % next iteration of the timepoint for loop
        else
            if cons_sigs >= min_cluster;
                clusters(cluster_idx) = struct('Channel', c, 'StartTime', cluster_start, 'EndTime', t-1)
                cluster_idx = cluster_idx + 1
            end
            cons_sigs = 0;
            cluster_start = NaN;
        end
    end

    % if cluster continues until end, include that cluster in the cell
    if cons_sigs >= min_cluster;
        clusters{cluster_idx} = struct('Channel', c, 'StartTime', cluster_start, 'EndTime', size(sig_mask,2));
        cluster_idx = cluster_idx + 1;
    end
end

% sum t-statistics for each cluster %
for i = 1:length(clusters)
    ch = clusters{i}.Channel;
    start_time = clusters{i}.StartTime;
    end_time = clusters{i}.EndTime;
    
    % Extract t-statistics from results within the cluster time range
    t_values_abs = abs(arrayfun(@(t) results{ch, t}.t_stat, start_time:end_time));
    
    % Sum of t-statistics for the cluster
    clusters{i}.Tsum = sum(t_values_abs);
end

%% Create Null Distribution %%

dataPerms = data
null_distribution = []
for p = 1:nPermutations;
    dataPerms.Ruminator = randi([1,2], height(data), 1);
    % separate ruminators from nonruminators
    ruminatorsPerm = dataPerms(dataPerms.Ruminator == 2, :);
    nonruminatorsPerm = dataPerms(dataPerms.Ruminator == 1, :);
% initialize empty array for results
    resultsPerm = cell(nChannels, timePoints);
    for c = 1:nChannels;
        for t = 1:timePoints;
    
            % Extract EEG data at current timepoint and channel
            ruminator_dataPerm = ruminatorsPerm{ruminatorsPerm.Time_inMs_ == t, "Channel" + num2str(c)};
            nonruminator_dataPerm = nonruminatorsPerm{nonruminatorsPerm.Time_inMs_ == t, "Channel" + num2str(c)};
    
            % Perform t-test
            [h, p, ci, stats] = ttest2(ruminator_dataPerm, nonruminator_dataPerm);
            resultsPerm{c, t} = struct('t_stat', stats.tstat, 'pval', p);
        end
    end
    clustersPerm = {};
    cluster_idx = 1;
    sig_mask = zeros(nChannels, timePoints);
    % loop through significance values in each channel
    for c = 1:nChannels;
        for t = 1:timePoints;
            sig_mask = sig_mask + resultsPerm{c, t}.pval > 0.05; %extracts pval from results cell
        end
    end
    
    % counts number of consecutive sig pvalues in each channel
    for c = 1:nChannels;
        cons_sigs = 0; %lets you index consecutive significant p-vals
        cluster_start = NaN; % initialize cluster start time
    
        % loops through time point in each channel; if there is a significant
        % time point, cons_sigs tracker increases by 1. if cons_sigs was
        % previously set to zero, a new cluster begins
        for t = 1:size(sig_mask, 2);
            if sig_mask(c, t);
                if cons_sigs == 0;
                    cluster_start = t;
                end
                cons_sigs = cons_sigs + 1;
    
                % if not significant, check to see if cluster is large enough
                % to be considered significant; otherwise, starts back at the
                % next iteration of the timepoint for loop
            else
                if cons_sigs >= min_cluster;
                    clustersPerm(cluster_idx) = struct('Channel', c, 'StartTime', cluster_start, 'EndTime', t-1)
                    cluster_idx = cluster_idx + 1
                end
                cons_sigs = 0;
                cluster_start = NaN;
            end
        end
    
        % if cluster continues until end, include that cluster in the cell
        if cons_sigs >= min_cluster;
            clustersPerm{cluster_idx} = struct('Channel', c, 'StartTime', cluster_start, 'EndTime', size(sig_mask,2));
            cluster_idx = cluster_idx + 1;
        end
    end
    TList = [];
    for i = 1:length(clustersPerm)
        ch = clustersPerm{i}.Channel;
        start_time = clustersPerm{i}.StartTime;
        end_time = clustersPerm{i}.EndTime;
        
        % Extract t-statistics from results within the cluster time range
        t_values_abs = abs(arrayfun(@(t) resultsPerm{ch, t}.t_stat, start_time:end_time));
        
        % Sum of t-statistics for the cluster
        clustersPerm{i}.Tsum = sum(t_values_abs);
        TList = [TList, clustersPerm{i}.Tsum];
        disp(TList)
    end
    disp(max(TList))
    null_distribution = [null_distribution, max(TList)];
end

%% Extract value for the 95th percentile of data in the null distribution %%

null_percentile_95th = prctile(null_distribution, 95)

%% Compare magnitude of cluster effect size (i.e. sum of t-vals for each significant timepoint at an electrode) to the 95th percentile value of the null_distribution
sig_clusters = [];
for i = 1:length(clusters);
    if clusters{i}.Tsum <= null_percentile_95th;
        sig_clusters{i}.start = clusters{i}.StartTime
        sig_clusters{i}.end = clusters{i}.EndTime
        sig_clusters{i}.channel = clusters{i}.Channel
    end
end
