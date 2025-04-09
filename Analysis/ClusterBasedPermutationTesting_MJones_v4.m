%% Cluster Based Permutation for ERP Analysis %%
% written for data collected with predictive T-Maze task using 32 electrode EEG montage
% preprocessed and exported from BrainVision Analyzer v 2.3.0
% Two primary components of this script:
    % 1) Parametric testing and clustering (Quantifying effect at every timepoint) - T-test at every millisecond of the ERP
        % T-test is used in this case for rumination data (i.e. either high
        % or low)
        % After testing, identify all effects greater than some parameter
        % (i.e. alpha < 0.05 or alpha < 0.01)
        % AVERAGE OR SUM all t-statistics or F-values at each electrode across all time points
            % This script specifies a cluster as having a significance values
            % (i.e. p-value) < 0.05 over the course of 10 or more consecutive
            % milliseconds.
                % These cluster parameters can be changed based on the user's preferences
    % 2) Permutation testing
        % Randomly redistribute or "scramble" data and reassign to two new
        % groups
        % repeat parametric testing and clustering steps nPermutations times to
        % create a null distribution; you will choose the largest effect
        % size from each permutation to add to the null distribution (i.e.
        % null distribution will consist of nPermutations total measures)
        % clusters with SUM OR AVERAGE of test statistics >= than the 95
        % percentile will be considered significant

% Mallory Jones 2025

%% Set up %%

dir = "C:\Users\mallo\OneDrive\Desktop\ClusterBasedPermutationTesting\";
dirOut = "C:\Users\mallo\OneDrive\Desktop\ClusterBasedPermutationTesting\";
data = readtable ("C:\Users\mallo\OneDrive\Desktop\ClusterBasedPermutationTesting\TrialData");
nPermutations = 20; %number of permutations; can be adjusted; 1000 generally accepted minimum
nChannels = 3; %number of channels; can be adjusted; depends on which channels were identified in preproc document
timePoints = 800; %total number of timepoints being investigated; can be adjusted based on length of ERP extracted from BrainVision

%% T-tests %%
% Runs a t-test comparing amplitude of EEG from ruminators versus
% non-ruminators at every millisecond of the ERP

% Set up key variables
groups = data.Ruminator; % creates a variable to indicate whether participant was/was not a ruminator
channels = data.Channel; % creates a variable to hold channel information
timepointData = table2array(data(:, 4:end)); % creates an array that holds all of the amplitudes at every channel at every timepoint
uniqueChannels = unique(channels); % makes a list of unique channels being examined; should match input selected in preproc script

% initialize empty array for results
results = cell(nChannels, timePoints);

for c = 1:length(uniqueChannels) % sets channel for each round of t-tests; each "round" of t-tests is conducted for one channel at a time 
    channelData = timepointData(strcmp(channels, uniqueChannels{c}), :); % extracts data only from channel currently being evaluated
    channelGroups = groups(strcmp(channels, uniqueChannels{c})); % creates channel specific variable to separate ruminator from nonruminators

    % separates ruminator channel data from nonruminator channel data
    ruminators = channelData(channelGroups == 1, :);
    nonruminators = channelData(channelGroups == 2, :);

    for t = 1:timePoints % performs t-test at every timepoint for each channel

    % Perform t-test
        [h, p, ci, stats] = ttest2(ruminators(:, t), nonruminators(:, t));
        results{c, t} = struct('t_stat', stats.tstat, 'pval', p);
    end
end

%% Determine clusters %%
% Define significance threshold
alpha = 0.05; % Can be adapted to be more strict
min_cluster = 1; % Minimum number of consecutive timepoints that must be lower than alpha to form cluster; can be adjusted

% Initialize empty cell for clusters and empty matrix for pval < 0.05
% logicals (i.e. binary true or false)
clusters = {};
cluster_idx = 1;
sig_mask = zeros(nChannels, timePoints);

% loop through significance values at all timepoints for each channel
for c = 1:nChannels;
    for t = 1:timePoints;
        sig_mask(c, t) = results{c, t}.pval < 0.05; % extracts pval from results cell
    end
end

% counts number of consecutive sig pvalues in each channel
for c = 1:nChannels;
    cons_sigs = 0; % lets you count how many consecutive significant p-vals
    cluster_start = NaN; % initialize cluster start time

    % loops through time point in each channel; if there is a significant
    % time point, cons_sigs tracker increases by 1.
    for t = 1:size(sig_mask, 2);
        if sig_mask(c, t); % evaluates whether p-value is significant or not; uses logical variable
            if cons_sigs == 0; % if cons_sigs was previously set to zero, a new cluster begins
                cluster_start = t;
            end
            cons_sigs = cons_sigs + 1;

            % if not significant, check to see if cluster is large enough
            % to be considered significant; otherwise, starts back at the
            % next iteration of the timepoint for loop
        else
            if cons_sigs >= min_cluster;
                clusters{cluster_idx} = struct('Channel', uniqueChannels{c}, 'StartTime', cluster_start, 'EndTime', t-1);
                cluster_idx = cluster_idx + 1;
            end
            cons_sigs = 0;
            cluster_start = NaN;
        end
    end

    % if cluster continues until end, include that cluster in the cell
    if cons_sigs >= min_cluster;
        clusters{cluster_idx} = struct('Channel', uniqueChannels{c}, 'StartTime', cluster_start, 'EndTime', size(sig_mask,2));
        cluster_idx = cluster_idx + 1;
    end
end

% averages t-statistics for each cluster %
for i = 1:length(clusters)
    ch_name = clusters{i}.Channel;
    start_time = clusters{i}.StartTime;
    end_time = clusters{i}.EndTime;

    ch_idx = find(strcmp(uniqueChannels, ch_name)); % links channel number back to the channel label
    
    % Extract t-statistics from results within the cluster time range
    t_values_abs = abs(arrayfun(@(t) results{ch_idx, t}.t_stat, start_time:end_time));
    
    % Sum of t-statistics for the cluster
    clusters{i}.Tavg = mean(t_values_abs);
end

%% Create Null Distribution %%
% repeats previous steps for randomized/nonsense data

dataPerms = data;
null_distribution = [];
for perm = 1:nPermutations;
    dataPerms.Ruminator = randi([1,2], height(data), 1); %randomizes group assignment for each participant

% initialize empty array for results
    groups = dataPerms.Ruminator;
    channels = dataPerms.Channel;
    timepointData = table2array(dataPerms(:, 4:end));
    
    resultsPerm = cell(nChannels, timePoints);
    
    for c = 1:length(uniqueChannels)
        channelData = timepointData(strcmp(channels, uniqueChannels{c}), :);
        channelGroups = groups(strcmp(channels, uniqueChannels{c}));
        ruminatorsPerm = channelData(channelGroups == 1, :);
        nonruminatorsPerm = channelData(channelGroups == 2, :);
    
        for t = 1:timePoints
    
        % Perform t-test
            [h, p, ci, stats] = ttest2(ruminatorsPerm(:, t), nonruminatorsPerm(:, t));
            resultsPerm{c, t} = struct('t_stat', stats.tstat, 'pval', p);
        end
    end

    clustersPerm = {};
    cluster_idxPerm = 1;
    sig_maskPerm = zeros(nChannels, timePoints);
    
    % loop through significance values in each channel
    for c = 1:nChannels;
        for t = 1:timePoints;
            sig_maskPerm(c, t) = resultsPerm{c, t}.pval < 0.05; % extracts pval from results cell
        end
    end
    
    % counts number of consecutive sig pvalues in each channel
    for c = 1:nChannels;
        cons_sigs = 0; %lets you index consecutive significant p-vals
        cluster_start = NaN; % initialize cluster start time
    
        % loops through time point in each channel; if there is a significant
        % time point, cons_sigs tracker increases by 1. if cons_sigs was
        % previously set to zero, a new cluster begins
        for t = 1:size(sig_maskPerm, 2);
            if sig_maskPerm(c, t);
                if cons_sigs == 0;
                    cluster_start = t;
                end
                cons_sigs = cons_sigs + 1;
    
                % if not significant, check to see if cluster is large enough
                % to be considered significant; otherwise, starts back at the
                % next iteration of the timepoint for loop
            else
                if cons_sigs >= min_cluster;
                    clustersPerm{cluster_idxPerm} = struct('Channel', uniqueChannels{c}, 'StartTime', cluster_start, 'EndTime', t-1);
                    cluster_idxPerm = cluster_idxPerm + 1;
                end
                cons_sigs = 0;
                cluster_start = NaN;
            end
        end
    
        % if cluster continues until end, include that cluster in the cell
        if cons_sigs >= min_cluster;
            clustersPerm{cluster_idxPerm} = struct('Channel', uniqueChannels{c}, 'StartTime', cluster_start, 'EndTime', size(sig_mask,2));
            cluster_idx = cluster_idxPerm + 1;
        end
    end
    
    % average t-statistics for each cluster %
    TList = [];
    for i = 1:length(clustersPerm)
        ch_name = clustersPerm{i}.Channel;
        start_time = clustersPerm{i}.StartTime;
        end_time = clustersPerm{i}.EndTime;
    
        ch_idxPerm = find(strcmp(uniqueChannels, ch_name)); % links channel number back to the channel label
        
        % Extract t-statistics from results within the cluster time range
        t_values_abs = abs(arrayfun(@(t) results{ch_idxPerm, t}.t_stat, start_time:end_time));
        
        % Sum of t-statistics for the cluster
        clustersPerm{i}.Tavg = mean(t_values_abs);
        TList = [TList, clustersPerm{i}.Tavg];
    end
    null_distribution = [null_distribution, max(TList)];
end

%% Extract value for the 95th percentile of data in the null distribution %%

null_percentile_95th = prctile(null_distribution, 95)

%% Compare magnitude of cluster effect size (i.e. sum of t-vals for each significant timepoint at an electrode) to the 95th percentile value of the null_distribution
sig_clusters = {};
for i = 1:length(clusters);
    if clusters{i}.Tavg >= null_percentile_95th;
        sig_clusters{end+1}.start = clusters{i}.StartTime;  % Append to end of cell array
        sig_clusters{end}.end = clusters{i}.EndTime;
        sig_clusters{end}.channel = clusters{i}.Channel;
    end
end
