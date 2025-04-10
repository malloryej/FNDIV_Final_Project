%% Preprocessing script for Cluster Based Permutation Testing script
% Designed for data preprocessed with CBPT_Preproc_MJones_V4 script
% takes in the file created by preproc script
% Conducts nPermutation number of t-tests to generate a null distribution
% REQUIREMENTS
    % TrialData file generated by preproc script

% Mallory Jones 2025

%% Set Parameters %%

dir = "C:\Users\mallo\OneDrive\Desktop\ClusterBasedPermutationTesting\Preproc";
dirOut = "C:\Users\mallo\OneDrive\Desktop\ClusterBasedPermutationTesting\Data\InputData\";
data = readtable ("C:\Users\mallo\OneDrive\Desktop\ClusterBasedPermutationTesting\Data\InputData\" + ...
    "PreprocData\TrialData");
nPermutations = 20; %number of permutations; can be adjusted; 1000 generally accepted minimum
nChannels = 3; %number of channels; can be adjusted; depends on which channels were identified in preproc document
timePoints = 800; %total number of timepoints being investigated; can be adjusted based on length of ERP extracted from BrainVision
channels = data.Channel;
uniqueChannels = unique(channels);

%% Create function to Identify Clusters in Data
function listlengths = find_clusters(significant_mask);
    listlengths = [];
    [num_rows, num_cols] = size(significant_mask);

    for row = 1:num_rows
        cluster_length = 0;

        for col = 1:num_cols
            if significant_mask(col) == 1
                cluster_length = cluster_length + 1;

            elseif cluster_length > 0 % only add to listlengths if cluster length > 0
                % ends current cluster; resets variables
                listlengths = [listlengths, cluster_length];
                cluster_length = 0;

            end
        end

        if cluster_length > 0
            listlengths = [listlengths, cluster_length];
            % handle cases where last timepoint is still in cluster

        end
    end
end

        
%% Create Null Distribution %%
% repeats previous steps for randomized/nonsense data

null_cluster_size = [];
for perm = 1:nPermutations;
    data.Ruminator = randi([1,2], height(data), 1); %randomizes group assignment for each participant and electrode

% initialize empty array for results
    groups = data.Ruminator;
    channels = data.Channel;
    timepointData = table2array(data(:, 4:end));
    
    resultsPerm = cell(nChannels, timePoints);
    
    for c = 1:nChannels
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

    sig_maskPerm = zeros(nChannels, timePoints); % create significance mask
    
    % loop through significance values in each channel
    for c = 1:length(uniqueChannels)
        for t = 1:timePoints;
            sig_maskPerm(c, t) = resultsPerm{c, t}.pval < 0.01; % extracts pval from results cell
        end
    end
    
    % find clusters - calls function we defined above
    clusters = find_clusters(sig_maskPerm);

    if isempty(clusters)
        max_cluster_size = 0;
    else
        max_cluster_size = max(clusters);
    end
    null_cluster_size = [null_cluster_size, max_cluster_size];
end

  

%% Extract value for the 95th percentile of data in the null distribution %%

min_cluster_size = prctile(null_cluster_size, 95)