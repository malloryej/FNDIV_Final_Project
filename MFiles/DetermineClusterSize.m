%% Preprocessing script for Cluster Based Permutation Testing script
% Conducts nPermutation number of t-tests to generate a null distribution
% 

% Mallory Jones 2025

%% Set Parameters %%

dir = "C:\Users\mallo\OneDrive\Desktop\ClusterBasedPermutationTesting\Preproc";
dirOut = "C:\Users\mallo\OneDrive\Desktop\ClusterBasedPermutationTesting\Data\InputData\";
allData = readtable ("C:\Users\mallo\OneDrive\Desktop\ClusterBasedPermutationTesting\Data\InputData\" + ...
    "PreprocData\TrialData");
nPermutations = 20; %number of permutations; can be adjusted; 1000 generally accepted minimum
nChannels = 3; %number of channels; can be adjusted; depends on which channels were identified in preproc document
timePoints = 800; %total number of timepoints being investigated; can be adjusted based on length of ERP extracted from BrainVision
channels = allData.Channel;
uniqueChannels = string(unique(channels));
varofInterest = 'Ruminator';

        
%% Create Null Distribution %%

null_dist = [];
for perm = 1:nPermutations
    dataPerm = allData;
    dataPerm.Group = randi([1,2], height(dataPerm), 1);
    resultsStructPerm = perform_tTests(dataPerm, channels);

    % find null, permuted clusters
    clustersPerm = find_clusters(resultsStructPerm, null_alpha);

    if isempty(clustersPerm) % if no cluster, add 0 as the max_cluster_size for that permutation
        max_cluster_size = 0;
    else
        max_cluster_size = max([clustersPerm.Cluster_length]);
    end
    null_dist = [null_dist, max_cluster_size];
end

  
%% Extract value for the 95th percentile of data in the null distribution %%
% min cluster is the minimum number of consecutive timepoints that must be
% significant for cluster to be significant

if max(null_dist) >= 1
    min_cluster_size = prctile(null_dist, 95);
else
    min_cluster_size = 1; % if all permutations yielded 0 for max cluster, set threshold to 1
end