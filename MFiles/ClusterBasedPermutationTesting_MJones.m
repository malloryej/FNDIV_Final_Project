%% Cluster Based Permutation for ERP Analysis %%
% written for data collected with predictive T-Maze task using 32 electrode EEG montage
% preprocessed and exported from BrainVision Analyzer v 2.3.0
% Two primary components of this script:
% 1) Permutation testing - randomly redistributes or "scrambles" data and reassigns to 
    % groups for nPermutations. Each permutation will result in a structured array of clusters. 
    % The largest cluster (i.e. the one with the most consecutively significant [at p <= 0.01] milliseconds) 
    % from each permutation will be added to the null distribution.
% 2) Parametric testing - T-test conducted at every millisecond of the ERP. 
    % After testing, clusters are computed. Clusters of size >= the 95th percentile of the null distribution 
    % will be considered significant at the 0.05 alpha level.

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
groups = data.Group; % creates a variable to indicate whether participant was/was not a ruminator; edit to match your category label
channels = data.Channel; % creates a variable to hold channel information

varNames = data.Properties.VariableNames;
timepointCols = varNames(startsWith(varNames, 'Time')); % creates an array that holds all of the amplitudes at every channel at every timepoint

uniqueChannels = unique(channels); % makes a list of unique channels being examined; should match input selected in preproc script

% initialize empty array for results
results = perform_tTests(data, uniqueChannels, 'Group')


%% Determine clusters %%
% Define significance threshold
alpha = 0.05; % Can be adapted to be more strict
min_cluster = 1; % Minimum number of consecutive timepoints; calculated with DetermineClusterSize script

% Initialize empty cell for clusters and empty matrix for pval < 0.05
clustersAll = struct('Channel', {}, 'TStat', {}, 'PVal', {}, 'Timepoint', {});

sig_mask = zeros(nChannels, timePoints); % create significance mask

% loop through significance values in each channel
for c = 1:nChannels
    for t = 1:timePoints;
        row = (c - 1) * timePoints + t;
        sig_mask(c, t) = results(row).PVal < 0.01;
    end
end

clusters = find_clusters(sig_mask);
for cluster = 1:length(clusters)
    results(row) = struct('Channel', chan, ...
        'TStat', stats.tstat, ...
        'PVal', p, ...
        'Timepoint', t);
    row = row + 1;
end