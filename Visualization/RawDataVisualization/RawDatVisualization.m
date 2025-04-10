%% Raw Data Visualization for CBPT script
% creates a grand average ERP waveform and superimposes the interquartile
    % range over the grand average ERP wave
% Charts individual waveforms from each participant at each electrode;
    % helps to identify potential outliers
% Creates histograms for amplitudes at assumed key points in the data (i.e. hypothesized
    % P3 and N2 peak as well as a baseline timepoint)
    % also helpful for identifying outliers
% REQUIREMENTS TrialData file generated by original preproc script

% Mallory Jones 2025


%% Set Parameters %%
dir = "C:\Users\mallo\OneDrive\Desktop\ClusterBasedPermutationTesting\Preproc";
dirOut = "C:\Users\mallo\OneDrive\Desktop\ClusterBasedPermutationTesting\Data\InputData\";
data = readtable ("C:\Users\mallo\OneDrive\Desktop\ClusterBasedPermutationTesting\Data\InputData\" + ...
    "PreprocData\TrialData");