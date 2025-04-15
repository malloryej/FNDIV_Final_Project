%% Preprocessing script for Cluster Based Permutation Testing script
% Designed for data pre-processed in BrainVision Analyzer version 2.3
    % preprocessing includes interpolation of channels *(where necessary), 
    % ocular correction, artifact rejection, segmentation, and averaging of
    % ERPs across trials
% takes in a folder with a csv file with ERP data for each participant
% reconstructs data from all individual suject files into one table
% formats one table to be compatable with ClusterBasedPermutationTesting
% script
% REQUIREMENTS
    % SubID "key" that decodes which subject ID is associated with either
        % condition (i.e. ruminator vs nonruminator
    % folder that contains files directly exported from BrainVision analyzer
        % should have one file for each participant that contains data for
        % ERP averaged over trials at every millisecond for each channel
        % recorded

% Mallory Jones 2025

%% Set folder path and initialize variables
folderPath = "C:\Users\mallo\OneDrive\Desktop\ClusterBasedPermutationTesting\Data\InputData\datFiles\datFiles_StimorInhib";
files = dir(fullfile(folderPath, '*.csv'));
SubID_condKey = readtable("C:\Users\mallo\OneDrive\Desktop\ClusterBasedPermutationTesting\Data\InputData\PreprocData\SubList_StimType.xlsx");

%% Create a list of subIDs
SubList = table();
for file = 1:length(files);
    fileName = files(file).name;
    SubID_Sub = erase(fileName, '.csv');
    Dat = string(SubID_Sub);
    tabDat = table(Dat);
    SubList = vertcat(SubList, tabDat);
end

SubList.Properties.VariableNames = "SubID";

%% Reconstruct single table with data from all timepoints for each participant
allData = table();
for sub = 1:height(SubList); % loops through subID list made in previous step
    subData = readtable(fullfile(folderPath, SubList.SubID(sub)));
    numRows = height(subData);
    subIDColumn = repmat(SubList.SubID(sub), numRows, 1); %repeats subject ID for number of rows associated with that participant's data
    subData.SubID = subIDColumn; % adds subID column to the new table
    subGroupColumn = repelem(SubID_condKey.StimType(sub), numRows, 1); % uses subID "key" to assign participants to a group 
    subData.Group = subGroupColumn; % adds group membership identifier
    for i = 1:width(subData)
        colData = subData{:,i};
        if any(strcmp(colData, 'Cz'));
            subData.Properties.VariableNames{i} = 'Channel';
            break
        end
    end
    %subData = subData(:, [end, 1:end-1]);
    %subData = subData(:, [end, 1:end-1]);
    allData = vertcat(allData, subData);
end

%% Update Variable Order and Names
colsToMove = {'SubID', 'Group'};  % columns you want at the beginning
allVars = allData.Properties.VariableNames;
newOrder = [colsToMove, setdiff(allVars, colsToMove, 'stable')];
allData = allData(:, newOrder);

% Define fixed columns
fixedCols = {'SubID', 'Group', 'Channel'};

% Identify timepoint columns (everything else)
timeCols = setdiff(allData.Properties.VariableNames, fixedCols, 'stable');

% Build full list of new variable names
varNames = [fixedCols, strcat("Time", string(1:length(timeCols)))];

% Assign new variable names
allData.Properties.VariableNames = varNames;

%% write table as csv

allData = allData(ismember(allData.Channel, ["FCz", "Fz", "Cz"]), :); %editable! Can add any channels by adding the string name for the channel

writetable(allData, "C:\Users\mallo\OneDrive\Desktop\ClusterBasedPermutationTesting\Data\InputData\" + ...
    "PreprocData\TrialData.csv") % Saves data as csv file; can change file name here