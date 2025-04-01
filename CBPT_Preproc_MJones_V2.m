%% Preprocessing script for Cluster Based Permutation Testing script

% Set folder path and initialize variables
folderPath = "C:\Users\mallo\OneDrive\Desktop\Cluster Based Permutation Testing\.csvFiles\";
files = dir(fullfile(folderPath, '*.csv'));
SubID_StimKey = readtable("C:\Users\mallo\OneDrive\Desktop\Cluster Based Permutation Testing\SubList_StimType.xlsx")


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

%% Create table with data from all timepoints for each participant
allData = table();
for sub = 1:height(SubList);
    subData = readtable(fullfile(folderPath, SubList.SubID(sub)));
    numRows = height(subData);
    subIDColumn = repmat(SubList.SubID(sub), numRows, 1);
    subData.subID = subIDColumn;
    subStimColumn = repelem(SubID_StimKey.StimType(sub), numRows, 1);
    subData.stimType = subStimColumn;
    subData = subData(:, [end, 1:end-1]);
    subData = subData(:, [end, 1:end-2]);
    allData = vertcat(allData, subData);
end


varNames = ["SubID", "Ruminator", "Channel"]; % Initiate list for timepoint variable names; starts by naming the first column "Channels"
count = 0;
for col = 1:(width(allData) - 3);
    count = count + 1;
    varNames = [varNames, strcat("T" + count)];
end

allData.Properties.VariableNames = varNames; % assigns variable names


writetable(allData, string(folderPath + 'TrialData.csv'))