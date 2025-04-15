%% Create function to perform tTests on null data

function resultsStruct = perform_tTests(data, chanList, interestVar)

    row = 1; % count where data will be added to structured result array

    % separate timepoint amplitudes from subID, channel info, and condition
        % helps with iterating across timepoints later
    varNames = data.Properties.VariableNames;
    timepointCols = varNames(startsWith(varNames, 'Time'));

    % initialize empty structured array for results
    resultsStruct = struct('Channel', {}, 'TStat', {}, 'PVal', {}, 'Timepoint', {});
    
    %perform t-tests
    for c = 1:length(chanList)
        chan = chanList(c);
        chanData = data(strcmp(data.Channel, chan), :);

        % Separate Groups
        group1 = chanData(chanData.(interestVar) == 1, :);
        group2 = chanData(chanData.(interestVar) == 2, :);
        group1Data = table2array(group1(:, timepointCols));
        group2Data = table2array(group2(:, timepointCols));

        for t = 1:length(timepointCols)
            [~, p, ~, stats] = ttest2(group1Data(:, t), group2Data(:, t));
            resultsStruct(row) = struct('Channel', chan, ...
                'TStat', stats.tstat, ...
                'PVal', p, ...
                'Timepoint', t);
            row = row + 1;
        end
    end
end