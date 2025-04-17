%% Final Data Visualization for CBPT script
% creates a grand average ERP waveform and superimposes the standard error
    % range over the grand average ERP wave
% shades regions where clusters are significant

%% Grand Average ERP Waveforms with standard error

% define time windows for shading
clustersShade = struct();  % initialize empty struct

for i = 1:length(clustersSig)
    ch = clustersSig(i).Channel;
    startT = clustersSig(i).Timepoint_start;
    endT = clustersSig(i).Timepoint_end;

    if isfield(clustersShade, ch)
        clustersShade.(ch) = [clustersShade.(ch); startT endT];  % append
    else
        clustersShade.(ch) = [startT endT];  % first entry
    end
end

% Align with timepoint 0 (i.e. -200 to 600 timepoints instead of 0 to 800)
for i = 1:length(fields)
    ch = fields{i};
    % Convert sample indices to actual time values
    idx = clustersShade.(ch);  % e.g., [123 150]
    clustersShade.(ch) = time(idx);  % replaces with [time(123) time(150)]
end

% generate ERP wave forms with SE and cluster shading

time = linspace(-200, 600, timePoints);

groupColors = [
    0.2 0.4 1.0;   % Blue-ish for Group 1
    1.0 0.3 0.3    % Reddish for Group 2
];
for chan = 1:length(channels)
    % Get channel name
    chanName = channels(chan);
    
    % Extract data for that channel
    chanDat = allData(strcmp(allData.Channel, chanName), :);
    
    % Open figure ONCE per channel
    figure('WindowStyle','docked');
    hold on;

    for group = 1:length(groups)
        % Extract data for this group within the channel
        groupCode = group;  % 1 or 2
        groupDat = chanDat(chanDat.Group == groupCode, :);
        
        % Find amplitude columns
        ampCols = contains(groupDat.Properties.VariableNames, 'Time');
        
        % Compute mean across subjects
        meanAmplitudes = mean(groupDat{:, ampCols}, 1);  % 1 x 800
        se = std(groupDat{:, ampCols}, 0, 1) / sqrt(size(groupDat{:, ampCols}, 1)); 
        fill([time, fliplr(time)], ...
            [meanAmplitudes - se, ...
            fliplr(meanAmplitudes + se)], ...
            groupColors(group, :), ...
            'FaceAlpha', 0.3, ...
            'EdgeColor', 'none', ...
            'DisplayName', char("SE " + groups(group)));
        
        % Plot
        plot(time, meanAmplitudes, ...
            'DisplayName', groups(group), ...
            'Color', groupColors(group, :), ...
            'LineWidth', 2);
    end

    % Add shaded cluster regions if defined for this channel
    if isfield(clustersShade, char(chanName))
        thisCluster = clustersShade.(char(chanName));

        yMin = -6;
        yMax = 6;
    
        for i = 1:size(thisCluster, 1)
            x1 = thisCluster(i, 1);
            x2 = thisCluster(i, 2);
            fill([x1 x2 x2 x1], [yMin, yMin, yMax, yMax],...
                [0.7 0.7 0.7], 'FaceAlpha', 0.4, 'EdgeColor', 'none', ...
                'HandleVisibility', 'off');
        end
    end
    xlabel('Time (ms)');
    ylabel('Mean Amplitude (µV)');
    title(sprintf('Average ERP - Channel %s', chanName));
    ylim([-6 6]);
    xlim([min(time) max(time)]);
    xline(0, 'k', 'LineWidth', 1.5, 'DisplayName', 'Time 0');
    legend('Location', 'best');
    grid on;
end