%% Raw Data Visualization for CBPT script
% creates a grand average ERP waveform and superimposes the interquartile
    % range over the grand average ERP wave
% Charts individual waveforms and raster plots from each participant at each electrode;


%% Grand Average ERP Waveforms with Interquartile Range
time = linspace(-200, 800, timePoints); % Change middle variable if your ERP is over a different time period

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
    
    for group = 1:length(groups)
        % Extract data for this group within the channel
        groupCode = group;  % 1 or 2
        groupDat = chanDat(chanDat.Group == groupCode, :);
        
        % Find amplitude columns
        ampCols = contains(groupDat.Properties.VariableNames, 'Time');
        
        % Compute mean across subjects
        meanAmplitudes = mean(groupDat{:, ampCols}, 1);  % 1 x 800
        q25 = prctile(groupDat{:, ampCols}, 25, 1);
        q75 = prctile(groupDat{:, ampCols}, 75, 1);
        subplot(3,1,1)
        hold on;
        fill([time, fliplr(time)], [q25, fliplr(q75)], ...
            groupColors(group, :), ...
            'FaceAlpha', 0.3, ...
            'EdgeColor', 'none', ...
            'DisplayName', char("IQR " + groups(group)));
        
        % Plot
        plot(time, meanAmplitudes, ...
            'DisplayName', groups(group), ...
            'Color', groupColors(group, :), ...
            'LineWidth', 2);
        
        if group == 1
            subplot(3,1,2)
        else
            subplot(3,1,3)
        end
        imagesc(groupDat{:,ampCols}, [-10,40])
    end

    subplot(3,1,1)
    xlabel('Time (ms)');
    ylabel('Mean Amplitude (µV)');
    title(sprintf('Average ERP - Channel %s', chanName));
    ylim([-4 20]); % Can adjust if necessary
    xlim([min(time) max(time)]);
    xline(0, 'k', 'LineWidth', 1.5, 'DisplayName', 'Time 0');
    legend('Location', 'best');
    grid on;
end