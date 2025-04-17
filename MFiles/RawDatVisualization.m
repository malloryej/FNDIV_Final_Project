%% Raw Data Visualization for CBPT script
% creates a grand average ERP waveform and superimposes the interquartile
    % range over the grand average ERP wave
% Charts individual waveforms from each participant at each electrode;
    % helps to identify potential outliers
% Creates histograms for amplitudes at assumed key points in the data (i.e. hypothesized
    % P3 and N2 peak as well as a baseline timepoint)
    % also helpful for identifying outliers


%% Grand Average ERP Waveforms with Interquartile Range
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
    figure('Visible', 'on');
    hold on;

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

%% Histogram
% Choosing a few key timepoints; helps to identify potential outliers
% x-axis will display mean ERP amplitudes at specified timepoint
% y-axis shows how many participants displayed that amplitude

% timepoints approximately related to the N2 component

N2_TPs = [200, 220, 250];

amps_atN2TP = groupDat(:, N2_TPs(1));
n = height(amps_atN2TP); 
bin = ceil(log2(n) +1); % applies Sturge's formula to determine appropriate number of bins


for chan = 1:length(channels)
    % Get channel name
    chanName = channels(chan);
    
    % Extract data for that channel
    chanDat = allData(strcmp(allData.Channel, chanName), :);

    for TP = 1:length(N2_TPs)
        figure;

        Time = N2_TPs(TP);
        colName = "Time" + N2_TPs(TP);
        amps_atN2TP = chanDat.(colName);
        n = height(amps_atN2TP); 
        bin = ceil(log2(n) +1); % applies Sturge's formula to determine appropriate number of bins

        histogram(amps_atN2TP, bin);
        title(['Histogram - ', char(chanName), ' @ ', num2str(Time), ' ms']);
        xlabel('Amplitude (µV)');
        ylabel('Subject Count');
        ylim([0 15]);
    end
end

% timepoints approximately related to the N2 component

P3_TPs = [300, 375, 450];
for chan = 1:length(channels)
    % Get channel name
    chanName = channels(chan);
    
    % Extract data for that channel
    chanDat = allData(strcmp(allData.Channel, chanName), :);

    for TP = 1:length(P3_TPs)
        figure;

        Time = P3_TPs(TP);
        colName = "Time" + P3_TPs(TP);
        amps_atN2TP = chanDat.(colName);
        n = height(amps_atN2TP); 
        bin = ceil(log2(n) +1); % applies Sturge's formula to determine appropriate number of bins

        histogram(amps_atN2TP, bin, ...
            'FaceColor', [1, 0.4, 0.6]);
        title(['Histogram - ', char(chanName), ' @ ', num2str(Time), ' ms']);
        xlabel('Amplitude (µV)');
        ylabel('Subject Count');
        ylim([0 15]);
    end
end