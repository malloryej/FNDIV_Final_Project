%% Create function to Identify Clusters in Data
function clusterResults = find_clusters(strucResults, alpha, min_size)
    % Output structure
    clusterResults = struct('Channel', {}, 'Timepoint_start', {}, 'Timepoint_end', {});

    % Get unique channels
    channels = unique({strucResults.Channel});
    clusterIdx = 1;

    for ch = 1:length(channels)
        thisChannel = channels{ch};

        % Filter results for this channel
        chResults = strucResults(strcmp({strucResults.Channel}, thisChannel));

        % Sort by timepoint (optional but safer)
        [~, sortIdx] = sort([chResults.Timepoint]);
        chResults = chResults(sortIdx);

        % Create logical vector of significance
        sigVec = [chResults.PVal] < alpha;
        tPoints = [chResults.Timepoint];

        % Find clusters of contiguous 1s
        t = 1;
        while t <= length(sigVec)
            if sigVec(t)
                startT = t;
                while t <= length(sigVec) && sigVec(t)
                    t = t + 1;
                end
                endT = t - 1;

                % Check if cluster length meets minimum size
                clusterLength = tPoints(endT) - tPoints(startT) + 1; % assumes timepoints are in ms and consecutive
                if clusterLength >= min_size
                    clusterResults(clusterIdx).Channel = thisChannel;
                    clusterResults(clusterIdx).Timepoint_start = tPoints(startT);
                    clusterResults(clusterIdx).Timepoint_end = tPoints(endT);
                    clusterIdx = clusterIdx + 1;
                end
            else
                t = t + 1;
            end
        end
    end
end
