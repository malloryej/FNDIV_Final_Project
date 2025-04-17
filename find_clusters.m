% Creates function to identify clusters

function clusterResults = find_clusters(strucResults, alpha)
    % Output structure
    clusterResults = struct('Channel', {}, ...
                            'Timepoint_start', {}, ...
                            'Timepoint_end', {}, ...
                            'Cluster_length', {});

    clusterIDX = 1;
    i = 1;

    while i <= length(strucResults)
        % Start of potential cluster
        if strucResults(i).PVal <= alpha
            currentChannel = strucResults(i).Channel;
            startTime = strucResults(i).Timepoint;
            endTime = startTime;

            % Continue while conditions are met
            j = i + 1;
            while j <= length(strucResults) && ...
                  strucResults(j).PVal <= alpha && ...
                  strucResults(j).Channel == currentChannel && ...
                  strucResults(j).Timepoint == endTime + 1

                endTime = strucResults(j).Timepoint;
                j = j + 1;
            end

            % Save cluster
            clusterResults(clusterIDX).Channel = currentChannel;
            clusterResults(clusterIDX).Timepoint_start = startTime;
            clusterResults(clusterIDX).Timepoint_end = endTime;
            clusterResults(clusterIDX).Cluster_length = endTime - startTime + 1;

            clusterIDX = clusterIDX + 1;
            i = j;  % continue from after this cluster
        else
            i = i + 1;
        end
    end
end
