%% Create function to Identify Clusters in Data
function listlengths = find_clusters(significant_mask)
    listlengths = [];
    [num_rows, num_cols] = size(significant_mask);

    for row = 1:num_rows
        cluster_length = 0;

        for col = 1:num_cols
            if significant_mask(col) == 1
                cluster_length = cluster_length + 1;

            elseif cluster_length > 0 % only add to listlengths if cluster length > 0
                % ends current cluster; resets variables
                listlengths = [listlengths, cluster_length];
                cluster_length = 0;

            end
        end

        if cluster_length > 0
            listlengths = [listlengths, cluster_length];
            % handle cases where last timepoint is still in cluster

        end
    end
end