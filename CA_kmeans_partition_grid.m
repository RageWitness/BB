function Partition = CA_kmeans_partition_grid(grid_xy, Config)
% CA_KMEANS_PARTITION_GRID  对网格坐标做 K-Means 分区
%
%   Partition.cluster_id  G x 1
%   Partition.centers     K x 2
%   Partition.members{k}  该 cluster 的网格索引

    Config = CA_fill_defaults(Config);
    K = Config.ca.kmeans.K;
    max_iter = Config.ca.kmeans.max_iter;
    n_rep = Config.ca.kmeans.n_replicates;
    G = size(grid_xy, 1);

    if G <= K
        Partition.cluster_id = (1:G)';
        Partition.centers = grid_xy;
        Partition.K = G;
        for k = 1:G
            Partition.members{k} = k;
        end
        return;
    end

    best_cost = Inf;
    best_cid = [];
    best_centers = [];

    for rep = 1:n_rep
        perm = randperm(G, K);
        centers = grid_xy(perm, :);

        for iter = 1:max_iter
            dists = zeros(G, K);
            for k = 1:K
                diffs = grid_xy - centers(k, :);
                dists(:, k) = sum(diffs.^2, 2);
            end
            [~, cid] = min(dists, [], 2);

            new_centers = zeros(K, 2);
            for k = 1:K
                mask = (cid == k);
                if any(mask)
                    new_centers(k, :) = mean(grid_xy(mask, :), 1);
                else
                    new_centers(k, :) = centers(k, :);
                end
            end

            if max(abs(new_centers(:) - centers(:))) < 1e-8
                centers = new_centers;
                break;
            end
            centers = new_centers;
        end

        cost = 0;
        for k = 1:K
            mask = (cid == k);
            if any(mask)
                diffs = grid_xy(mask, :) - centers(k, :);
                cost = cost + sum(sum(diffs.^2));
            end
        end

        if cost < best_cost
            best_cost = cost;
            best_cid = cid;
            best_centers = centers;
        end
    end

    Partition.cluster_id = best_cid;
    Partition.centers = best_centers;
    Partition.K = K;
    for k = 1:K
        Partition.members{k} = find(best_cid == k);
    end

    fprintf('[CA] K-Means: K=%d, G=%d, cost=%.1f\n', K, G, best_cost);
end
