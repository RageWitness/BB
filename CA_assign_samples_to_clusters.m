function CADataset = CA_assign_samples_to_clusters(CADataset, Partition, ~)
% CA_ASSIGN_SAMPLES_TO_CLUSTERS  按样本位置到 cluster center 距离归簇

    centers = Partition.centers;
    K = size(centers, 1);
    B = numel(CADataset.band);

    for b = 1:B
        X = CADataset.band(b).X;
        N = size(X, 1);
        if N == 0, continue; end

        cid = zeros(N, 1);
        for i = 1:N
            dists = sum((centers - X(i, :)).^2, 2);
            [~, cid(i)] = min(dists);
        end
        CADataset.band(b).cluster_id = cid;
    end
end
