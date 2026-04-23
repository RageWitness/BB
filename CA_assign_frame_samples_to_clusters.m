function CADataset = CA_assign_frame_samples_to_clusters(CADataset, Partition, ~)
% CA_ASSIGN_FRAME_SAMPLES_TO_CLUSTERS  帧样本按位置到 cluster center 归簇

    centers = Partition.centers;
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
