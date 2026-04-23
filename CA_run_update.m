function [OnlineSpatialFP_new, CAResult] = CA_run_update(OnlineSpatialFP, EventList, Config)
% CA_RUN_UPDATE  CA 总入口：样本收集 → 残差 → KMGPR → lambda → 回写

    Config = CA_fill_defaults(Config);

    CAResult = struct();
    CAResult.config = Config.ca;

    % --- 1. 样本收集 ---
    CASamples = CA_collect_reliable_samples(EventList, Config);
    n_total = numel(CASamples);
    if n_total > 0
        n_valid = sum([CASamples.valid]);
    else
        n_valid = 0;
    end
    CAResult.n_samples_total = n_total;
    CAResult.n_samples_valid = n_valid;
    CAResult.samples = CASamples;

    % 空样本安全返回
    if n_valid == 0
        OnlineSpatialFP_new = OnlineSpatialFP;
        CAResult.status = 'no_valid_ca_samples';
        CAResult.dataset   = [];
        CAResult.partition = [];
        CAResult.models    = [];
        CAResult.delta     = [];
        CAResult.update_log = struct( ...
            'n_updated_entries', 0, ...
            'mean_abs_delta', 0, ...
            'max_abs_delta', 0, ...
            'status', 'no_valid_ca_samples');
        fprintf('[CA] 无有效样本，跳过更新\n');
        return;
    end

    % --- 2. 残差数据集 ---
    CADataset = CA_build_residual_dataset(CASamples, OnlineSpatialFP, Config);

    % --- 3. K-Means 分区 ---
    Partition = CA_kmeans_partition_grid(OnlineSpatialFP.grid_xy, Config);

    % --- 4. 样本归簇 ---
    CADataset = CA_assign_samples_to_clusters(CADataset, Partition, Config);

    % --- 5. GPR 训练 ---
    Models = CA_train_kmgpr_models(CADataset, Partition, Config);

    % --- 6. 预测 delta ---
    DeltaA = CA_predict_kmgpr_delta(Models, CADataset, OnlineSpatialFP, Partition, Config);

    % --- 7. lambda map ---
    DeltaA = CA_compute_lambda_map(DeltaA, Config);

    % --- 8. 回写 ---
    [OnlineSpatialFP_new, UpdateLog] = CA_apply_delta_to_online_map(OnlineSpatialFP, DeltaA, Config);

    CAResult.status     = 'ok';
    CAResult.dataset    = CADataset;
    CAResult.partition  = Partition;
    CAResult.models     = Models;
    CAResult.delta      = DeltaA;
    CAResult.update_log = UpdateLog;
end
