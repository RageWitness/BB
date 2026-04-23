function [OnlineSpatialFP_new, CAResult] = CA_run_update_frames( ...
    OnlineSpatialFP, SourceContext, Y_lin_all, Y_dBm_all, Config)
% CA_RUN_UPDATE_FRAMES  按帧 CA 主入口
%
%   SourceContext + Y_lin_all → 帧样本 → 残差 → 加权 KMGPR → lambda → 回写

    Config = CA_fill_defaults(Config);

    CAResult = struct();
    CAResult.config = Config.ca;

    % --- 1. 帧样本收集 ---
    CAFrameSamples = CA_collect_reliable_frame_samples( ...
        SourceContext, Y_lin_all, Y_dBm_all, Config);

    n_total = numel(CAFrameSamples);
    if n_total > 0
        n_valid = sum([CAFrameSamples.valid]);
    else
        n_valid = 0;
    end
    CAResult.n_samples_total = n_total;
    CAResult.n_samples_valid = n_valid;
    CAResult.frame_samples   = CAFrameSamples;

    % per-band n_eff
    [~, B, ~] = size(Y_lin_all);
    n_eff_per_band = zeros(1, B);
    if n_valid > 0
        for b = 1:B
            bm = [CAFrameSamples.band_id] == b & [CAFrameSamples.valid];
            if any(bm)
                ws = [CAFrameSamples(bm).sample_weight];
                n_eff_per_band(b) = sum(ws)^2 / sum(ws.^2);
            end
        end
    end
    CAResult.n_eff_per_band = n_eff_per_band;

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
        fprintf('[CA-Frame] 无有效帧样本，跳过更新\n');
        return;
    end

    % --- 2. 残差数据集 ---
    CADataset = CA_build_frame_residual_dataset( ...
        CAFrameSamples, OnlineSpatialFP, Config);

    % --- 3. K-Means 分区 ---
    Partition = CA_kmeans_partition_grid(OnlineSpatialFP.grid_xy, Config);

    % --- 4. 样本归簇 ---
    CADataset = CA_assign_frame_samples_to_clusters(CADataset, Partition, Config);

    % --- 5. 加权 GPR 训练 ---
    Models = CA_train_weighted_kmgpr_models(CADataset, Partition, Config);

    % --- 6. 预测 delta ---
    DeltaA = CA_predict_kmgpr_delta(Models, CADataset, OnlineSpatialFP, Partition, Config);

    % --- 7. lambda map ---
    DeltaA = CA_compute_lambda_map(DeltaA, Config);

    % --- 8. 回写 ---
    [OnlineSpatialFP_new, UpdateLog] = CA_apply_delta_to_online_map( ...
        OnlineSpatialFP, DeltaA, Config);

    CAResult.status     = 'ok';
    CAResult.dataset    = CADataset;
    CAResult.partition  = Partition;
    CAResult.models     = Models;
    CAResult.delta      = DeltaA;
    CAResult.update_log = UpdateLog;
end
