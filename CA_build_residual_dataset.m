function CADataset = CA_build_residual_dataset(CASamples, OnlineSpatialFP, Config)
% CA_BUILD_RESIDUAL_DATASET  用 IDW 插值旧图，构造残差 Y = obs_centered - old_map
%
%   CADataset.band(b).X:  N x 2  样本位置
%   CADataset.band(b).Y:  N x M  残差（允许 NaN）

    Config = CA_fill_defaults(Config);
    B = OnlineSpatialFP.B;
    M = OnlineSpatialFP.M;
    grid_xy = OnlineSpatialFP.grid_xy;

    for b = 1:B
        CADataset.band(b).X          = [];
        CADataset.band(b).Y          = [];
        CADataset.band(b).grid_idx   = [];
        CADataset.band(b).sample_ids = [];
        CADataset.band(b).source_uid = {};
        CADataset.band(b).cluster_id = [];
    end

    valid_mask = false(1, numel(CASamples));
    for i = 1:numel(CASamples)
        valid_mask(i) = CASamples(i).valid;
    end
    valid_idx = find(valid_mask);

    if isempty(valid_idx)
        return;
    end

    for b = 1:B
        xs = []; ys = []; sids = []; uids = {};
        for j = 1:numel(valid_idx)
            s = CASamples(valid_idx(j));
            if s.band_id ~= b, continue; end

            x_r = s.position_xy;
            z_obs = s.obs_centered_dBm;  % M x 1

            z_old = idw_interpolate(grid_xy, ...
                OnlineSpatialFP.band(b).centered_dBm, x_r, Config);  % M x 1

            residual = z_obs - z_old;  % NaN where obs invalid

            xs   = [xs; x_r];          %#ok<AGROW>
            ys   = [ys; residual'];    %#ok<AGROW>
            sids = [sids; s.sample_id]; %#ok<AGROW>
            uids{end+1} = s.source_uid; %#ok<AGROW>
        end

        CADataset.band(b).X          = xs;
        CADataset.band(b).Y          = ys;       % N x M, NaN ok
        CADataset.band(b).sample_ids = sids;
        CADataset.band(b).source_uid = uids(:);
        CADataset.band(b).cluster_id = zeros(size(sids));
    end

    n_total = sum(arrayfun(@(bd) size(bd.X, 1), CADataset.band));
    fprintf('[CA] 残差数据集: %d 个有效样本 across %d bands\n', n_total, B);
end


function z_interp = idw_interpolate(grid_xy, Z_map, x_r, Config)
% IDW_INTERPOLATE  反距离加权插值
%   grid_xy: G x 2,  Z_map: M x G,  x_r: 1 x 2
%   返回 M x 1

    interp_cfg = Config.ca.interp;
    mode = interp_cfg.mode;
    M = size(Z_map, 1);

    if strcmp(mode, 'nearest')
        dists = sqrt(sum((grid_xy - x_r).^2, 2));
        [~, idx] = min(dists);
        z_interp = Z_map(:, idx);
        return;
    end

    % IDW (default)
    K_i   = interp_cfg.k;
    p_exp = interp_cfg.power;
    eps_v = interp_cfg.eps;

    dists = sqrt(sum((grid_xy - x_r).^2, 2));

    if min(dists) < eps_v
        [~, idx] = min(dists);
        z_interp = Z_map(:, idx);
        return;
    end

    [~, sort_idx] = sort(dists);
    nn_idx = sort_idx(1:min(K_i, numel(sort_idx)));
    d_nn = dists(nn_idx);

    w = 1 ./ (d_nn .^ p_exp + eps_v);
    w = w / sum(w);

    z_interp = Z_map(:, nn_idx) * w;  % M x K * K x 1 = M x 1
end
