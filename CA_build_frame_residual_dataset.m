function CADataset = CA_build_frame_residual_dataset( ...
    CAFrameSamples, OnlineSpatialFP, Config)
% CA_BUILD_FRAME_RESIDUAL_DATASET  从帧样本构造残差数据集（含权重 W）
%
%   CADataset.band(b).X:  N x 2
%   CADataset.band(b).Y:  N x M  (NaN ok)
%   CADataset.band(b).W:  N x 1
%   CADataset.band(b).frame_id:   N x 1
%   CADataset.band(b).source_uid: cell N x 1
%   CADataset.band(b).cluster_id: N x 1

    Config = CA_fill_defaults(Config);
    B = OnlineSpatialFP.B;
    grid_xy = OnlineSpatialFP.grid_xy;

    for b = 1:B
        CADataset.band(b).X          = [];
        CADataset.band(b).Y          = [];
        CADataset.band(b).W          = [];
        CADataset.band(b).frame_id   = [];
        CADataset.band(b).source_uid = {};
        CADataset.band(b).cluster_id = [];
    end

    valid_mask = false(1, numel(CAFrameSamples));
    for i = 1:numel(CAFrameSamples)
        valid_mask(i) = CAFrameSamples(i).valid;
    end
    valid_idx = find(valid_mask);

    if isempty(valid_idx), return; end

    for b = 1:B
        xs = []; ys = []; ws = []; fids = []; uids = {};
        for j = 1:numel(valid_idx)
            s = CAFrameSamples(valid_idx(j));
            if s.band_id ~= b, continue; end

            x_r = s.position_xy;
            z_obs = s.obs_centered_dBm;  % M x 1

            z_old = idw_interpolate(grid_xy, ...
                OnlineSpatialFP.band(b).centered_dBm, x_r, Config);

            residual = z_obs - z_old;  % NaN where obs invalid

            xs   = [xs; x_r];              %#ok<AGROW>
            ys   = [ys; residual'];        %#ok<AGROW>
            ws   = [ws; s.sample_weight];   %#ok<AGROW>
            fids = [fids; s.frame_id];      %#ok<AGROW>
            uids{end+1} = s.source_uid;     %#ok<AGROW>
        end

        CADataset.band(b).X          = xs;
        CADataset.band(b).Y          = ys;
        CADataset.band(b).W          = ws;
        CADataset.band(b).frame_id   = fids;
        CADataset.band(b).source_uid = uids(:);
        CADataset.band(b).cluster_id = zeros(size(ws));
    end

    n_total = sum(arrayfun(@(bd) size(bd.X, 1), CADataset.band));
    fprintf('[CA-Frame] 残差数据集: %d 个帧样本 across %d bands\n', n_total, B);
end


function z_interp = idw_interpolate(grid_xy, Z_map, x_r, Config)
    interp_cfg = Config.ca.interp;
    mode = interp_cfg.mode;

    if strcmp(mode, 'nearest')
        dists = sqrt(sum((grid_xy - x_r).^2, 2));
        [~, idx] = min(dists);
        z_interp = Z_map(:, idx);
        return;
    end

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

    z_interp = Z_map(:, nn_idx) * w;
end
