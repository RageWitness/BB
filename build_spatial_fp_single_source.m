function SpatialFP = build_spatial_fp_single_source(APs, Bands, GridValid, Config)
% BUILD_SPATIAL_FP_SINGLE_SOURCE  构建每个频带的单源空间指纹库
%
%   SpatialFP = build_spatial_fp_single_source(APs, Bands, GridValid, Config)
%
%   对每个频带 b、每个有效参考点 g，放置参考单源（0 dBm），
%   用信道模型（含 Monte Carlo 阴影采样）计算各 AP 的 RSS 指纹。
%
%   输出：
%       SpatialFP.B                    - 频带数
%       SpatialFP.M                    - AP 数
%       SpatialFP.G                    - 有效网格点数
%       SpatialFP.grid_xy              - (G x 2)
%       SpatialFP.ref_power_dBm        - (1 x B) 各频带参考发射功率
%       SpatialFP.band(b).F_dBm        - (M x G)
%       SpatialFP.band(b).F_lin        - (M x G)
%       SpatialFP.band(b).std_dBm      - (M x G)
%       SpatialFP.band(b).std_lin      - (M x G)
%       SpatialFP.band(b).mean_dBm     - (1 x G)
%       SpatialFP.band(b).centered_dBm - (M x G)
%       SpatialFP.band(b).mean_lin     - (1 x G)
%       SpatialFP.band(b).centered_lin - (M x G)
%       SpatialFP.band(b).fc_Hz
%       SpatialFP.band(b).bw_Hz
%       SpatialFP.band(b).model

    B = Bands.B;
    M = APs.num;
    G = GridValid.Nvalid;

    % 默认参考功率
    if isfield(Config, 'm25') && isfield(Config.m25, 'ref_power_dBm')
        ref_power = Config.m25.ref_power_dBm;
    else
        ref_power = 100 * ones(1, B);  % 默认 100 dBm
    end

    SpatialFP.B = B;
    SpatialFP.M = M;
    SpatialFP.G = G;
    SpatialFP.grid_xy       = GridValid.xy;
    SpatialFP.ref_power_dBm = ref_power;

    fprintf('[M2.5] 构建 SpatialFP: M=%d AP, G=%d 网格点, B=%d 频带\n', M, G, B);

    for b = 1:B
        F_dBm   = zeros(M, G);
        F_lin   = zeros(M, G);
        Std_dBm = zeros(M, G);
        Std_lin = zeros(M, G);

        for g = 1:G
            [rss_dBm, rss_lin, rss_std_dBm, rss_std_lin] = ...
                simulate_reference_fp_point_m1bridge( ...
                    GridValid.xy(g, :), b, ref_power(b), APs, Bands, Config);
            F_dBm(:, g)   = rss_dBm;
            F_lin(:, g)   = rss_lin;
            Std_dBm(:, g) = rss_std_dBm;
            Std_lin(:, g) = rss_std_lin;
        end

        % 组装单频带结构
        band_fp.F_dBm   = F_dBm;
        band_fp.F_lin   = F_lin;
        band_fp.std_dBm = Std_dBm;
        band_fp.std_lin = Std_lin;
        band_fp.fc_Hz   = Bands.fc_Hz(b);
        band_fp.bw_Hz   = Bands.bw_Hz(b);
        band_fp.model   = Bands.model{b};

        % 预计算均值和中心化指纹（dBm + 线性 两域）
        band_fp = precompute_centered_fp_wknn(band_fp);

        SpatialFP.band(b) = band_fp;

        fprintf('  Band %d (%s): F_dBm [%.1f, %.1f], F_lin [%.2e, %.2e]\n', ...
            b, Bands.name{b}, min(F_dBm(:)), max(F_dBm(:)), ...
            min(F_lin(:)), max(F_lin(:)));
    end

    fprintf('[M2.5] SpatialFP 构建完成\n');
end
