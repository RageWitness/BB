function PL_dB = compute_pathloss_from_cache(g, freq_hz, GeometryCache)
% COMPUTE_PATHLOSS_FROM_CACHE  从几何缓存 + 频点参数计算 (M x 1) 路径损耗
%
%   输入：
%     g             - 网格点全局索引
%     freq_hz       - 频率
%     GeometryCache - 由 precompute_geometry_cache 生成
%
%   输出：
%     PL_dB - (M x 1) 路径损耗

    p = get_channel_params(freq_hz);

    d_row      = GeometryCache.d(g, :)';        % M x 1
    Nw_row     = GeometryCache.Nw(g, :)';
    Nb_row     = GeometryCache.Nb(g, :)';
    dc_row     = GeometryCache.dc(g, :)';
    Nc_row     = GeometryCache.Nc(g, :)';
    Nwdif_row  = GeometryCache.Nw_dif(g, :)';
    diff_ok    = GeometryCache.diff_ok(g, :)';

    M = numel(d_row);

    % 多斜率函数
    d_safe  = max(d_row, p.d0);
    PL_ms_d = pl_multislope(d_safe, p);

    PL_direct = PL_ms_d + Nw_row * p.Lw;

    % 绕射
    PL_diff = inf(M, 1);
    has_diff = diff_ok & (Nw_row > 0);
    if any(has_diff)
        dc_safe = max(dc_row(has_diff), p.d0);
        PL_diff(has_diff) = pl_multislope(dc_safe, p) ...
            + Nc_row(has_diff)    * p.Lc ...
            + Nwdif_row(has_diff) * p.Lw ...
            + max(Nb_row(has_diff) - 1, 0) * p.Deltaic;
    end

    PL_dB = PL_direct;
    los_mask = (Nw_row == 0);
    PL_dB(los_mask) = PL_ms_d(los_mask);

    diff_better = (Nw_row > 0) & (PL_diff < PL_direct);
    PL_dB(diff_better) = PL_diff(diff_better);
end


function pl = pl_multislope(d, p)
    pl = p.PL0 ...
         + p.N1 * log10(min(d, p.db) / p.d0) ...
         + p.N2 * log10(max(d, p.db) / p.db);
end
