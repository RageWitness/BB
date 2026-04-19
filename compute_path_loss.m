function out = compute_path_loss(src, ap, freq_hz, buildings, opts)
% COMPUTE_PATH_LOSS  统一信道模型 (多斜率 + MWM + 简化绕射)
%
%   out.PL_total      - 最终路径损耗 (dB)
%   out.PL_direct     - 直达路径总损耗
%   out.PL_diff       - 简化绕射路径总损耗 (无效时为 inf)
%   out.PL_ms_direct  - 直达多斜率部分
%   out.chosen_mode   - 'los' | 'direct' | 'diffraction'
%   out.wall_count    - Nw_i (直达穿墙数)
%   out.corner_count  - Nc_i (绕射拐角数)
%   out.recv_power    - Pt + Gt + Gi - PL_total (若提供 opts.tx_power)
%   out.params        - 调用的频点参数

    if nargin < 5, opts = struct(); end
    if ~isfield(opts, 'tx_power'), opts.tx_power = []; end
    if ~isfield(opts, 'tx_gain'),  opts.tx_gain  = 0; end
    if ~isfield(opts, 'rx_gain'),  opts.rx_gain  = 0; end
    if ~isfield(opts, 'VG'),       opts.VG       = []; end

    p = get_channel_params(freq_hz);

    src = src(:)';
    ap  = ap(:)';
    d_i = norm(ap - src);
    d_i_safe = max(d_i, p.d0);

    PL_ms_direct = pl_multislope(d_i_safe, p);

    [Nw_i, Nb_i] = count_wall_crossings(src, ap, buildings);
    PL_direct = PL_ms_direct + Nw_i * p.Lw;

    PL_diff = inf;
    Nc_i = 0;
    if Nw_i > 0
        pi_d = find_diffraction_path(src, ap, buildings, opts.VG);
        if pi_d.ok
            dc = max(pi_d.total_length, p.d0);
            PL_ms_dc = pl_multislope(dc, p);
            PL_diff = PL_ms_dc ...
                    + pi_d.num_corners * p.Lc ...
                    + pi_d.wall_crossings_along_path * p.Lw ...
                    + max(Nb_i - 1, 0) * p.Deltaic;
            Nc_i = pi_d.num_corners;
        end
    end

    if Nw_i == 0
        out.PL_total = PL_ms_direct;
        out.chosen_mode = 'los';
    else
        if PL_diff < PL_direct
            out.PL_total = PL_diff;
            out.chosen_mode = 'diffraction';
        else
            out.PL_total = PL_direct;
            out.chosen_mode = 'direct';
        end
    end

    out.PL_direct    = PL_direct;
    out.PL_diff      = PL_diff;
    out.PL_ms_direct = PL_ms_direct;
    out.wall_count   = Nw_i;
    out.corner_count = Nc_i;
    out.params       = p;

    if ~isempty(opts.tx_power)
        out.recv_power = opts.tx_power + opts.tx_gain + opts.rx_gain - out.PL_total;
    else
        out.recv_power = NaN;
    end
end


function pl = pl_multislope(d, p)
    pl = p.PL0 ...
         + p.N1 * log10(min(d, p.db) / p.d0) ...
         + p.N2 * log10(max(d, p.db) / p.db);
end
