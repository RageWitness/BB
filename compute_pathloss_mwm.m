function [PL_dB, sigma_dB] = compute_pathloss_mwm(src_xy, ap_xy_all, freq_hz, buildings)
% COMPUTE_PATHLOSS_MWM  对所有 AP 批量计算统一 MWM 路径损耗
%   使用 persistent 静态角点可见图缓存，避免每次重建

    persistent VG_cache bld_sig
    cur_sig = get_buildings_signature(buildings);
    if isempty(VG_cache) || ~isequal(bld_sig, cur_sig)
        VG_cache = build_static_corner_graph(buildings);
        bld_sig = cur_sig;
    end

    M = size(ap_xy_all, 1);
    PL_dB = zeros(M, 1);
    p = get_channel_params(freq_hz);
    sigma_dB = p.sigma;

    for m = 1:M
        out = compute_path_loss(src_xy, ap_xy_all(m, :), freq_hz, buildings, ...
            struct('VG', VG_cache));
        PL_dB(m) = out.PL_total;
    end
end


function sig = get_buildings_signature(buildings)
    if isempty(buildings)
        sig = [0 0 0 0];
        return;
    end
    sig = zeros(numel(buildings), 4);
    for k = 1:numel(buildings)
        sig(k, :) = [buildings(k).xmin, buildings(k).xmax, ...
                     buildings(k).ymin, buildings(k).ymax];
    end
end

