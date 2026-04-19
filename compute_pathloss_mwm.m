function [PL_dB, sigma_dB] = compute_pathloss_mwm(src_xy, ap_xy_all, freq_hz, buildings)
% COMPUTE_PATHLOSS_MWM  对所有 AP 批量计算统一 MWM 路径损耗
%
%   src_xy     - (1x2) 源
%   ap_xy_all  - (M x 2) 所有 AP
%   freq_hz    - 频率
%   buildings  - 建筑列表
%
%   返回：
%     PL_dB    - (M x 1) 路径损耗
%     sigma_dB - 频点对应的残差 sigma

    M = size(ap_xy_all, 1);
    PL_dB = zeros(M, 1);
    p = get_channel_params(freq_hz);
    sigma_dB = p.sigma;

    for m = 1:M
        out = compute_path_loss(src_xy, ap_xy_all(m, :), freq_hz, buildings, struct());
        PL_dB(m) = out.PL_total;
    end
end
