function params = get_channel_params(freq_hz)
% GET_CHANNEL_PARAMS  返回统一 MWM 信道模型的频点参数
%
%   输出字段：
%     PL0, Lw, Lc, Deltaic, sigma, db, N1, N2, d0, f_MHz

    f_MHz = freq_hz / 1e6;

    % 全局结构参数（所有频点共用）/原始N1 N2：20，40
    params.d0 = 1;
    params.db = 20;
    params.N1 = 20;
    params.N2 = 30;
    params.f_MHz = f_MHz;

    % PL0(f) = 20*log10(f_MHz) - 28
    params.PL0 = 20 * log10(f_MHz) - 28;

    % 频点相关参数表（按 f_MHz 匹配，容差 0.5 MHz）freq,lw,lc,delta,sigma
        % 96,   3.0, 5.0, 2.0, 4.0;
        % 89,   2.8, 4.8, 1.8, 4.0;
        % 2400, 6.0, 10.0, 4.0, 6.0;
        % 2593, 6.4, 10.5, 4.2, 6.5 下面这个为wl版
    table = {
        96,   3.6, 5.4, 2.0, 4.0;
        89,   3.0, 5.8, 1.8, 4.0;
        2400, 6.8, 10.0, 4.0, 6.0;
        2593, 7, 10.5, 4.2, 6.5
    };

    matched = false;
    for k = 1:size(table, 1)
        if abs(f_MHz - table{k,1}) < 0.5
            params.Lw      = table{k, 2};
            params.Lc      = table{k, 3};
            params.Deltaic = table{k, 4};
            params.sigma   = table{k, 5};
            matched = true;
            break;
        end
    end

    if ~matched
        % 兜底：按低/高频粗分
        if f_MHz < 500
            params.Lw = 3.0; params.Lc = 5.0; params.Deltaic = 2.0; params.sigma = 4.0;
        else
            params.Lw = 6.0; params.Lc = 10.0; params.Deltaic = 4.0; params.sigma = 6.0;
        end
    end
end
