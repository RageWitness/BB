function PL_dB = compute_pathloss_cost231wi(dist_m, fc_Hz, params)
% COMPUTE_PATHLOSS_COST231WI  COST231 Walfisch-Ikegami 路径损耗模型
%
%   PL_dB = compute_pathloss_cost231wi(dist_m, fc_Hz, params)
%
%   适用频率范围：800 MHz ~ 2000 MHz（可外推至 2.6 GHz）
%   适用距离范围：20 m ~ 5000 m
%
%   输入：
%       dist_m  - (M x 1) 源到各 AP 的距离 (m)
%       fc_Hz   - 载频 (Hz)
%       params  - 环境参数结构体：
%                   .h_base_m      AP/基站天线高度 (m)
%                   .h_mobile_m    移动端/辐射源高度 (m)
%                   .h_roof_m      建筑物平均屋顶高度 (m)
%                   .w_street_m    街道宽度 (m)
%                   .w_building_m  建筑物间距 (m)
%                   .phi_deg       街道与传播方向夹角 (度)
%                   .city_type     'medium' 或 'metropolitan'
%
%   输出：
%       PL_dB   - (M x 1) 路径损耗 (dB)
%
%   参考：COST231 Final Report, Chapter 4.4 (Walfisch-Ikegami model)

    f_MHz = fc_Hz / 1e6;
    d_km  = dist_m / 1e3;

    h_b   = params.h_base_m;       % 基站天线高度
    h_m   = params.h_mobile_m;     % 移动端高度
    h_r   = params.h_roof_m;       % 平均屋顶高度
    w     = params.w_street_m;     % 街道宽度
    b     = params.w_building_m;   % 建筑物间距
    phi   = params.phi_deg;        % 道路方向角 (度)

    M = numel(dist_m);
    PL_dB = zeros(M, 1);

    % 距离下限保护
    d_km = max(d_km, 0.02);  % 最小 20m

    for i = 1:M
        d = d_km(i);

        % ===== 判断 LOS / NLOS =====
        % 简化判断：若距离 < 某阈值且基站高于屋顶，认为 LOS
        if h_b > h_r && d < 0.05  % 50m 内且基站高于屋顶
            % --- LOS 模型 ---
            PL_dB(i) = 42.6 + 26 * log10(d) + 20 * log10(f_MHz);
        else
            % --- NLOS 模型 ---

            % 1. 自由空间损耗 L0
            L0 = 32.4 + 20 * log10(d) + 20 * log10(f_MHz);

            % 2. 屋顶-街道绕射与散射损耗 Lrts
            Lori = compute_Lori(phi);
            delta_h_m = h_r - h_m;
            Lrts = -16.9 - 10 * log10(w) + 10 * log10(f_MHz) + ...
                   20 * log10(delta_h_m) + Lori;

            % 3. 多屏绕射损耗 Lmsd
            delta_h_b = h_b - h_r;

            % Lbsh
            if h_b > h_r
                Lbsh = -18 * log10(1 + delta_h_b);
            else
                Lbsh = 0;
            end

            % ka
            if h_b > h_r
                ka = 54;
            elseif d >= 0.5
                ka = 54 - 0.8 * abs(delta_h_b);
            else
                ka = 54 - 0.8 * abs(delta_h_b) * d / 0.5;
            end

            % kd
            if h_b > h_r
                kd = 18;
            else
                kd = 18 - 15 * abs(delta_h_b) / h_r;
            end

            % kf
            if strcmp(params.city_type, 'metropolitan')
                kf = -4 + 1.5 * (f_MHz / 925 - 1);
            else
                kf = -4 + 0.7 * (f_MHz / 925 - 1);
            end

            Lmsd = Lbsh + ka + kd * log10(d) + kf * log10(f_MHz) - ...
                   9 * log10(b);

            % 4. 组合
            if Lrts + Lmsd > 0
                PL_dB(i) = L0 + Lrts + Lmsd;
            else
                PL_dB(i) = L0;
            end
        end
    end
end


function Lori = compute_Lori(phi)
% COMPUTE_LORI  街道方向修正因子
    if phi >= 0 && phi < 35
        Lori = -10 + 0.354 * phi;
    elseif phi >= 35 && phi < 55
        Lori = 2.5 + 0.075 * (phi - 35);
    else
        Lori = 4.0 - 0.114 * (phi - 55);
    end
end
