function EventList = score_event_type_m3(EventList, Config)
% SCORE_EVENT_TYPE_M3  对每个事件计算四类源评分并判定类型
%
%   EventList = score_event_type_m3(EventList, Config)
%
%   评分体系（显式评分器）：
%     A. trusted_fixed:  频带覆盖 + 高功率 + 稳定性 + 位置匹配
%     B. prior_pos_known: 时间匹配 + 位置匹配 + 功率近0dBm + 稳定性
%     C. prior_time_known: 时间匹配 + 功率近0dBm + 稳定性
%     D. ordinary_target: 残差项 + 目标功率范围 + 持续时间
%
%   输入/输出：
%       EventList(e) — 增加 score_trusted, score_prior_pos, score_prior_time,
%                      score_target, type_hat 字段

    % --- 填充评分参数默认值 ---
    sc = fill_score_defaults(Config);

    n_events = numel(EventList);

    for e = 1:n_events
        ev = EventList(e);

        %% A. trusted_fixed 分数
        % 频带覆盖分数
        S_cover = sum(ev.band_coverage_vec) / numel(ev.band_coverage_vec);

        % 高功率分数：与 trusted 标称功率 15 dBm 的接近度
        S_power_high = exp(-((ev.power_level_est - sc.trusted_power_ref)^2) / ...
                           (sc.sigma_trusted_power^2));

        % 稳定性分数
        S_stable = exp(-(ev.power_stability_est) / sc.kappa_trusted);

        % 位置匹配分数（直接用已提取的）
        S_pos_anchor = ev.position_prior_match_score;

        ev.score_trusted = sc.alpha(1) * S_cover + ...
                           sc.alpha(2) * S_power_high + ...
                           sc.alpha(3) * S_stable + ...
                           sc.alpha(4) * S_pos_anchor;

        %% B. prior_pos_known 分数
        % 时间匹配
        S_time = ev.schedule_match_score;

        % 功率近 0 dBm
        S_power_0 = exp(-((ev.power_level_est - sc.prior_power_ref)^2) / ...
                        (sc.sigma_prior_power^2));

        ev.score_prior_pos = sc.beta(1) * S_time + ...
                             sc.beta(2) * S_pos_anchor + ...
                             sc.beta(3) * S_power_0 + ...
                             sc.beta(4) * S_stable;

        %% C. prior_time_known 分数
        ev.score_prior_time = sc.gamma(1) * S_time + ...
                              sc.gamma(2) * S_power_0 + ...
                              sc.gamma(3) * S_stable;

        %% D. ordinary_target 分数
        % 残差项：与前三类最大分数的互补
        max_other = max([ev.score_trusted, ev.score_prior_pos, ev.score_prior_time]);
        S_residual = 1 - max_other;

        % 目标功率范围分数：功率落在 [P_min, P_max] 内
        if ev.power_level_est >= sc.target_power_range(1) - sc.target_power_margin && ...
           ev.power_level_est <= sc.target_power_range(2) + sc.target_power_margin
            S_power_target = 1.0;
        else
            % 软边界
            d_low  = max(0, sc.target_power_range(1) - ev.power_level_est);
            d_high = max(0, ev.power_level_est - sc.target_power_range(2));
            d_out  = max(d_low, d_high);
            S_power_target = exp(-(d_out^2) / (sc.sigma_target_power^2));
        end

        % 持续时间分数
        S_duration = min(1, ev.duration / sc.L_ref);

        ev.score_target = sc.delta(1) * S_residual + ...
                          sc.delta(2) * S_power_target + ...
                          sc.delta(3) * S_duration;

        %% 判定类型
        scores = [ev.score_trusted, ev.score_prior_pos, ...
                  ev.score_prior_time, ev.score_target];
        type_labels = {'trusted_fixed', 'prior_pos_known', ...
                       'prior_time_known', 'ordinary_target'};
        [~, best_idx] = max(scores);
        ev.type_hat = type_labels{best_idx};

        EventList(e) = ev;
    end

    % 打印统计
    if n_events > 0
        types = {EventList.type_hat};
        fprintf('[M3] 事件分类完成: trusted=%d, prior_pos=%d, prior_time=%d, target=%d\n', ...
            sum(strcmp(types, 'trusted_fixed')), ...
            sum(strcmp(types, 'prior_pos_known')), ...
            sum(strcmp(types, 'prior_time_known')), ...
            sum(strcmp(types, 'ordinary_target')));
    end
end


%% ==================== 局部函数 ====================

function sc = fill_score_defaults(Config)
% FILL_SCORE_DEFAULTS  填充评分参数默认值

    if isfield(Config, 'm3') && isfield(Config.m3, 'score')
        s = Config.m3.score;
    else
        s = struct();
    end

    % --- trusted_fixed 评分权重 ---
    if isfield(s, 'alpha')
        sc.alpha = s.alpha;
    else
        sc.alpha = [0.30, 0.25, 0.15, 0.30];  % [cover, power_high, stable, pos_match]
    end

    % trusted 参考功率
    if isfield(s, 'trusted_power_ref')
        sc.trusted_power_ref = s.trusted_power_ref;
    else
        sc.trusted_power_ref = 15;  % dBm
    end
    if isfield(s, 'sigma_trusted_power')
        sc.sigma_trusted_power = s.sigma_trusted_power;
    else
        sc.sigma_trusted_power = 10;  % dB
    end

    % 稳定性高斯核
    if isfield(s, 'kappa_trusted')
        sc.kappa_trusted = s.kappa_trusted;
    else
        sc.kappa_trusted = 5;  % dB^2
    end

    % --- prior_pos_known 评分权重 ---
    if isfield(s, 'beta')
        sc.beta = s.beta;
    else
        sc.beta = [0.35, 0.30, 0.20, 0.15];  % [time, pos, power_0, stable]
    end

    % prior 参考功率
    if isfield(s, 'prior_power_ref')
        sc.prior_power_ref = s.prior_power_ref;
    else
        sc.prior_power_ref = 0;  % dBm
    end
    if isfield(s, 'sigma_prior_power')
        sc.sigma_prior_power = s.sigma_prior_power;
    else
        sc.sigma_prior_power = 8;  % dB
    end

    % --- prior_time_known 评分权重 ---
    if isfield(s, 'gamma')
        sc.gamma = s.gamma;
    else
        sc.gamma = [0.50, 0.30, 0.20];  % [time, power_0, stable]
    end

    % --- ordinary_target 评分权重 ---
    if isfield(s, 'delta')
        sc.delta = s.delta;
    else
        sc.delta = [0.40, 0.35, 0.25];  % [residual, power_range, duration]
    end

    % 目标功率范围
    if isfield(s, 'target_power_range')
        sc.target_power_range = s.target_power_range;
    else
        sc.target_power_range = [-15, -5];  % dBm
    end
    if isfield(s, 'target_power_margin')
        sc.target_power_margin = s.target_power_margin;
    else
        sc.target_power_margin = 3;  % dB
    end
    if isfield(s, 'sigma_target_power')
        sc.sigma_target_power = s.sigma_target_power;
    else
        sc.sigma_target_power = 8;
    end

    % 参考持续帧数
    if isfield(s, 'L_ref')
        sc.L_ref = s.L_ref;
    else
        sc.L_ref = 10;  % 帧
    end
end
