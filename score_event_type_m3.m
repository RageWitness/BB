function EventList = score_event_type_m3(EventList, Config)
% SCORE_EVENT_TYPE_M3  对每个事件计算四类源评分并判定类型
%
%   当 Config.m3.relative_power.enable = true 时，优先使用
%   相对噪声底功率特征评分；关闭时回退到旧的绝对功率评分。

    sc = fill_score_defaults(Config);
    n_events = numel(EventList);

    % 预分配新字段
    for e = 1:n_events
        EventList(e).score_trusted    = 0;
        EventList(e).score_prior_pos  = 0;
        EventList(e).score_prior_time = 0;
        EventList(e).score_target     = 0;
        EventList(e).type_hat         = '';
    end

    for e = 1:n_events
        ev = EventList(e);

        if sc.use_relative_power
            ev = score_relative_mode(ev, sc);
        else
            ev = score_absolute_mode(ev, sc);
        end

        scores = [ev.score_trusted, ev.score_prior_pos, ...
                  ev.score_prior_time, ev.score_target];
        type_labels = {'trusted_fixed', 'prior_pos_known', ...
                       'prior_time_known', 'ordinary_target'};
        [~, best_idx] = max(scores);
        ev.type_hat = type_labels{best_idx};

        EventList(e) = ev;
    end

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

function ev = score_relative_mode(ev, sc)
% SCORE_RELATIVE_MODE  基于相对噪声底特征的评分

    % 公共项
    S_cover = sum(ev.band_coverage_vec) / numel(ev.band_coverage_vec);
    S_time  = ev.schedule_match_score;
    S_pos   = ev.position_prior_match_score;
    S_stable = exp(-(ev.power_stability_est) / sc.kappa_trusted);

    ex_mean = get_field_default(ev, 'power_excess_mean_dB', ev.power_level_est);
    ex_top1 = get_field_default(ev, 'power_excess_top1_dB', ev.power_level_est);
    ex_topk = get_field_default(ev, 'power_excess_topk_dB', ev.power_level_est);
    n_valid = get_field_default(ev, 'n_valid_ap', 0);

    % A. trusted
    S_sig_trusted = 0.5 * sigmoid((ex_mean - sc.trusted_excess_mean_ref) / sc.trusted_excess_scale) + ...
                    0.5 * sigmoid((ex_top1 - sc.trusted_excess_top1_ref) / sc.trusted_excess_scale);
    S_valid_trusted = min(1, n_valid / max(sc.trusted_n_valid_ref, 1));
    ev.score_trusted = sc.alpha_rel(1) * S_cover + ...
                       sc.alpha_rel(2) * S_sig_trusted + ...
                       sc.alpha_rel(3) * S_valid_trusted + ...
                       sc.alpha_rel(4) * S_stable + ...
                       sc.alpha_rel(5) * S_pos;

    % B. prior_pos
    S_sig_prior = exp(-((ex_topk - sc.prior_excess_ref)^2) / (sc.prior_excess_sigma^2));
    S_valid_prior = min(1, n_valid / max(sc.prior_n_valid_ref, 1));
    ev.score_prior_pos = sc.beta_rel(1) * S_time + ...
                         sc.beta_rel(2) * S_pos + ...
                         sc.beta_rel(3) * S_sig_prior + ...
                         sc.beta_rel(4) * S_stable + ...
                         sc.beta_rel(5) * S_valid_prior;

    % C. prior_time
    ev.score_prior_time = sc.gamma_rel(1) * S_time + ...
                          sc.gamma_rel(2) * S_sig_prior + ...
                          sc.gamma_rel(3) * S_stable + ...
                          sc.gamma_rel(4) * S_valid_prior;

    % D. target
    max_other = max([ev.score_trusted, ev.score_prior_pos, ev.score_prior_time]);
    S_residual = max(0, 1 - max_other);
    S_sig_target = sigmoid((ex_top1 - sc.target_excess_top1_ref) / sc.target_excess_scale);
    S_valid_target = exp(-((n_valid - sc.target_n_valid_ref)^2) / ...
                        (2 * max(sc.target_n_valid_sigma, 1e-6)^2));
    S_duration = min(1, ev.duration / sc.L_ref);
    ev.score_target = sc.delta_rel(1) * S_residual + ...
                      sc.delta_rel(2) * S_sig_target + ...
                      sc.delta_rel(3) * S_valid_target + ...
                      sc.delta_rel(4) * S_duration;
end


function ev = score_absolute_mode(ev, sc)
% SCORE_ABSOLUTE_MODE  旧版绝对功率评分（兼容保留）

    % A. trusted
    S_cover = sum(ev.band_coverage_vec) / numel(ev.band_coverage_vec);
    S_power_high = exp(-((ev.power_level_est - sc.trusted_power_ref)^2) / ...
                       (sc.sigma_trusted_power^2));
    S_stable = exp(-(ev.power_stability_est) / sc.kappa_trusted);
    S_pos_anchor = ev.position_prior_match_score;
    ev.score_trusted = sc.alpha(1) * S_cover + ...
                       sc.alpha(2) * S_power_high + ...
                       sc.alpha(3) * S_stable + ...
                       sc.alpha(4) * S_pos_anchor;

    % B. prior_pos
    S_time = ev.schedule_match_score;
    S_power_0 = exp(-((ev.power_level_est - sc.prior_power_ref)^2) / ...
                    (sc.sigma_prior_power^2));
    ev.score_prior_pos = sc.beta(1) * S_time + ...
                         sc.beta(2) * S_pos_anchor + ...
                         sc.beta(3) * S_power_0 + ...
                         sc.beta(4) * S_stable;

    % C. prior_time
    ev.score_prior_time = sc.gamma(1) * S_time + ...
                          sc.gamma(2) * S_power_0 + ...
                          sc.gamma(3) * S_stable;

    % D. target
    max_other = max([ev.score_trusted, ev.score_prior_pos, ev.score_prior_time]);
    S_residual = 1 - max_other;
    if ev.power_level_est >= sc.target_power_range(1) - sc.target_power_margin && ...
       ev.power_level_est <= sc.target_power_range(2) + sc.target_power_margin
        S_power_target = 1.0;
    else
        d_low  = max(0, sc.target_power_range(1) - ev.power_level_est);
        d_high = max(0, ev.power_level_est - sc.target_power_range(2));
        d_out  = max(d_low, d_high);
        S_power_target = exp(-(d_out^2) / (sc.sigma_target_power^2));
    end
    S_duration = min(1, ev.duration / sc.L_ref);
    ev.score_target = sc.delta(1) * S_residual + ...
                      sc.delta(2) * S_power_target + ...
                      sc.delta(3) * S_duration;
end


function sc = fill_score_defaults(Config)
% FILL_SCORE_DEFAULTS  填充评分参数默认值

    if isfield(Config, 'm3') && isfield(Config.m3, 'score')
        s = Config.m3.score;
    else
        s = struct();
    end

    if isfield(Config, 'm3') && isfield(Config.m3, 'relative_power')
        rp = Config.m3.relative_power;
    else
        rp = struct();
    end
    if isfield(rp, 'enable')
        sc.use_relative_power = rp.enable;
    else
        sc.use_relative_power = true;
    end

    % ---- relative 模式参数 ----
    sc.alpha_rel = get_field_default(s, 'alpha_rel', [0.20, 0.30, 0.20, 0.10, 0.20]);
    sc.beta_rel  = get_field_default(s, 'beta_rel',  [0.35, 0.30, 0.15, 0.10, 0.10]);
    sc.gamma_rel = get_field_default(s, 'gamma_rel', [0.45, 0.25, 0.15, 0.15]);
    sc.delta_rel = get_field_default(s, 'delta_rel', [0.45, 0.25, 0.15, 0.15]);

    sc.trusted_excess_mean_ref = get_field_default(s, 'trusted_excess_mean_ref', 10);
    sc.trusted_excess_top1_ref = get_field_default(s, 'trusted_excess_top1_ref', 15);
    sc.trusted_excess_scale    = get_field_default(s, 'trusted_excess_scale', 4);
    sc.trusted_n_valid_ref     = get_field_default(s, 'trusted_n_valid_ref', 6);

    sc.prior_excess_ref   = get_field_default(s, 'prior_excess_ref', 8);
    sc.prior_excess_sigma = get_field_default(s, 'prior_excess_sigma', 6);
    sc.prior_n_valid_ref  = get_field_default(s, 'prior_n_valid_ref', 4);

    sc.target_excess_top1_ref = get_field_default(s, 'target_excess_top1_ref', 5);
    sc.target_excess_scale    = get_field_default(s, 'target_excess_scale', 3);
    sc.target_n_valid_ref     = get_field_default(s, 'target_n_valid_ref', 3);
    sc.target_n_valid_sigma   = get_field_default(s, 'target_n_valid_sigma', 2);

    % relative/absolute 共用参数
    sc.kappa_trusted = get_field_default(s, 'kappa_trusted', 5);
    sc.L_ref         = get_field_default(s, 'L_ref', 10);

    % ---- absolute 旧逻辑参数（兼容保留）----
    sc.alpha = get_field_default(s, 'alpha', [0.30, 0.25, 0.15, 0.30]);
    sc.trusted_power_ref = get_field_default(s, 'trusted_power_ref', 15);
    sc.sigma_trusted_power = get_field_default(s, 'sigma_trusted_power', 10);

    sc.beta = get_field_default(s, 'beta', [0.35, 0.30, 0.20, 0.15]);
    sc.prior_power_ref = get_field_default(s, 'prior_power_ref', 0);
    sc.sigma_prior_power = get_field_default(s, 'sigma_prior_power', 8);

    sc.gamma = get_field_default(s, 'gamma', [0.50, 0.30, 0.20]);
    sc.delta = get_field_default(s, 'delta', [0.40, 0.35, 0.25]);
    sc.target_power_range = get_field_default(s, 'target_power_range', [50, 65]);
    sc.target_power_margin = get_field_default(s, 'target_power_margin', 3);
    sc.sigma_target_power = get_field_default(s, 'sigma_target_power', 8);
end


function y = sigmoid(x)
    y = 1 ./ (1 + exp(-x));
end


function v = get_field_default(s, fn, def)
    if isstruct(s) && isfield(s, fn)
        v = s.(fn);
    else
        v = def;
    end
end
