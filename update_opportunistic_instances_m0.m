function [M0State, OppCandidatesPerBand] = update_opportunistic_instances_m0( ...
    M0State, SourceTemplates, GridValid, Config)
% UPDATE_OPPORTUNISTIC_INSTANCES_M0  机遇源到达 + 三类位置先验 + 精确功率先验
%
%   每帧每频带至多一个机遇源。Poisson 到达；到达时同时生成真实参数与先验信息，
%   并保证两者一致。

    B  = Config.m0.num_bands;
    dt = Config.m0.dt;
    OppCfg = Config.m0.source.opportunistic;
    tpl    = SourceTemplates.opportunistic;

    OppCandidatesPerBand = struct('candidates', cell(1, B));

    for b = 1:B
        was_active = M0State.opp_active_band(b);

        if was_active
            M0State.opp_life_band(b) = M0State.opp_life_band(b) - 1;
            if M0State.opp_life_band(b) <= 0
                M0State.opp_active_band(b) = false;
                M0State.opp_inst_band(b) = 0;
            end
        end

        if M0State.opp_active_band(b)
            inst = find_instance(M0State, M0State.opp_inst_band(b));
            if ~isempty(inst)
                OppCandidatesPerBand(b).candidates = build_opp_candidate(inst, b, true);
            end
        else
            % 频率：opportunistic 是"全局"到达率，按频带均分（或按用户的 'uniform' 策略）
            lambda_b = tpl.lambda_arrival / B;
            p = 1 - exp(-lambda_b * dt);
            if rand() < p
                % 1) 决定位置先验类型并联动生成真实位置
                [loc_prior, true_pos] = sample_location_prior_and_truth( ...
                    OppCfg, GridValid, Config.m0.buildings);

                % 2) 决定功率先验
                power_sample = sample_power_in_range(OppCfg.power_range_dBm);
                if rand() < OppCfg.power_prior_prob.exact
                    power_prior   = struct('type','exact','value', power_sample);
                    true_tx_power = power_sample;
                else
                    power_prior   = struct('type','none','value', []);
                    true_tx_power = power_sample;
                end

                life = sample_life(tpl.life_mode, tpl.life_param, OppCfg.life_range, Config);

                inst = struct();
                inst.instance_id      = M0State.next_instance_id;
                inst.template_key     = tpl.template_key;
                inst.source_type      = 'opportunistic';
                inst.band_list        = b;
                inst.true_position    = true_pos;
                inst.true_tx_power    = true_tx_power;
                inst.tx_power_by_band = [];
                inst.location_prior   = loc_prior;
                inst.power_prior      = power_prior;
                inst.life_remaining   = life;
                inst.is_active        = true;

                M0State.next_instance_id = M0State.next_instance_id + 1;
                M0State.instances = append_instance(M0State.instances, inst);
                M0State.opp_active_band(b) = true;
                M0State.opp_life_band(b)   = life;
                M0State.opp_inst_band(b)   = inst.instance_id;

                OppCandidatesPerBand(b).candidates = build_opp_candidate(inst, b, false);
            end
        end
    end
end


function cand = build_opp_candidate(inst, b, is_continuing)
    cand.instance_id      = inst.instance_id;
    cand.template_key     = inst.template_key;
    cand.source_type      = 'opportunistic';
    cand.source_subtype   = 'opportunistic';
    cand.band_list        = b;
    cand.band_id          = b;
    cand.true_position    = inst.true_position;
    cand.true_tx_power    = inst.true_tx_power;
    cand.tx_power_by_band = [];
    cand.is_continuing    = is_continuing;
    cand.template_idx     = 1;
    cand.is_persistent    = false;
    cand.location_prior   = inst.location_prior;
    cand.power_prior      = inst.power_prior;
end


function [loc_prior, true_pos] = sample_location_prior_and_truth(OppCfg, GridValid, Buildings)
% 三选一：exact / region / gaussian
    pE = OppCfg.location_prior_prob.exact;
    pR = OppCfg.location_prior_prob.region;
    pG = OppCfg.location_prior_prob.gaussian;
    psum = pE + pR + pG;
    pE = pE / psum; pR = pR / psum;

    u = rand();
    if u < pE
        % --- exact ---
        idx = randi(GridValid.Nvalid);
        true_pos = GridValid.xy(idx, :);
        loc_prior = struct('type','exact','value', true_pos);
    elseif u < pE + pR
        % --- region: 选一栋建筑，bbox 内采样 ---
        N_bld = numel(Buildings);
        k = randi(N_bld);
        bld = Buildings(k);
        bbox = [bld.xmin, bld.xmax, bld.ymin, bld.ymax];
        true_pos = sample_in_bbox(bbox);
        loc_prior = struct('type','region', ...
            'value', struct('building_id', k, 'bbox', bbox));
    else
        % --- gaussian: 中心从合法网格挑，sigma 各向同性 ---
        idx = randi(GridValid.Nvalid);
        mu  = GridValid.xy(idx, :);
        sig = OppCfg.gaussian.sigma_default;
        true_pos = sample_gaussian_with_clip(mu, sig, GridValid);
        loc_prior = struct('type','gaussian', ...
            'value', struct('mu', mu, 'sigma', sig));
    end
end


function p = sample_in_bbox(bbox)
    p = [bbox(1) + (bbox(2)-bbox(1))*rand(), ...
         bbox(3) + (bbox(4)-bbox(3))*rand()];
end


function pos = sample_gaussian_with_clip(mu, sigma, GridValid)
    for tryK = 1:20
        cand = mu + sigma * randn(1, 2);
        % 在合法网格中找最近点检查距离
        d2 = sum((GridValid.xy - cand).^2, 2);
        [dmin, imin] = min(d2);
        if dmin < (1.5 * Config_grid_step(GridValid))^2
            pos = cand;
            return;
        end
    end
    % 回退：取最近合法网格点
    d2 = sum((GridValid.xy - mu).^2, 2);
    [~, imin] = min(d2);
    pos = GridValid.xy(imin, :);
end


function s = Config_grid_step(GridValid)
    if isfield(GridValid, 'step') && ~isempty(GridValid.step)
        s = GridValid.step;
    else
        s = 3;  % 与 Config.fp.grid_step 默认值对齐
    end
end


function p = sample_power_in_range(range_dBm)
    p = range_dBm(1) + (range_dBm(2) - range_dBm(1)) * rand();
end


function life = sample_life(mode, param, life_range, Config)
    switch mode
        case 'geom'
            life = geornd(param) + 1;
        case 'uniform'
            life = randi([life_range(1), life_range(2)]);
        case 'exp'
            life = max(1, ceil(exprnd(param) / Config.m0.dt));
        otherwise
            life = 10;
    end
end


function inst = find_instance(M0State, inst_id)
    inst = [];
    for k = 1:numel(M0State.instances)
        if M0State.instances(k).instance_id == inst_id
            inst = M0State.instances(k);
            return;
        end
    end
end


function arr = append_instance(arr, inst)
    if isempty(arr)
        arr = inst;
    else
        arr(end+1) = inst;
    end
end
