function [M0State, Candidates] = update_persistent_cal_instances_m0( ...
    M0State, SourceTemplates, Config)
% UPDATE_PERSISTENT_CAL_INSTANCES_M0  持续存在的标校源（always-on 保底）
%
%   输出：Candidates 为结构体数组，长度 0 或 1
%     每帧若启用则恒输出一个候选，作为最低优先级保底源。

    Candidates = [];
    if ~Config.m0.source.persistent_enable
        return;
    end

    tpl = SourceTemplates.persistent_cal;

    c = struct();
    c.instance_id      = -1;                         % 持续源使用固定 id
    c.template_key     = tpl.template_key;
    c.source_type      = 'persistent_cal';
    c.source_subtype   = 'persistent_cal';
    c.band_list        = tpl.band_list;
    c.band_id          = tpl.band_id;
    c.true_position    = tpl.fixed_pos_xy;
    c.true_tx_power    = tpl.tx_power_dBm;
    c.tx_power_by_band = [];
    c.is_continuing    = true;
    c.template_idx     = 0;
    c.is_persistent    = true;

    c.location_prior = struct('type','exact','value', tpl.fixed_pos_xy);
    c.power_prior    = struct('type','exact','value', tpl.tx_power_dBm);

    Candidates = c;
end
