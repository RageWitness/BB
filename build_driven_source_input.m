function dsi = build_driven_source_input(apb, Config)
% BUILD_DRIVEN_SOURCE_INPUT  将 SourceEvent 封装为统一标签驱动输入结构
%
%   dsi = build_driven_source_input(apb, Config)
%
%   输入：
%     apb    - FrameState_t.active_per_band(b)，来自 finalize_frame_state_m0
%     Config - 全局配置（含 Config.debug.expose_true_source_state）
%
%   输出：
%     dsi - DrivenSourceInput 结构，供定位/标校链路使用

    LC = source_label_constants();

    if ~apb.has_source
        dsi = empty_dsi(apb.band_id, apb.timestamp);
        return;
    end

    dsi.label     = LC.type_to_label(apb.source_type);
    dsi.band_id   = apb.band_id;
    dsi.band_list = apb.band_list;
    dsi.timestamp = apb.timestamp;

    dsi.location_prior = pack_location_prior(apb);
    dsi.power_prior    = pack_power_prior(apb);

    dsi.meta.source_type_name      = apb.source_type;
    dsi.meta.is_calibration_source = strcmp(apb.source_type, 'broadband_cal') || ...
                                     strcmp(apb.source_type, 'persistent_cal');
    dsi.meta.is_target_source      = strcmp(apb.source_type, 'target');
    dsi.meta.is_opportunistic_source = strcmp(apb.source_type, 'opportunistic');

    expose = isfield(Config, 'debug') && ...
             isfield(Config.debug, 'expose_true_source_state') && ...
             Config.debug.expose_true_source_state;
    if expose
        dsi.debug_info.true_position    = apb.true_position;
        dsi.debug_info.true_tx_power    = apb.true_tx_power;
        dsi.debug_info.tx_power_by_band = apb.tx_power_by_band;
        dsi.debug_info.instance_id      = apb.instance_id;
        dsi.debug_info.template_key     = apb.template_key;
    else
        dsi.debug_info.true_position    = [];
        dsi.debug_info.true_tx_power    = [];
        dsi.debug_info.tx_power_by_band = [];
        dsi.debug_info.instance_id      = [];
        dsi.debug_info.template_key     = '';
    end
end


function lp = pack_location_prior(apb)
    lp = apb.location_prior;
    if strcmp(lp.type, 'region') && isstruct(lp.value)
        if isfield(lp.value, 'building_id') && ~isfield(lp.value, 'region_id')
            lp.value.region_id   = lp.value.building_id;
            lp.value.region_name = sprintf('building_%d', lp.value.building_id);
        end
        if isfield(lp.value, 'bbox') && ~isfield(lp.value, 'region_geometry')
            lp.value.region_geometry = lp.value.bbox;
        end
    end
end


function pp = pack_power_prior(apb)
    pp = apb.power_prior;
    if strcmp(apb.source_type, 'broadband_cal') && ~isempty(apb.tx_power_by_band)
        pp.type  = 'exact_by_band';
        pp.value = apb.tx_power_by_band;
    end
end


function dsi = empty_dsi(band_id, timestamp)
    dsi.label     = 0;
    dsi.band_id   = band_id;
    dsi.band_list = band_id;
    dsi.timestamp = timestamp;

    dsi.location_prior = struct('type', 'none', 'value', []);
    dsi.power_prior    = struct('type', 'none', 'value', []);

    dsi.meta.source_type_name        = '';
    dsi.meta.is_calibration_source   = false;
    dsi.meta.is_target_source        = false;
    dsi.meta.is_opportunistic_source = false;

    dsi.debug_info.true_position    = [];
    dsi.debug_info.true_tx_power    = [];
    dsi.debug_info.tx_power_by_band = [];
    dsi.debug_info.instance_id      = [];
    dsi.debug_info.template_key     = '';
end
