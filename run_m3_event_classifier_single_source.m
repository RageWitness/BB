function [EventList, GroupList] = run_m3_event_classifier_single_source( ...
    Y_dBm_all, Y_lin_all, SpatialFP, SignatureLib, SourceTemplates, Config)
% RUN_M3_EVENT_CLASSIFIER_SINGLE_SOURCE
% Simplified M3 path:
% 1) detect per-band events
% 2) extract event features
% 3) assign unified non-typed localization routing
%
% Source-type judgment modules are disconnected.

    %#ok<INUSD>
    fprintf('\n============================================\n');
    fprintf('  M3 Event Detection (type-classification removed)\n');
    fprintf('============================================\n\n');

    tic;

    EventListRaw = detect_band_events_m3(Y_dBm_all, Y_lin_all, Config);
    if isempty(EventListRaw)
        fprintf('[M3] No events detected\n');
        EventList = [];
        GroupList = [];
        return;
    end

    EventList = extract_event_features_m3( ...
        EventListRaw, Y_dBm_all, Y_lin_all, SpatialFP, SignatureLib, SourceTemplates, Config);

    EventList = assign_untyped_routing_m3(EventList);
    GroupList = [];

    elapsed = toc;
    fprintf('\n[M3] Done: %d events, %.2f s (no source-type classification)\n', ...
        numel(EventList), elapsed);
end


%% ==================== local helpers ====================

function EventList = assign_untyped_routing_m3(EventList)
% All detected events are routed to localization directly.

    for e = 1:numel(EventList)
        EventList(e).type_hat = 'untyped_source';
        EventList(e).route_action = 'localize_only';
        EventList(e).linked_template_key = '';
        EventList(e).upgrade_hint = 'none';

        % Keep score fields for compatibility with existing logs/plots.
        EventList(e).score_trusted = 0;
        EventList(e).score_prior_pos = 0;
        EventList(e).score_prior_time = 0;
        EventList(e).score_target = 1;
    end
end

