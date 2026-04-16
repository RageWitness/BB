function GroupList = route_group_m3(GroupList, Config)
% ROUTE_GROUP_M3  group 级路由
%
%   trusted_fixed    -> calibrate_direct
%   prior_pos_known  -> calibrate_direct
%   prior_time_known -> localize_then_calibrate
%   ordinary_target  -> localize_only
%   unknown          -> hold

    if isempty(GroupList)
        return;
    end

    GroupList = ensure_group_fields_route(GroupList);
    hold_enable = true;
    if isfield(Config, 'm3') && isfield(Config.m3, 'route') && ...
       isfield(Config.m3.route, 'hold_enable')
        hold_enable = logical(Config.m3.route.hold_enable);
    end

    for g = 1:numel(GroupList)
        grp = GroupList(g);
        switch grp.type_hat_group
            case 'trusted_fixed'
                grp.route_action_group = 'calibrate_direct';
            case 'prior_pos_known'
                grp.route_action_group = 'calibrate_direct';
            case 'prior_time_known'
                grp.route_action_group = 'localize_then_calibrate';
            case 'ordinary_target'
                grp.route_action_group = 'localize_only';
            otherwise
                if hold_enable
                    grp.route_action_group = 'hold';
                else
                    grp.route_action_group = 'localize_only';
                end
                if isempty(grp.hold_reason)
                    grp.hold_reason = 'route_hold';
                end
        end
        GroupList(g) = grp;
    end

    if ~isempty(GroupList)
        types = {GroupList.type_hat_group};
        fprintf('[M3] group 分类: trusted=%d, prior_pos=%d, prior_time=%d, ordinary=%d, unknown=%d\n', ...
            sum(strcmp(types, 'trusted_fixed')), ...
            sum(strcmp(types, 'prior_pos_known')), ...
            sum(strcmp(types, 'prior_time_known')), ...
            sum(strcmp(types, 'ordinary_target')), ...
            sum(strcmp(types, 'unknown')));
    end
end


function GroupList = ensure_group_fields_route(GroupList)
% ENSURE_GROUP_FIELDS_ROUTE  统一 route 字段，避免异构赋值
    if ~isfield(GroupList, 'route_action_group')
        [GroupList.route_action_group] = deal('hold');
    end
    if ~isfield(GroupList, 'hold_reason')
        [GroupList.hold_reason] = deal('');
    end
    if ~isfield(GroupList, 'type_hat_group')
        [GroupList.type_hat_group] = deal('unknown');
    end
end
