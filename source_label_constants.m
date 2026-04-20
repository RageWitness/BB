function C = source_label_constants()
% SOURCE_LABEL_CONSTANTS  统一标签常量表（全局唯一定义）
%
%   C = source_label_constants()
%
%   标签值：
%     1 = persistent_cal   持续存在的标校源
%     2 = broadband_cal    宽带标校源
%     3 = target           待定位源
%     4 = opportunistic    机遇源 / 先验源

    C.PERSISTENT_CAL = 1;
    C.BROADBAND_CAL  = 2;
    C.TARGET         = 3;
    C.OPPORTUNISTIC  = 4;

    C.ALL_LABELS = [1, 2, 3, 4];

    C.name_map = containers.Map( ...
        {1, 2, 3, 4}, ...
        {'persistent_cal', 'broadband_cal', 'target', 'opportunistic'});

    C.route_map = containers.Map( ...
        {1, 2, 3, 4}, ...
        {'calibrate_direct', 'calibrate_direct', 'localize_only', 'localize_then_calibrate'});

    C.type_to_label = containers.Map( ...
        {'persistent_cal', 'broadband_cal', 'target', 'opportunistic'}, ...
        {1, 2, 3, 4});
end
