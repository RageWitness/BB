function [SourceTemplates, M0State, M0Logs, GridValid, Config, APs] = init_m0_nonoverlap(Config_override)
% INIT_M0_NONOVERLAP  M0 init entry

    Config = default_config_m0();
    if nargin >= 1 && ~isempty(Config_override)
        Config = merge_struct(Config, Config_override);
    end

    APs = generate_aps(Config);
    GridValid = build_valid_grid_nonoverlap(Config, APs);
    SourceTemplates = build_source_templates_nonoverlap(Config, GridValid);
    M0State = init_runtime_state_m0_nonoverlap(Config, SourceTemplates);

    M0Logs.TruthLogTarget = struct( ...
        'frame_id',     {}, ...
        'band_id',      {}, ...
        'instance_id',  {}, ...
        'template_key', {}, ...
        'true_pos_xy',  {}, ...
        'tx_power_dBm', {});

    M0Logs.TruthLogAll = struct( ...
        'frame_id',     {}, ...
        'band_id',      {}, ...
        'source_type',  {}, ...
        'instance_id',  {}, ...
        'template_key', {}, ...
        'true_pos_xy',  {}, ...
        'tx_power_dBm', {});

    fprintf('[M0] ====== M0 initialized ======\n');
end


function Config = default_config_m0()
% DEFAULT_CONFIG_M0  default configuration

    % area
    Config.area.x_range = [0, 300];
    Config.area.y_range = [0, 300];

    % AP layout
    % layout:
    %   'edge_internal' : generate edge + inner APs from counts
    %   'uniform_grid'  : legacy num_x * num_y grid
    %   'explicit'      : use Config.ap.pos_xy directly (N x 2)
    Config.ap.layout = 'edge_internal';
    Config.ap.n_edge = 4;
    Config.ap.n_inner = 1;
    Config.ap.edge_margin_m = 10;
    Config.ap.inner_span_ratio = 0.20;
    Config.ap.num_x = 4;  % legacy fallback
    Config.ap.num_y = 4;  % legacy fallback
    Config.ap.pos_xy = []; % used when layout='explicit'

    % fingerprint grid
    Config.fp.grid_step = 3;
    Config.fp.ap_exclusion_radius = 1.5;

    % M0 core
    Config.m0.num_bands = 4;
    Config.m0.dt = 1;
    Config.m0.T_total = 500;
    Config.m0.priority_order = {'trusted_fixed', 'prior_pos_known', ...
                                'prior_time_known', 'ordinary_target'};

    % trusted
    Config.m0.trusted.TrustedNum = 2;
    Config.m0.trusted.lambda_arrival = 0.005;
    Config.m0.trusted.life_mode = 'geom';
    Config.m0.trusted.life_param = 0.05;

    % prior
    Config.m0.prior.Sprior = [2, 2, 1, 1];
    Config.m0.prior.Tprior = [1, 1, 1, 1];
    Config.m0.prior.schedule_mode = 'periodic';
    Config.m0.prior.period_frames = [50, 60, 70, 80];
    Config.m0.prior.duration_frames = [10, 10, 8, 8];
    Config.m0.prior.phase_frames = [0, 5, 10, 15];

    % ordinary target
    Config.m0.target.Nocoop = [3, 3, 2, 2];
    Config.m0.target.lambda_arrival = [0.02, 0.02, 0.015, 0.015];
    Config.m0.target.life_mode = 'geom';
    Config.m0.target.life_param = 0.08;
    Config.m0.target.life_range = [5, 30];
    Config.m0.target.power_range_dBm = [50, 65; ...
                                        50, 65; ...
                                        50, 65; ...
                                        50, 65];
    Config.m0.target.position_mode = 'uniform';

    % M3 config
    Config.m3.grouping.enable = true;
    Config.m3.grouping.max_gap_frames = 2;
    Config.m3.relative_power.enable = true;
    Config.m3.relative_power.topk = 3;
    Config.m3.relative_power.tau_valid_ap_dB = 3;

    Config.m3.trusted_hard.enable = true;
    Config.m3.trusted_hard.power_threshold_lin = 1e5;
    Config.m3.trusted_hard.require_all_bands = true;
    Config.m3.trusted_hard.require_all_aps = true;
    Config.m3.trusted_hard.min_duration_frames = 2;
    Config.m3.trusted_hard.aggregation = 'mean';
    Config.m3.trusted_hard.mode = 'strict';
    Config.m3.trusted_hard.min_fraction = 0.95;

    Config.m3.prior_pos.min_position_score = 0.60;
    Config.m3.prior_pos.min_margin = 0.10;
    Config.m3.prior_time.min_time_score = 0.60;
    Config.m3.prior_time.min_margin = 0.10;

    Config.m3.ordinary.min_ratio_sum = 3;
    Config.m3.ordinary.min_ratio_top1 = 5;
    Config.m3.ordinary.min_valid_ap = 2;
    Config.m3.ordinary.min_duration_frames = 3;

    Config.m3.linear_power.topk = 3;
    Config.m3.linear_power.alpha_noise = 2;
    Config.m3.route.hold_enable = true;

    % M4 config
    Config.m4.distance_mode = 'shape_scale_masked';
    Config.m4.masked_shape.enable = true;
    Config.m4.masked_shape.tau_low_dB = 3;
    Config.m4.masked_shape.tau_high_dB = 10;
    Config.m4.masked_shape.eps = 1e-12;
    Config.m4.lambda_shape = 0.7;
    Config.m4.lambda_resid = 0.3;
end


function APs = generate_aps(Config)
% GENERATE_APS  build AP positions

    layout = get_struct_field(Config.ap, 'layout', 'uniform_grid');
    x_range = Config.area.x_range;
    y_range = Config.area.y_range;

    if strcmp(layout, 'explicit')
        pos_xy = get_struct_field(Config.ap, 'pos_xy', []);
        if isempty(pos_xy) || size(pos_xy, 2) ~= 2
            error('[M0] explicit layout requires Config.ap.pos_xy as N x 2');
        end
        APs.pos_xy = pos_xy;
        APs.num = size(pos_xy, 1);
        APs.role = repmat({'explicit'}, APs.num, 1);
        fprintf('[M0] AP count: %d (layout=explicit)\n', APs.num);
        return;
    end

    if strcmp(layout, 'edge_internal') || strcmp(layout, 'edge_internal_5')
        n_edge = get_struct_field(Config.ap, 'n_edge', 2);
        n_inner = get_struct_field(Config.ap, 'n_inner', 3);
        if (n_edge < 1) || (n_inner < 0)
            error('[M0] edge_internal requires n_edge>=1 and n_inner>=0');
        end

        edge_margin = get_struct_field(Config.ap, 'edge_margin_m', 10);
        span_ratio = get_struct_field(Config.ap, 'inner_span_ratio', 0.20);

        edge_pos = generate_edge_positions(x_range, y_range, edge_margin, n_edge);
        inner_pos = generate_inner_positions(x_range, y_range, span_ratio, n_inner);

        APs.pos_xy = [edge_pos; inner_pos];
        APs.num = size(APs.pos_xy, 1);
        APs.role = [repmat({'edge'}, n_edge, 1); repmat({'inner'}, n_inner, 1)];
        fprintf('[M0] AP count: %d (layout=edge_internal, edge=%d, inner=%d)\n', ...
            APs.num, n_edge, n_inner);
        return;
    end

    nx = Config.ap.num_x;
    ny = Config.ap.num_y;
    margin_x = (x_range(2) - x_range(1)) / (2 * nx);
    margin_y = (y_range(2) - y_range(1)) / (2 * ny);
    ax = linspace(x_range(1) + margin_x, x_range(2) - margin_x, nx);
    ay = linspace(y_range(1) + margin_y, y_range(2) - margin_y, ny);
    [AX, AY] = meshgrid(ax, ay);

    APs.pos_xy = [AX(:), AY(:)];
    APs.num = nx * ny;
    APs.role = repmat({'grid'}, APs.num, 1);
    fprintf('[M0] AP count: %d (layout=uniform_grid, %dx%d)\n', APs.num, nx, ny);
end


function edge_pos = generate_edge_positions(x_range, y_range, edge_margin, n_edge)
% GENERATE_EDGE_POSITIONS  edge AP positions

    xmin = x_range(1) + edge_margin;
    xmax = x_range(2) - edge_margin;
    ymin = y_range(1) + edge_margin;
    ymax = y_range(2) - edge_margin;
    xc = 0.5 * (x_range(1) + x_range(2));
    yc = 0.5 * (y_range(1) + y_range(2));

    if xmin >= xmax || ymin >= ymax
        error('[M0] edge_margin_m is too large for current area');
    end

    if n_edge == 4
        edge_pos = [xc, ymax; ...
                    xmax, yc; ...
                    xc, ymin; ...
                    xmin, yc];
        return;
    end

    base = floor(n_edge / 4);
    remn = mod(n_edge, 4);
    cnt = base * ones(1, 4);
    cnt(1:remn) = cnt(1:remn) + 1;

    edge_pos = zeros(n_edge, 2);
    k = 1;

    xs_top = side_samples(xmin, xmax, cnt(1));
    for i = 1:numel(xs_top)
        edge_pos(k, :) = [xs_top(i), ymax]; k = k + 1;
    end

    ys_right = side_samples(ymax, ymin, cnt(2));
    for i = 1:numel(ys_right)
        edge_pos(k, :) = [xmax, ys_right(i)]; k = k + 1;
    end

    xs_bottom = side_samples(xmax, xmin, cnt(3));
    for i = 1:numel(xs_bottom)
        edge_pos(k, :) = [xs_bottom(i), ymin]; k = k + 1;
    end

    ys_left = side_samples(ymin, ymax, cnt(4));
    for i = 1:numel(ys_left)
        edge_pos(k, :) = [xmin, ys_left(i)]; k = k + 1;
    end
end


function inner_pos = generate_inner_positions(x_range, y_range, span_ratio, n_inner)
% GENERATE_INNER_POSITIONS  inner AP positions

    if n_inner <= 0
        inner_pos = zeros(0, 2);
        return;
    end

    xc = 0.5 * (x_range(1) + x_range(2));
    yc = 0.5 * (y_range(1) + y_range(2));

    if n_inner == 1
        inner_pos = [xc, yc];
        return;
    end

    dx = span_ratio * (x_range(2) - x_range(1));
    dy = span_ratio * (y_range(2) - y_range(1));
    n_ring = n_inner - 1;
    th = linspace(0, 2 * pi, n_ring + 1);
    th(end) = [];
    ring = [xc + dx * cos(th(:)), yc + dy * sin(th(:))];
    inner_pos = [xc, yc; ring];
end


function vals = side_samples(a, b, n)
% SIDE_SAMPLES  sample n points on segment [a,b]
    if n <= 0
        vals = [];
    else
        vals = linspace(a, b, n + 2);
        vals = vals(2:end-1);
    end
end


function S = merge_struct(S, T)
% MERGE_STRUCT  recursive struct merge
    fnames = fieldnames(T);
    for i = 1:numel(fnames)
        fn = fnames{i};
        if isfield(S, fn) && isstruct(S.(fn)) && isstruct(T.(fn))
            S.(fn) = merge_struct(S.(fn), T.(fn));
        else
            S.(fn) = T.(fn);
        end
    end
end


function val = get_struct_field(S, field_name, default_val)
% GET_STRUCT_FIELD  read field with default
    if isfield(S, field_name)
        val = S.(field_name);
    else
        val = default_val;
    end
end
