function GridValid = build_valid_grid_nonoverlap(Config, APs)
% BUILD_VALID_GRID_NONOVERLAP  生成避开 AP 重合位置的有效网格集合
%
%   GridValid = build_valid_grid_nonoverlap(Config, APs)
%
%   有效网格点满足：min_i ||r_g - a_i||_2 > r_excl
%
%   输入：
%       Config  - 配置结构体，需包含：
%                   Config.area.x_range  = [xmin, xmax]  (m)
%                   Config.area.y_range  = [ymin, ymax]  (m)
%                   Config.fp.grid_step  (m)
%                   Config.fp.ap_exclusion_radius  (m)
%       APs     - AP 结构体，需包含：
%                   APs.pos_xy  (Nap x 2)
%
%   输出：
%       GridValid.xy       - (Nvalid x 2) 有效网格坐标
%       GridValid.mask     - (Ny x Nx) logical，标记有效位置
%       GridValid.x_vec    - 网格 x 轴坐标向量
%       GridValid.y_vec    - 网格 y 轴坐标向量
%       GridValid.Nvalid   - 有效网格点数

    step = Config.fp.grid_step;
    r_excl = Config.fp.ap_exclusion_radius;

    x_vec = Config.area.x_range(1) : step : Config.area.x_range(2);
    y_vec = Config.area.y_range(1) : step : Config.area.y_range(2);

    [X, Y] = meshgrid(x_vec, y_vec);
    Nx = length(x_vec);
    Ny = length(y_vec);

    % 计算每个网格点到所有 AP 的最小距离
    all_xy = [X(:), Y(:)];  % (Ny*Nx) x 2
    Nap = size(APs.pos_xy, 1);

    min_dist = inf(size(all_xy, 1), 1);
    for i = 1:Nap
        d = sqrt((all_xy(:,1) - APs.pos_xy(i,1)).^2 + ...
                 (all_xy(:,2) - APs.pos_xy(i,2)).^2);
        min_dist = min(min_dist, d);
    end

    valid_flag = min_dist > r_excl;  % (Ny*Nx) x 1

    GridValid.x_vec = x_vec;
    GridValid.y_vec = y_vec;
    GridValid.mask  = reshape(valid_flag, Ny, Nx);
    GridValid.xy    = all_xy(valid_flag, :);
    GridValid.Nvalid = size(GridValid.xy, 1);

    fprintf('[M0] 有效网格点: %d / %d (排除半径 %.1f m)\n', ...
        GridValid.Nvalid, Ny*Nx, r_excl);
end
