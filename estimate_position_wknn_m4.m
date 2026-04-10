function est_pos_xy = estimate_position_wknn_m4(grid_xy, neighbor_idx, weights)
% ESTIMATE_POSITION_WKNN_M4  加权 K 近邻位置估计
%
%   est_pos_xy = estimate_position_wknn_m4(grid_xy, neighbor_idx, weights)
%
%   公式：x_hat = sum(w_g * r_g) / sum(w_g)
%   其中 weights 已归一化，所以直接加权求和。
%
%   输入：
%       grid_xy      - (G x 2) 参考点坐标
%       neighbor_idx - (1 x K) 近邻索引
%       weights      - (1 x K) 归一化权重
%
%   输出：
%       est_pos_xy   - (1 x 2) 估计位置

    % 取近邻坐标 (K x 2)
    nn_pos = grid_xy(neighbor_idx, :);

    % 加权求和
    est_pos_xy = weights * nn_pos;  % (1 x K) * (K x 2) = (1 x 2)
end
