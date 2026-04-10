function [neighbor_idx, neighbor_dist, weights] = select_knn_neighbors_m4( ...
    dist_vec, K, eps_val)
% SELECT_KNN_NEIGHBORS_M4  从距离向量中选取 K 个最近邻并计算权重
%
%   [neighbor_idx, neighbor_dist, weights] = select_knn_neighbors_m4(
%       dist_vec, K, eps_val)
%
%   算法：
%     1. 对 dist_vec 升序排序，取前 K 个
%     2. 权重：w_g = 1 / (d_g + eps)
%     3. 归一化权重
%
%   输入：
%       dist_vec  - (1 x G) 与每个参考点的距离
%       K         - 近邻数
%       eps_val   - 防除零常数（默认 1e-6）
%
%   输出：
%       neighbor_idx  - (1 x K) 最近邻索引
%       neighbor_dist - (1 x K) 最近邻距离
%       weights       - (1 x K) 归一化权重

    if nargin < 3 || isempty(eps_val)
        eps_val = 1e-6;
    end

    G = numel(dist_vec);
    K = min(K, G);  % 防止 K > G

    [sorted_dist, sorted_idx] = sort(dist_vec, 'ascend');

    neighbor_idx  = sorted_idx(1:K);
    neighbor_dist = sorted_dist(1:K);

    % 逆距离权重
    w = 1 ./ (neighbor_dist + eps_val);

    % 归一化
    weights = w / sum(w);
end
