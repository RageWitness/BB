function oh = encode_time_pattern_type(time_pattern_type)
% ENCODE_TIME_PATTERN_TYPE  时间行为模式 → 3维 one-hot
%
%   oh = encode_time_pattern_type(time_pattern_type)
%
%   编码：
%     'scheduled' → [1 0 0]
%     'poisson'   → [0 1 0]
%     'unknown'   → [0 0 1]

    switch time_pattern_type
        case 'scheduled'
            oh = [1, 0, 0];
        case 'poisson'
            oh = [0, 1, 0];
        case 'unknown'
            oh = [0, 0, 1];
        otherwise
            oh = [0, 0, 1];  % 默认归为 unknown
    end
end
