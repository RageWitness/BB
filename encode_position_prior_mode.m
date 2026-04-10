function oh = encode_position_prior_mode(mode)
% ENCODE_POSITION_PRIOR_MODE  位置先验模式 → 3维 one-hot
%
%   oh = encode_position_prior_mode(mode)
%
%   编码：
%     'fixed_known'                 → [1 0 0]
%     'time_known_position_unknown' → [0 1 0]
%     'unknown'                     → [0 0 1]

    switch mode
        case 'fixed_known'
            oh = [1, 0, 0];
        case 'time_known_position_unknown'
            oh = [0, 1, 0];
        case 'unknown'
            oh = [0, 0, 1];
        otherwise
            oh = [0, 0, 1];
    end
end
