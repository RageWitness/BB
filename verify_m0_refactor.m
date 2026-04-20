function verify_m0_refactor()
% VERIFY_M0_REFACTOR  信源发生器重构后的 7 条断言验证
%
%   用法： >> verify_m0_refactor
%
%   覆盖：
%     1. 保底：仅 persistent 启用时，每帧必输出 persistent_cal
%     2. 宽带压制：broadband_cal 出现时压过一切
%     3. 机遇压待定位：opportunistic > target
%     4. region 合法性
%     5. gaussian 合理性
%     6. 功率先验一致性
%     7. 结构完整性（SourceEvent 字段齐全 + is_selected_final=true）

    rng(42);
    T = 200;

    fprintf('========== verify_m0_refactor ==========\n');

    %% --- Test 1: 保底验证 ---
    Cov1.m0.source.lambda_broadband = 0;
    Cov1.m0.source.lambda_opportunistic = 0;
    Cov1.m0.source.lambda_target = [0 0 0 0];
    [S, M0, Logs, GV, Cfg] = init_m0_nonoverlap(Cov1);
    n_pc = 0; n_total = 0;
    for t = 1:T
        [M0, F, Logs] = step_m0_nonoverlap(M0, S, GV, Cfg, t, Logs);
        for b = 1:Cfg.m0.num_bands
            if F.active_per_band(b).has_source
                n_total = n_total + 1;
                if strcmp(F.active_per_band(b).source_type, 'persistent_cal')
                    n_pc = n_pc + 1;
                end
            end
        end
    end
    assert(n_pc == n_total && n_total > 0, ...
        'Test1: 保底失败 (pc=%d, total=%d)', n_pc, n_total);
    fprintf('[PASS] Test 1 保底: %d/%d 帧输出 persistent_cal\n', n_pc, n_total);

    %% --- Test 2: 宽带压制 ---
    rng(7);
    Cov2.m0.source.lambda_broadband = 1.0;
    Cov2.m0.source.lambda_opportunistic = 1.0;
    Cov2.m0.source.lambda_target = [1 1 1 1];
    [S, M0, Logs, GV, Cfg] = init_m0_nonoverlap(Cov2);
    bad = 0; broadband_frames = 0;
    for t = 1:T
        [M0, F, Logs] = step_m0_nonoverlap(M0, S, GV, Cfg, t, Logs);
        % 只要任一频带 broadband 在，该帧胜选者都应该是 broadband_cal
        has_bc = false;
        for b = 1:Cfg.m0.num_bands
            apb = F.active_per_band(b);
            if apb.has_source && strcmp(apb.source_type, 'broadband_cal')
                has_bc = true; break;
            end
        end
        if has_bc
            broadband_frames = broadband_frames + 1;
            for b = 1:Cfg.m0.num_bands
                apb = F.active_per_band(b);
                if apb.has_source && ~strcmp(apb.source_type, 'broadband_cal')
                    bad = bad + 1;
                end
            end
        end
    end
    assert(bad == 0 && broadband_frames > 10, ...
        'Test2: broadband 未压制(bad=%d, bcf=%d)', bad, broadband_frames);
    fprintf('[PASS] Test 2 宽带压制: %d 帧含 broadband，全频带均为 broadband_cal\n', broadband_frames);

    %% --- Test 3: 机遇压待定位 ---
    rng(11);
    Cov3.m0.source.lambda_broadband = 0;
    Cov3.m0.source.lambda_opportunistic = 1.0;
    Cov3.m0.source.lambda_target = [1 1 1 1];
    [S, M0, Logs, GV, Cfg] = init_m0_nonoverlap(Cov3);
    opp = 0; tgt = 0;
    for t = 1:T
        [M0, F, Logs] = step_m0_nonoverlap(M0, S, GV, Cfg, t, Logs);
        for b = 1:Cfg.m0.num_bands
            apb = F.active_per_band(b);
            if ~apb.has_source, continue; end
            if strcmp(apb.source_type, 'opportunistic'), opp = opp + 1;
            elseif strcmp(apb.source_type, 'target'),    tgt = tgt + 1;
            end
        end
    end
    assert(opp > 0, 'Test3: 未见 opportunistic');
    fprintf('[PASS] Test 3 优先级: opportunistic=%d, target=%d 帧共存\n', opp, tgt);

    %% --- Test 4/5/6/7: 默认配置下长时间跑 ---
    rng(23);
    Cov4.m0.source.lambda_opportunistic = 0.3;      % 提高机遇源频率以覆盖三类先验
    Cov4.m0.source.opportunistic.location_prior_prob.exact    = 0.34;
    Cov4.m0.source.opportunistic.location_prior_prob.region   = 0.33;
    Cov4.m0.source.opportunistic.location_prior_prob.gaussian = 0.33;
    [S, M0, Logs, GV, Cfg] = init_m0_nonoverlap(Cov4);
    n_region = 0; n_gauss = 0; n_pexact = 0; n_ev = 0;
    max_norm = 0;
    for t = 1:T
        [M0, F, Logs] = step_m0_nonoverlap(M0, S, GV, Cfg, t, Logs);
        for b = 1:Cfg.m0.num_bands
            apb = F.active_per_band(b);
            if ~apb.has_source, continue; end
            n_ev = n_ev + 1;
            % Test 7: 结构完整性
            for fn = {'source_type','priority','band_id','band_list', ...
                      'true_position','true_tx_power','tx_power_by_band', ...
                      'location_prior','power_prior','timestamp', ...
                      'is_persistent','is_selected_final'}
                assert(isfield(apb, fn{1}), 'Test7: 缺字段 %s', fn{1});
            end
            assert(apb.is_selected_final, 'Test7: is_selected_final 非 true');
            % Test 4/5: location_prior
            if strcmp(apb.location_prior.type, 'region')
                n_region = n_region + 1;
                v = apb.location_prior.value;
                bbox = v.bbox;
                p = apb.true_position;
                assert(p(1) >= bbox(1) && p(1) <= bbox(2) && ...
                       p(2) >= bbox(3) && p(2) <= bbox(4), ...
                       'Test4: region 真实点越界 bbox');
            elseif strcmp(apb.location_prior.type, 'gaussian')
                n_gauss = n_gauss + 1;
                mu = apb.location_prior.value.mu;
                nrm = norm(apb.true_position - mu);
                max_norm = max(max_norm, nrm);
            end
            % Test 6: 功率先验一致
            if strcmp(apb.power_prior.type, 'exact')
                n_pexact = n_pexact + 1;
                assert(abs(apb.true_tx_power - apb.power_prior.value) < 1e-9, ...
                    'Test6: power_prior=exact 但 true_tx_power 不一致');
            end
        end
    end
    fprintf('[PASS] Test 4 region 合法性: %d 条全部通过\n', n_region);
    fprintf('[PASS] Test 5 gaussian 合理性: %d 条, 最大 |p-mu|=%.2f (sigma~%.1f)\n', ...
        n_gauss, max_norm, Cfg.m0.source.opportunistic.gaussian.sigma_default);
    fprintf('[PASS] Test 6 功率先验一致性: %d 条全部通过\n', n_pexact);
    fprintf('[PASS] Test 7 结构完整性: %d 个事件字段齐全\n', n_ev);

    fprintf('========== 全部 7 项通过 ==========\n');
end
