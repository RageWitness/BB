function test_mwm_channel()
% TEST_MWM_CHANNEL  统一 MWM 信道模型最小可复现验证
%   覆盖：
%     1) 无遮挡 LOS
%     2) 单建筑直阻 → direct vs diffraction
%     3) AP 在建筑内
%     4) 路径擦边/穿角点不重复计墙
%     5) 四个频点参数表

    fprintf('\n========== test_mwm_channel ==========\n');

    % 单元测试 1：参数表
    fprintf('\n[1] get_channel_params 四频点：\n');
    fs = [89e6, 96e6, 2.4e9, 2593e6];
    for i = 1:numel(fs)
        p = get_channel_params(fs(i));
        fprintf('  f=%6.0f MHz | PL0=%6.2f Lw=%.2f Lc=%.2f Dc=%.2f sigma=%.2f\n', ...
            p.f_MHz, p.PL0, p.Lw, p.Lc, p.Deltaic, p.sigma);
    end
    p_low  = get_channel_params(96e6);
    p_high = get_channel_params(2593e6);
    assert(p_low.Lw < p_high.Lw, '低频墙损应小于高频');
    p_24  = get_channel_params(2.4e9);
    assert(p_high.PL0 > p_24.PL0, '2593MHz 的 PL0 应高于 2.4GHz');
    fprintf('  PASS: 频点参数关系正确\n');

    % 单元测试 2：无建筑遮挡
    fprintf('\n[2] LOS 无遮挡：\n');
    out = compute_path_loss([10,10], [50,10], 96e6, struct([]));
    p = get_channel_params(96e6);
    pl_ref = p.PL0 + p.N1 * log10(min(40,p.db)/p.d0);
    fprintf('  PL=%.2f, 期望=%.2f, mode=%s, Nw=%d\n', out.PL_total, pl_ref, out.chosen_mode, out.wall_count);
    assert(abs(out.PL_total - pl_ref) < 1e-6 && out.wall_count == 0, 'LOS 应只等于多斜率');
    fprintf('  PASS\n');

    % 单元测试 3：一栋建筑挡中间
    fprintf('\n[3] 单建筑横挡：\n');
    bld = struct('xmin',20,'xmax',30,'ymin',-5,'ymax',5,'material','concrete');
    out = compute_path_loss([10,0], [50,0], 2.4e9, bld);
    fprintf('  Nw=%d, mode=%s, PL_dir=%.2f, PL_diff=%.2f, PL=%.2f\n', ...
        out.wall_count, out.chosen_mode, out.PL_direct, out.PL_diff, out.PL_total);
    assert(out.wall_count == 2, '单建筑直穿应为 2 堵墙');
    assert(~isinf(out.PL_diff), '应能找到绕射路径');
    assert(out.PL_total <= out.PL_direct + 1e-6, '应取最小');
    fprintf('  PASS\n');

    % 单元测试 4：AP 在建筑内
    fprintf('\n[4] AP 位于建筑内部：\n');
    bld2 = struct('xmin',40,'xmax',60,'ymin',40,'ymax',60,'material','concrete');
    out = compute_path_loss([10,50], [50,50], 96e6, bld2);
    fprintf('  Nw=%d (期望 1: 进入终点建筑), mode=%s, PL=%.2f\n', ...
        out.wall_count, out.chosen_mode, out.PL_total);
    assert(out.wall_count == 1, 'AP 在建筑内, src 在外, 应有 1 次进入边界穿越');
    fprintf('  PASS\n');

    % 单元测试 5：擦边
    fprintf('\n[5] 路径恰好沿建筑底边：\n');
    bld3 = struct('xmin',20,'xmax',30,'ymin',0,'ymax',10,'material','concrete');
    out = compute_path_loss([10,0], [50,0], 96e6, bld3);
    fprintf('  Nw=%d (期望 0), PL=%.2f\n', out.wall_count, out.PL_total);
    assert(out.wall_count == 0, '擦边不应计墙');
    fprintf('  PASS\n');

    % 集成测试：默认建筑布局 + 真实 AP，建库一次
    fprintf('\n[6] 集成：默认建筑 + 全频点单点建库：\n');
    blds = get_default_buildings();
    fprintf('  建筑数: %d\n', numel(blds));
    src = [150, 150];
    aps = [
        65.2850   13.7581;
        14.0990   186.6792;
        192.2744  255.9371;
        203.5478  25.5561;
        273.7034  98.0424;
        67.5742   113.4034;
        288.0000  215.4216;
        75.6478   274.5550
    ];
    for i = 1:numel(fs)
        [PL, sigma] = compute_pathloss_mwm(src, aps, fs(i), blds);
        fprintf('  f=%6.0f MHz: PL=[%.1f, %.1f] (range %.1f), sigma=%.2f\n', ...
            fs(i)/1e6, min(PL), max(PL), max(PL)-min(PL), sigma);
        assert(all(PL > 0 & PL < 300), 'PL 范围异常');
    end
    fprintf('  PASS\n');

    fprintf('\n========== ALL TESTS PASSED ==========\n');
end
