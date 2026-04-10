function [SpatialFP, SignatureLib] = init_m25_single_source_fp( ...
    APs, Bands, GridValid, Config, SourceTemplates)
% INIT_M25_SINGLE_SOURCE_FP  M2.5 模块总入口
%
%   [SpatialFP, SignatureLib] = init_m25_single_source_fp(
%       APs, Bands, GridValid, Config, SourceTemplates)
%
%   功能：
%     1. 构建 SpatialFP  — 单源空间指纹库（供 M4 WKNN 定位）
%     2. 构建 SignatureLib — 静态模板特征库（供 M3/M5/M6 接口）
%     3. 校验输出完整性
%
%   输入：
%       APs             - AP 结构体
%       Bands           - 频带结构体
%       GridValid       - 有效网格
%       Config          - 配置
%       SourceTemplates - M0 的模板库
%
%   输出：
%       SpatialFP       - 空间指纹库
%       SignatureLib    - 静态模板特征库

    fprintf('\n============================================\n');
    fprintf('  M2.5 指纹库与模板库构建\n');
    fprintf('============================================\n\n');

    %% 1. 构建 SpatialFP
    tic;
    SpatialFP = build_spatial_fp_single_source(APs, Bands, GridValid, Config);
    t_fp = toc;
    fprintf('[M2.5] SpatialFP 耗时 %.2f s\n\n', t_fp);

    %% 2. 构建 SignatureLib
    SignatureLib = build_signature_lib_single_source(SourceTemplates, Config);

    %% 3. 校验
    fprintf('\n');
    validate_m25_outputs_single_source(SpatialFP, SignatureLib);

    fprintf('\n============================================\n');
    fprintf('  M2.5 构建完成\n');
    fprintf('============================================\n');
end
