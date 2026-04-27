function SpatialFP = CB_recompute_fingerprint_derivatives_from_dBm(SpatialFP, Config)
% CB_RECOMPUTE_FINGERPRINT_DERIVATIVES_FROM_DBM  Rebuild derived FP fields.

    if nargin < 2
        Config = struct();
    end

    fp_mode = 'both';
    if isfield(Config, 'm25') && isfield(Config.m25, 'fingerprint_mode')
        fm = Config.m25.fingerprint_mode;
        if strcmp(fm, 'legacy_only')
            fp_mode = 'legacy';
        end
    end
    if isfield(Config, 'm25') && isfield(Config.m25, 'keep_legacy_shape') && ...
            ~logical(Config.m25.keep_legacy_shape)
        fp_mode = 'rf_minmax';
    end

    for b = 1:SpatialFP.B
        F_dBm = SpatialFP.band(b).F_dBm;
        SpatialFP.band(b).F_lin = 10.^(F_dBm / 10);

        % In this project rf_raw is matched against observed dBm RSS.
        SpatialFP.band(b).RF_raw = F_dBm;

        SpatialFP.band(b) = precompute_centered_fp_wknn(SpatialFP.band(b), fp_mode);
        SpatialFP.band(b).RF_raw = F_dBm;
        SpatialFP.band(b).rf_raw_domain = 'dBm';
    end
end
