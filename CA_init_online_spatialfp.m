function OnlineSpatialFP = CA_init_online_spatialfp(OfflineSpatialFP)
% CA_INIT_ONLINE_SPATIALFP  从离线指纹库初始化在线副本

    OnlineSpatialFP = OfflineSpatialFP;
    OnlineSpatialFP.is_online_map   = true;
    OnlineSpatialFP.update_round    = 0;
    OnlineSpatialFP.update_history  = [];
end
