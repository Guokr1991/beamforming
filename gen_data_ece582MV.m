%% point target (PSF) simulation and MV beamform
clear all; close all; clc;

load('point_target.mat')
[rf_out, x, z] = linearScanMV(rf,acq_params,bf_params,[],0);
save('point_target_MV128fft.mat','rf_out','x','z','rf');

zpos = 45;
xpos = 0;
sector = 16;
[rf, rf_raw, bf_params, acq_params] = fieldLinearScan(xpos,zpos,sector);
save('point_target.mat','acq_params','bf_params','rf','rf_raw');
[rf_out, x, z] = linearScanMV_mod2(rf,acq_params,bf_params,[],0);
save('point_target_MV128fft.mat','rf_out','x','z','rf');

%% point grid simulation and MV beamform
clear all; close all; clc;

zpos = [35:10:45 30:5:60 42.5:2.5:47.5];
xpos = [-5*ones(1,3) zeros(1,6)  5*ones(1,3)];
sector = 20;
[rf, rf_raw, bf_params, acq_params] = fieldLinearScan(xpos,zpos,sector);
save('point_grid.mat','acq_params','bf_params','rf','rf_raw');
[rf_out, x, z] = linearScanMV(rf,acq_params,bf_params,[],0);
save('point_grid_MV128fft.mat','rf_out','x','z','rf');