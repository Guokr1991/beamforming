% % point target (PSF) simulation and MV beamform
% clear all; close all; clc;
% fprintf('Starting single point target simulation. \n')
% zpos = 45;
% xpos = 0;
% sector = 16;
% [rf, rf_raw, bf_params, acq_params] = fieldLinearScan(xpos,zpos,sector);
% save('point_target.mat','acq_params','bf_params','rf','rf_raw');
% [rf_out, x, z] = linearScanMV(rf,acq_params,bf_params,[],0);
% save('point_target_MV128fft.mat','rf_out','x','z','rf');
% fprintf('Finished. \n')
% 
% %% large grid simulation and MV beamform
% clear all; close all; clc;
% fprintf('Starting point grid simulation. \n')
% 
% % zpos = [35:10:55 30:5:60 42.5:2.5:47.5];
% % xpos = [-5*ones(1,3) zeros(1,7)  5*ones(1,3)];
% zpos = [42.5 45 47.5];
% xpos = zeros(1,length(zpos));
% sector = 14;
% [rf, rf_raw, bf_params, acq_params] = fieldLinearScan(xpos,zpos,sector);
% save('point_gridv.mat','acq_params','bf_params','rf','rf_raw');
% [rf_out, x, z] = linearScanMV(rf,acq_params,bf_params,[],0);
% save('point_gridv_MV128fft.mat','rf_out','x','z','rf');
% fprintf('Finished. \n')
% 
% % horizontal grid simulation and MV beamform
% clear all; close all; clc;
% fprintf('Starting point grid simulation. \n')
% 
% xpos = [-1 0 1];
% zpos = 45.*ones(1,length(xpos));
% sector = 14;
% [rf, rf_raw, bf_params, acq_params] = fieldLinearScan(xpos,zpos,sector);
% save('point_gridh1.mat','acq_params','bf_params','rf','rf_raw');
% [rf_out,~, x, z] = linearScanDR(rf,acq_params,bf_params);
% save('point_gridh1_DR.mat','rf_out','x','z'); clear rf_out x z
% [rf_out, x, z] = linearScanMVfast(rf,acq_params,bf_params,[],0);
% save('point_gridh1_MV128fft.mat','rf_out','x','z');
% fprintf('Finished. \n')
% 
% %% horizontal grid simulation and MV beamform
% clear all; close all; clc;
% fprintf('Starting point grid simulation. \n')
% 
% xpos = [-0.5 0 0.5];
% zpos = 45.*ones(1,length(xpos));
% sector = 14;
% [rf, rf_raw, bf_params, acq_params] = fieldLinearScan(xpos,zpos,sector);
% save('point_gridh05.mat','acq_params','bf_params','rf','rf_raw');
% [rf_out,~, x, z] = linearScanDR(rf,acq_params,bf_params);
% save('point_gridh05_DR.mat','rf_out','x','z'); clear rf_out x z
% [rf_out, x, z] = linearScanMVfast(rf,acq_params,bf_params,[],0);
% save('point_gridh05_MV128fft.mat','rf_out','x','z');
% fprintf('Finished. \n')
% 
% %% horizontal grid simulation and MV beamform
% clear all; close all; clc;
% fprintf('Starting point grid simulation. \n')
% 
% xpos = [-5 0 5];
% zpos = 45.*ones(1,length(xpos));
% sector = 16;
% [rf, rf_raw, bf_params, acq_params] = fieldLinearScan(xpos,zpos,sector);
% save('point_gridh5.mat','acq_params','bf_params','rf','rf_raw');
% [rf_out,~, x, z] = linearScanDR(rf,acq_params,bf_params);
% save('point_gridh5_DR.mat','rf_out','x','z'); clear rf_out x z
% [rf_out, x, z] = linearScanMVfast(rf,acq_params,bf_params,[],0);
% save('point_gridh5_MV128fft.mat','rf_out','x','z');
% fprintf('Finished. \n')
% 
% %% lesion simulation and MV beamform
% clear all; close all; clc;
% fprintf('Starting lesion simulation. \n')
% sector = 16;
% [rf, rf_raw, bf_params, acq_params] = fieldLinearScan([],[],sector);
% save('les.mat','acq_params','bf_params','rf','rf_raw');
% [rf_out,~, x, z] = linearScanDR(rf,acq_params,bf_params);
% save('les_DR.mat','rf_out','x','z'); clear rf_out x z
% [rf_out, x, z] = linearScanMVfast(rf,acq_params,bf_params,[],0);
% save('les_MV128fft.mat','rf_out','x','z');
% fprintf('Finished. \n')
% 
% %% steering simulation
% clear all; close all; clc;
% fprintf('Starting point steering simulation. \n')
% xpos = [0 5 -2];
% zpos = [40 42 45];
% sector = 16;
% [rf, rf_raw, bf_params, acq_params] = fieldLinearScan(xpos,zpos,sector);
% save('steer_example.mat','acq_params','bf_params','rf','rf_raw');
% [rf_out, rf_steer, x, z] = linearScanDR(rf_raw,acq_params,bf_params);
% save('steer_example_DR.mat','rf_out','rf_steer','x','z'); clear rf_out x z
% fprintf('Finished. \n')

% %% point grid simulation and MV beamform
% clear all; close all; clc;
% fprintf('Starting point grid simulation. \n')
% 
% xpos = [-2 0 2 -2 0 2 -2 0 2];
% zpos = [42.5.*ones(1,3) 45.*ones(1,3) 47.5.*ones(1,3)];
% sector = 14;
% [rf, rf_raw, bf_params, acq_params] = fieldLinearScan(xpos,zpos,sector);
% save('point_grid.mat','acq_params','bf_params','rf','rf_raw');
% [rf_out,~, x, z] = linearScanDR(rf,acq_params,bf_params);
% save('point_grid_DR.mat','rf_out','x','z'); clear rf_out x z
% [rf_out,~, x, z] = linearScanDR(rf,acq_params,bf_params,'hann');
% save('point_grid_DRhann.mat','rf_out','x','z'); clear rf_out x z
% [rf_out, x, z] = linearScanMVfast(rf,acq_params,bf_params,[],0);
% save('point_grid_MV128fft.mat','rf_out','x','z');
% fprintf('Finished. \n')

%% SSA MV beamform
clear all; close all; clc;
load('./ece582_dataset/SSA_points.mat');
fprintf('Data loaded. Starting SSA beamforming. \n')
[rf_out, x, z] = linearScanMVfast(rf,acq_params,bf_params,[],1);
save('ssapoints_MV128fft.mat','rf_out','x','z');
fprintf('Finished. \n')

%% SSA MV beamform
clear all; close all; clc;
load('./ece582_dataset/SSA_fetal.mat');
fprintf('Data loaded. Starting SSA beamforming. \n')
[rf_out, x, z] = linearScanMVfast(rf,acq_params,bf_params,[],1);
save('ssafetal_MV128fft.mat','rf_out','x','z');
fprintf('Finished. \n')
