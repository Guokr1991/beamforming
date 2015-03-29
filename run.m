clear all; close all; clc;

load('ref_data.mat')
[rf_out, x, z] = linearScanMV_mod2(rf,acq_params,bf_params,[],1);
save('MV128fft.mat','rf_out','x','z','rf');