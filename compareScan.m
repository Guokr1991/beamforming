close all; clear all; clc

figure;
load('phtm_hanningDR.mat')
subplot(121)
rf2bmode(rf_out, x, z);
title('Hanning')
clear rf_out x z
load('phtm_boxcarDR.mat')
subplot(122)
rf2bmode(rf_out, x, z);
title('Boxcar')

figure;
load('sim_hanningDR.mat')
subplot(121)
rf2bmode(rf_out, x, z);
title('Hanning')
clear rf_out x z
load('sim_boxcarDR.mat')
subplot(122)
rf2bmode(rf_out, x, z);
title('Boxcar')