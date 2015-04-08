
addpath ../
%% plane wave fronts and steering (on axis wavefronts are directly incident)
close all; clear all;
load steer_example.mat
t = 1e6.*(acq_params.t0:acq_params.t0+size(rf_raw,1)-1).*1/acq_params.fs;
el_n = 1:size(rf,2);
tx_on = round(size(rf,3)/2);

tmp = rf_raw(1:249,:,tx_on);
rf_norm(1:249,:) = tmp./max(tmp(:));
tmp = rf_raw(250:700,:,tx_on);
rf_norm(250:700,:) = tmp./max(tmp(:));
tmp = rf_raw(701:end,:,tx_on);
rf_norm(701:size(rf_raw,1),:) = tmp./max(tmp(:));

figure(1)
subplot(121)
I1 = imagesc(el_n,t,rf_norm); colormap gray; axis square
xlabel('Element no.'), ylim([min(t) 62]);
ylabel('t (\mus)')
title('Raw Channel Data')

load steer_example_DR.mat
clear rf_norm tmp
tmp = rf_steer(1:249,:,tx_on);
rf_norm(1:249,:) = tmp./max(tmp(:));
tmp = rf_steer(250:700,:,tx_on);
rf_norm(250:700,:) = tmp./max(tmp(:));
tmp = rf_steer(701:end,:,tx_on);
rf_norm(701:size(rf,1),:) = tmp./max(tmp(:));

figure(1)
subplot(122)
I2 = imagesc(el_n,t,rf_norm); colormap gray; axis square
xlabel('Element no.'),ylim([min(t) 62]);
ylabel('t (\mus)')
title('Steered Channel Data')
%% lesion comparison figures
clear all; close all;

DR = 60;
zlims = [30 50]./1000;

load les_DR.mat
idx = find(z>=40/1000,1);
figure(2)
subplot(131)
env = rf2bmode(rf_out,DR,x,z,[],zlims);
title('DS (Boxcar)')

tmp = env(idx,:);
figure(3); hold all
plot(1000*x,20.*log10(tmp/max(tmp)),'k'); axis square

load les_DRhann.mat
figure(2)
subplot(132)
env = rf2bmode(rf_out,DR,x,z,[],zlims);
title('DS (Hanning)')

tmp = env(idx,:);
figure(3)
plot(1000*x,20.*log10(tmp/max(tmp)),'b'); axis square

load les_MV128fft.mat
figure(2)
subplot(133)
env = rf2bmode(rf_out,DR,x,z,[],zlims);
title('Minimum Variance')

tmp = env(idx,:);
figure(3)
plot(1000*x,20.*log10(tmp/max(tmp)),'r'); axis square
legend('DS (Boxcar)','DS (Hanning)','MV')
xlabel('Lateral Distance [mm]'); ylabel('Power [dB]')

%% PSF figures
clear all; close all

DR = 60;
xlims = [-4 4]./1000;

load point_target_MV128fft.mat
zlims = [min(z) max(z)];
figure(3)
subplot(313)
[env,~,x,z] = rf2bmode(rf_out,DR,x,z,xlims,zlims);
title('Minimum Variance')

idx = find(z>=45/1000,1);
tmp = env(idx,:);
figure(4); hold all
plot(1000*x,20.*log10(tmp/max(tmp)),'r'); axis square

load point_target_DR.mat
figure(3)
subplot(311)
[env,~,x,z] = rf2bmode(rf_out,DR,x,z,xlims,zlims);
title('DS (Boxcar)')

idx = find(z>=45/1000,1);
tmp = env(idx,:);
figure(4)
plot(1000*x,20.*log10(tmp/max(tmp)),'k'); axis square

load point_target_DRhann.mat
figure(3)
subplot(312)
[env,~,x,z] = rf2bmode(rf_out,DR,x,z,xlims,zlims);
title('DS (Hanning)')

idx = find(z>=45/1000,1);
tmp = env(idx,:);
figure(4)
plot(1000*x,20.*log10(tmp/max(tmp)),'b'); axis square
legend('MV','DS (Boxcar)','DS (Hanning)')
xlabel('Lateral Distance [mm]'); ylabel('Power [dB]')
hold off

%% point 1 & 5 cm horizontal point grid plots 

clear all; close all;

DR = 60;
xlims = [];
% [-4 4]./1000;
% zlims = [30 50]./1000;
label = {'05','1'};
fig = 4;
for i = 1:length(label)
    fig = fig+1;
    load(['point_gridh' label{i} '_MV128fft.mat'])
    zlims = [min(z) max(z)];
    figure(fig)
    subplot(313)
    [env,~,x,z] = rf2bmode(rf_out,DR,x,z,xlims,zlims);
    title('Minimum Variance')

    idx = find(z>=45/1000,1);
    tmp = env(idx,:);
    figure(fig+1); hold all
    plot(1000*x,20.*log10(tmp/max(tmp)),'r'); axis square

    load(['point_gridh' label{i} '_DR.mat'])
    figure(fig)
    subplot(311)
    [env,~,x,z] = rf2bmode(rf_out,DR,x,z,xlims,zlims);
    title('DS (Boxcar)')

    idx = find(z>=45/1000,1);
    tmp = env(idx,:);
    figure(fig+1)
    plot(1000*x,20.*log10(tmp/max(tmp)),'k'); axis square

    load(['point_gridh' label{i} '_DRhann.mat'])
    figure(fig)
    subplot(312)
    [env,~,x,z] = rf2bmode(rf_out,DR,x,z,xlims,zlims);
    title('DS (Hanning)')

    idx = find(z>=45/1000,1);
    tmp = env(idx,:);
    figure(fig+1)
    plot(1000*x,20.*log10(tmp/max(tmp)),'b'); axis square
    legend('MV','DS (Boxcar)','DS (Hanning)')
    xlabel('Lateral Distance [mm]'); ylabel('Power [dB]')
    hold off
    fig = fig+1;
end

%% vertical point grid plots with horizontal
 
clear all; close all;

DR = 60;
xlims = [];
fig = 9;

load(['point_gridv_MV128fft.mat'])
zlims = [min(z) max(z)];
figure(fig)
subplot(233)
[env,~,x,z] = rf2bmode(rf_out,DR,x,z,xlims,zlims);
title('Minimum Variance')

load(['point_gridv_DR.mat'])
figure(fig)
subplot(231)
[env,~,x,z] = rf2bmode(rf_out,DR,x,z,xlims,zlims);
title('DS (Boxcar)')

load(['point_gridv_DRhann.mat'])
figure(fig)
subplot(232)
[env,~,x,z] = rf2bmode(rf_out,DR,x,z,xlims,zlims);
title('DS (Hanning)')

load(['point_gridh5_MV128fft.mat'])
zlims = [min(z) max(z)];
figure(fig)
subplot(236)
[env,~,x,z] = rf2bmode(rf_out,DR,x,z,xlims,zlims);

load(['point_gridh5_DR.mat'])
figure(fig)
subplot(234)
[env,~,x,z] = rf2bmode(rf_out,DR,x,z,xlims,zlims);

load(['point_gridh5_DRhann.mat'])
figure(fig)
subplot(235)
[env,~,x,z] = rf2bmode(rf_out,DR,x,z,xlims,zlims);

