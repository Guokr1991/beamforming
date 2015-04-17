
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
clear all; close all; clc

DR = 50;
zlims = [30 50]./1000;

load les_DR.mat
idx = find(z>=40/1000,1);
figure(2)
subplot(131)
[env,~,xpos,zpos] = rf2bmode(rf_out,DR,x,z,[],zlims);
title('DS (Boxcar)')

tmp = env(idx,:);
figure(3); hold all
plot(1000*x,20.*log10(tmp/max(tmp)),'k'); axis square

% masks for contrast calculations
xi = []; zi = [];
for ii = 1:length(xpos)
    kk = find((xpos(ii)^2+(zpos-0.04).^2)<0.005^2);
    xi = [xi ii.*ones(1,length(kk))];
    zi = [zi kk];
    clear kk
end

mask_les = zeros(size(env));
mask_bg = ones(size(env));
for nn = 1:length(xi)
    mask_bg(zi(nn),xi(nn)) = 0;
    mask_les(zi(nn),xi(nn)) = 1;
end

tmp = mask_les.*env;
env_les = tmp(find(tmp ~= 0));
tmp = mask_bg.*env;
env_bg = tmp(find(tmp ~= 0));
CNR(1) = calcContrast(env_les, env_bg);

load les_DRhann.mat
figure(2)
subplot(132)
env = rf2bmode(rf_out,DR,x,z,[],zlims);
title('DS (Hanning)')

tmp = env(idx,:);
figure(3)
plot(1000*x,20.*log10(tmp/max(tmp)),'b'); axis square
tmp = mask_les.*env;
env_les = tmp(find(tmp ~= 0));
tmp = mask_bg.*env;
env_bg = tmp(find(tmp ~= 0));
CNR(2) = calcContrast(env_les, env_bg);

load les_MV128fft.mat
figure(2)
subplot(133)
env = rf2bmode(rf_out,DR,x,z,[],zlims);
title('MV')

tmp = env(idx,:);
figure(3)
plot(1000*x,20.*log10(tmp/max(tmp)),'r'); axis square
legend('DS (Boxcar)','DS (Hanning)','MV','location','southeast')
xlabel('Lateral Distance [mm]'); ylabel('Power [dB]')

tmp = mask_les.*env;
env_les = tmp(find(tmp ~= 0));
tmp = mask_bg.*env;
env_bg = tmp(find(tmp ~= 0));
CNR(3) = calcContrast(env_les, env_bg);

fprintf('MV: CNR = %.2f \nDS (Boxcar): CNR = %.2f \nDS (Hanning): CNR = %.2f \n',...
    CNR(3),CNR(1),CNR(2));

%% lesion comparison figures with CF
clear all; close all; clc

DR = 50;
zlims = [30 50]./1000;

load les_DR.mat
idx = find(z>=40/1000,1);
figure(2)
subplot(141)
env = rf2bmode(rf_out,DR,x,z,[],zlims);
title('DS (Boxcar)')

tmp = env(idx,:);
figure(3); hold all
plot(1000*x,20.*log10(tmp/max(tmp)),'k'); axis square

load les_DRhann.mat
figure(2)
subplot(142)
env = rf2bmode(rf_out,DR,x,z,[],zlims);
title('DS (Hanning)')

tmp = env(idx,:);
figure(3)
plot(1000*x,20.*log10(tmp/max(tmp)),'b'); axis square

load les_MV128fft.mat
figure(2)
subplot(143)
env = rf2bmode(rf_out,DR,x,z,[],zlims);
title('MV')
tmp = env(idx,:);
figure(3)
plot(1000*x,20.*log10(tmp/max(tmp)),'g'); axis square

load les_MV128fft_CF.mat
figure(2)
subplot(144)
env = rf2bmode(rf_cf,DR,x,z,[],zlims);
title('MV + CF')

tmp = env(idx,:);
figure(3)
plot(1000*x,20.*log10(tmp/max(tmp)),'r'); axis square
legend('DS (Boxcar)','DS (Hanning)','MV','MV + CF','Location','southwest')
xlabel('Lateral Distance [mm]'); ylabel('Power [dB]')

%% PSF figures
clear all; close all; clc

DR = 50;
xlims = [-4 4]./1000;
ylims = [-120 0];
factor = 20; % interpolation factor
z0 = 45;

load point_target_MV128fft.mat
zlims = [min(z) max(z)];
figure(3)
subplot(313)
[env,~,x,z] = rf2bmode(rf_out,DR,x,z,xlims,zlims);
title('MV')

idx = find(z>=z0/1000,1);
tmp = env(idx,:);
tmpdb = 20.*log10(tmp/max(tmp));
tmpdb = interp(tmpdb,factor);
x = interp(x,factor);
fwhm(1) = 2*abs(x(find(tmpdb>=-6,1)))*1e3;
lobes = sort(findpeaks(tmpdb),'descend')
PSL(1) = lobes(find(lobes<-10,1));
figure(4); hold all
plot(1000*x,tmpdb,'r'); axis square; xlim(xlims.*1000); ylim(ylims);

load point_target_DR.mat
figure(3)
subplot(311)
[env,~,x,z] = rf2bmode(rf_out,DR,x,z,xlims,zlims);
title('DS (Boxcar)')

idx = find(z>=z0/1000,1);
tmp = env(idx,:);
tmpdb = 20.*log10(tmp/max(tmp));
tmpdb = interp(tmpdb,factor);
x = interp(x,factor);
fwhm(2) = 2*abs(x(find(tmpdb>=-6,1)))*1e3;
lobes = sort(findpeaks(tmpdb),'descend');
PSL(2) = lobes(find(lobes<-10,1));
figure(4)
plot(1000*x,tmpdb,'k'); axis square; xlim(xlims.*1000); ylim(ylims);

load point_target_DRhann.mat
figure(3)
subplot(312)
[env,~,x,z] = rf2bmode(rf_out,DR,x,z,xlims,zlims);
title('DS (Hanning)')

idx = find(z>=z0/1000,1);
tmp = env(idx,:);
tmpdb = 20.*log10(tmp/max(tmp));
tmpdb = interp(tmpdb,factor);
x = interp(x,factor);
fwhm(3) = 2*abs(x(find(tmpdb>=-6,1)))*1e3;
lobes = sort(findpeaks(tmpdb),'descend');
PSL(3) = lobes(find(lobes<-10,1));
figure(4)
plot(1000*x,tmpdb,'b'); axis square; xlim(xlims.*1000); ylim(ylims);
legend('MV','DS (Boxcar)','DS (Hanning)','Location','southwest')
xlabel('Lateral Distance [mm]'); ylabel('Power [dB]')
hold off

fprintf('MV: FWHM = %.3f mm \nDS (Boxcar): FWHM = %.3f mm \nDS (Hanning): FWHM = %.3f mm \n',...
    fwhm(1),fwhm(2),fwhm(3));

fprintf('MV: PSL = %d dB \nDS (Boxcar): PSL = %d dB \nDS (Hanning): PSL = %d dB \n',...
    round(PSL(1)),round(PSL(2)),round(PSL(3)));
%% point 1 & 5 cm horizontal point grid plots 

clear all; close all;

DR = 50;
ylims = [-110 0];
% [-4 4]./1000;
% zlims = [30 50]./1000;
label = {'05','1','5'};
fig = 4;
for i = 1:length(label)
    fig = fig+1;
    load(['point_gridh' label{i} '_MV128fft.mat'])
    zlims = [min(z) 0.047].*1000;
    figure(fig)
    subplot(313)
    [env,~,x,z] = rf2bmode(rf_out,DR,x,z);
    ylim(zlims)
    title('MV')

    idx = find(z>=45/1000,1);
    tmp = env(idx,:);
    figure(fig+1); hold all
    plot(1000*x,20.*log10(tmp/max(tmp)),'r'); 

    load(['point_gridh' label{i} '_DR.mat'])
    figure(fig)
    subplot(311)
    [env,~,x,z] = rf2bmode(rf_out,DR,x,z);
    ylim(zlims)
    title('DS (Boxcar)')

    idx = find(z>=45/1000,1);
    tmp = env(idx,:);
    figure(fig+1)
    plot(1000*x,20.*log10(tmp/max(tmp)),'k');

    load(['point_gridh' label{i} '_DRhann.mat'])
    figure(fig)
    subplot(312)
    [env,~,x,z] = rf2bmode(rf_out,DR,x,z);
    ylim(zlims)
    title('DS (Hanning)')

    idx = find(z>=45/1000,1);
    tmp = env(idx,:);
    figure(fig+1)
    plot(1000*x,20.*log10(tmp/max(tmp)),'b'); ylim(ylims); xlim([-6 6]);
    if strcmp(label(i),'5'), xlim([-8 8]); end
    legend('MV','DS (Boxcar)','DS (Hanning)','Location','southeast')
    xlabel('Lateral Distance [mm]'); ylabel('Power [dB]')
    hold off
    fig = fig+1;
end

%% vertical point grid plots with horizontal grid (separate grids)
 
clear all; close all;

DR = 50;
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

%% horizontal and vertical point grid 
clear all; close all;

DR = 50;
xlims = [];
fig = 10;

load(['point_grid_MV128fft.mat'])
zlims = [min(z) max(z)];
figure(fig)
subplot(133)
[env,~,x,z] = rf2bmode(rf_out,DR,x,z,xlims,zlims);
title('Minimum Variance')

load(['point_grid_DR.mat'])
figure(fig)
subplot(131)
[env,~,x,z] = rf2bmode(rf_out,DR,x,z,xlims,zlims);
title('DS (Boxcar)')

load(['point_grid_DRhann.mat'])
figure(fig)
subplot(132)
[env,~,x,z] = rf2bmode(rf_out,DR,x,z,xlims,zlims);
title('DS (Hanning)')


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% SSA WIRE TARGET FIGURES %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SSA MV side-by-side comparison (wire targets)
clear all; close all; 

DR = 50;

figure
load SSA_points_MV128fft_Mp25.mat
xmv = x;
zmv = z;
subplot(212)
rf2bmode(rf_out,DR,xmv,zmv);
title('Minimum Variance')

load SSA_points.mat
rf_cont = squeeze(sum(rf,2));
subplot(211)
rf2bmode(rf_cont,DR,x,z);
xlim(1000.*[min(xmv) max(xmv)]);
ylim(1000.*[min(zmv) max(zmv)]);
title('Conventional')

%% SSA beam pattern comparisons
clear all; close all

% depth = 76;
% latlims = [5 15];
depth = 79.1;
latlims = [35 50];

figure
hold on
load SSA_points.mat
rf_cont = squeeze(sum(rf,2));
env =rf2bmode(rf_cont, [], x, z, [], [], 0);
zi = find(z*1000>=depth,1);
xi = find(x*1000>=latlims(1) & x*1000<=latlims(2));
tmp = env(zi,xi);
plot(x(xi),20.*log10(tmp./max(tmp(:))),'b');

load SSA_points_MV128fft_Mp25.mat
env = rf2bmode(rf_out, [], x, z, [], [], 0);
zi = find(z*1000>=depth,1);
xi = find(x*1000>=latlims(1) & x*1000<=latlims(2));
tmp = env(zi,xi);
plot(x(xi),20.*log10(tmp./max(tmp(:))),'r');
hold off

legend('Conventional','MV')
xlabel('Lateral Distance [mm]'); ylabel('Power [dB]')

%% SSA 6 & 12 dB contour comparisons
clear all; close all

% zlims = [74 80]./1e3;
% xlims = [5 15]./1e3;
% zlims = [113 119]./1e3;
% xlims = [8 16]./1e3;

zlims = [76 84]./1e3;
xlims = [36 46]./1e3;

figure
load SSA_points.mat
rf_cont = squeeze(sum(rf,2));
subplot(121)
hold on
[env,~,xwin, zwin] = rf2bmode(rf_cont, 50, x, z, xlims, zlims);
envdb = 20.*log10(env./max(env(:)));
contour(1e3.*xwin,1e3.*zwin,envdb,[-6 -6],'b'); axis image; 
contour(1e3.*xwin,1e3.*zwin,envdb,[-12 -12],'r');
contour(1e3.*xwin,1e3.*zwin,envdb,[-20 -20],'g');
set(gca,'YDir','reverse');
hold off
title('Conventional')

load SSA_points_MV128fft_Mp25.mat
subplot(122)
hold on
[env,~,xwin, zwin] = rf2bmode(rf_out, 50, x, z, xlims, zlims);
envdb = 20.*log10(env./max(env(:)));
contour(1e3.*xwin,1e3.*zwin,envdb,[-6 -6],'b'); axis image; 
contour(1e3.*xwin,1e3.*zwin,envdb,[-12 -12],'r');
contour(1e3.*xwin,1e3.*zwin,envdb,[-20 -20],'g');
set(gca,'YDir','reverse');
hold off
title('MV')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% SSA FETAL PHANTOM FIGURES %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SSA MV side-by-side comparison (phantom target)
clear all; close all; 

DR = 45;

figure
load SSA_fetal_MV128fft_Mp25.mat
xmv = x;
zmv = z;
subplot(122)
rf2bmode(rf_out,DR,xmv,zmv);
title('MV SSA')

load SSA_fetal.mat
rf_control = squeeze(sum(rf,2));
subplot(121)
rf2bmode(rf_control,DR,x,z);
xlim(1000.*[min(xmv) max(xmv)]);
ylim(1000.*[min(zmv) max(zmv)]);
title('Conventional SSA')

figure
subplot(122)
rf2bmode(rf_control,DR,x,z);
xlim(1000.*[min(xmv) max(xmv)]);
ylim(1000.*[min(zmv) max(zmv)]);
title('Swept Synthetic Aperture')

load ./SSA_datasets/focusedSingle_fetal.mat
subplot(121)
rf2bmode(rf_focused,DR,x,z);
xlim(1000.*[min(xmv) max(xmv)]);
ylim(1000.*[min(zmv) max(zmv)]);
title('Single Position Conventional')

%% SSA MV side-by-side comparison with Mp (phantom target)
clear all; close all; 

DR = 45;

figure
load SSA_fetal_MV128fft_Mp5.mat
xmv = x;
zmv = z;
subplot(142)
rf2bmode(rf_out,DR,xmv,zmv);
title('Sub-array average (Mp = 5)')

load SSA_fetal_MV128fft_Mp10.mat
subplot(143)
rf2bmode(rf_out,DR,xmv,zmv);
title('Sub-array average (Mp = 10)')

load SSA_fetal_MV128fft_Mp25.mat
subplot(144)
rf2bmode(rf_out,DR,xmv,zmv);
title('Sub-array average (Mp = 25)')

load SSA_fetal.mat
rf_cont = squeeze(sum(rf,2));
subplot(141)
rf2bmode(rf_cont,DR,x,z);
xlim(1000.*[min(xmv) max(xmv)]);
ylim(1000.*[min(zmv) max(zmv)]);
title('Conventional')

%% Zoomed SSA MV side-by-side Mp comparison (phantom target)
clear all; close all; 

DR = 45;

xlims = [10 30];
ylims = [68 88];

figure
load SSA_fetal_MV128fft_Mp5.mat
xmv = x;
zmv = z;
subplot(142)
rf2bmode(rf_out,DR,xmv,zmv);
xlim(xlims);
ylim(ylims);
title('Sub-array average (Mp = 5)')

load SSA_fetal_MV128fft_Mp10.mat
subplot(143)
rf2bmode(rf_out,DR,xmv,zmv);
xlim(xlims);
ylim(ylims);
title('Sub-array average (Mp = 10)')

load SSA_fetal_MV128fft_Mp25.mat
subplot(144)
rf2bmode(rf_out,DR,xmv,zmv);
xlim(xlims);
ylim(ylims);
title('Sub-array average (Mp = 25)')

load SSA_fetal.mat
rf_cont = squeeze(sum(rf,2));
subplot(141)
rf2bmode(rf_cont,DR,x,z);
xlim(xlims);
ylim(ylims);
title('Conventional')

%% SSA MV side-by-side comparison with CF (phantom target)
clear all; close all

DR = 45;
CFDR = 70;

figure
load SSA_fetal_MV128fft_Mp25.mat
xmv = x;
zmv = z;
subplot(132)
rf2bmode(rf_out,DR,xmv,zmv);
title('MV (Mp = 25)')

load SSA_fetal_MV128fft_Mp25_CF.mat
subplot(133)
rf2bmode(rf_cf,CFDR,xmv,zmv);
title('MV + CF (Mp = 25, -70 dB)')

load SSA_fetal.mat
rf_cont = squeeze(sum(rf,2));
subplot(131)
rf2bmode(rf_cont,DR,x,z);
xlim(1000.*[min(xmv) max(xmv)]);
ylim(1000.*[min(zmv) max(zmv)]);
title('Conventional')

%% SSA MV side-by-side comparison with CF (phantom target)
clear all; close all

DR = 45;
CFDR = 70;

figure
load SSA_fetal_MV128fft_Mp25.mat
xmv = x;
zmv = z;
subplot(141)
rf2bmode(rf_out,DR,xmv,zmv);
title('MV (Mp = 25)')

load SSA_fetal_gauss.mat
subplot(143)
rf2bmode(rf_win,DR,x,z);
title('MV (Gaussian Window)')
xlim(1000.*[min(xmv) max(xmv)]);
ylim(1000.*[min(zmv) max(zmv)]);

load SSA_fetal_gauss_CF.mat
subplot(144)
env = rf2bmode(rf_cf,CFDR,x,z);
title('Gaussian Window + CF (-70 dB)')
xlim(1000.*[min(xmv) max(xmv)]);
ylim(1000.*[min(zmv) max(zmv)]);

load SSA_fetal_MV128fft_Mp25_CF.mat
subplot(142)
rf2bmode(rf_cf,CFDR,xmv,zmv);
title('MV + CF (Mp = 25, -70 dB)')


