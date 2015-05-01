%% lesion comparison figures
clear all; close all; clc

DR = 50;
zlims = [30 50]./1000;

load les_DR.mat
idx = find(z>=40/1000,1);
figure(1)
subplot(131)
[env,~,xpos,zpos] = rf2bmode(rf_out,DR,x,z,[],zlims);
title('DS (Boxcar)'), grid off
ax = gca;
c = ax.Color;
ax.Color = 'white';
ax.FontSize = 16;
xlabel('x [mm]')
ylabel('z [mm]')

tmp = env(idx,:);
figure(2); hold all
plot(1000*x,20.*log10(tmp/max(tmp)),'k--'); axis square

load les_DRhann.mat
figure(1)
subplot(132)
env = rf2bmode(rf_out,DR,x,z,[],zlims);
title('DS (Hanning)'), grid off
ax = gca;
c = ax.Color;
ax.Color = 'white';
ax.FontSize = 16;
xlabel('x [mm]')
ylabel('z [mm]')

load les_MV128fft.mat
figure(1)
subplot(133)
env = rf2bmode(rf_out,DR,x,z,[],zlims);
title('MV'), grid off
set(gcf,'position',[100 100 700 300])
ax = gca;
c = ax.Color;
ax.Color = 'white';
ax.FontSize = 16;
xlabel('x [mm]')
ylabel('z [mm]')

tmp = env(idx,:);
figure(2)
plot(1000*x,20.*log10(tmp/max(tmp)),'k-'); axis square
h = legend('DS (Boxcar)','DS (Hanning)','MV','location','southeast');
set(h,'FontSize',12);
xlabel('Lateral Distance [mm]'); ylabel('Power [dB]')


%%
clear all; close all; clc

DR = 50;
xlims = [-2.5 2.5]./1000;
ylims = [-120 0];
factor = 20; % interpolation factor
z0 = 45;

load point_target_MV128fft.mat
zlims = [min(z) 46.5/1000];
figure(1)
subplot(313)
[env,~,x,z] = rf2bmode(rf_out,DR,x,z,xlims,zlims);
ylabel('Depth [mm]')
title('MV'), grid off
ax = gca;
ax.YTick = [44.5:0.5:50];
xlabel('x [mm]')
ylabel('z [mm]')

load point_target_DR.mat
figure(1)
subplot(311)
[env,~,x,z] = rf2bmode(rf_out,DR,x,z,xlims,zlims);
ylabel('Depth [mm]')
title('DS (Boxcar)'), grid off
ax = gca;
ax.YTick = [44.5:0.5:50];
xlabel('x [mm]')
ylabel('z [mm]')

load point_target_DRhann.mat
figure(1)
subplot(312)
[env,~,x,z] = rf2bmode(rf_out,DR,x,z,xlims,zlims);
ylabel('Depth [mm]')
title('DS (Hanning)'), grid off
ax = gca;
ax.YTick = [44.5:0.5:50];
xlabel('x [mm]')
ylabel('z [mm]')
set(gcf,'position',[100 100 400 400])

%% SSA MV side-by-side comparison with CF (phantom target)
clear all; close all

DR = 45;
CFDR = 70;

figure(1)


load SSA_fetal_MV128fft_Mp25.mat
xmv = x;
zmv = z;
subplot(132)
% figure(1)
rf2bmode(rf_out,DR,xmv,zmv);
title('MV'), grid off
% set(gcf,'position',[100 100 400 400])
ax = gca;
ax.FontSize = 12;
xlabel('x [mm]')
ylabel('z [mm]')


% load ./SSA_datasets/focusedSingle_fetal.mat
% % subplot(221)
% % % figure(2)
% rf2bmode(rf_focused,DR,x,z);
% xlim(1000.*[min(xmv) max(xmv)]);
% ylim(1000.*[min(zmv) max(zmv)]);
% title('Single Position SA'), grid off
% set(gcf,'position',[100 100 400 400])
% ax = gca;
% ax.FontSize = 20;


load SSA_fetal_MV128fft_Mp25_CF.mat
subplot(133)
% figure(3)
rf2bmode(rf_cf,CFDR,xmv,zmv);
title('MV + CF'), grid off
set(gcf,'position',[100 100 700 300])
ax = gca;
ax.FontSize = 12;
xlabel('x [mm]')
ylabel('z [mm]')


load SSA_fetal.mat
rf_cont = squeeze(sum(rf,2));
subplot(131)
% figure(4)
rf2bmode(rf_cont,DR,x,z);
xlim(1000.*[min(xmv) max(xmv)]);
ylim(1000.*[min(zmv) max(zmv)]);
title('DS'), grid off

% set(gcf,'position',[100 100 400 400])
ax = gca;
ax.FontSize = 12;
xlabel('x [mm]')
ylabel('z [mm]')

%%
clear all; close all; 

DR = 45;

xlims = [10 30];
ylims = [68 88];

load SSA_fetal_MV128fft_Mp25.mat
xmv = x;
zmv = z;
subplot(122)
rf2bmode(rf_out,DR,xmv,zmv);
xlim(xlims);
ylim(ylims);
grid off
axis off
title('MV','Fontsize',14,'Color','w')

load SSA_fetal.mat
rf_cont = squeeze(sum(rf,2));
subplot(121)
rf2bmode(rf_cont,DR,x,z);
xlim(xlims);
ylim(ylims);
grid off
axis off
title('DS','Fontsize',14,'Color','w')