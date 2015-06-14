% Minimum Variance Beamforming lab
% BME 265 - Lab 6
% Nick Bottenus
% 3/15/12

clear;clc;close all

%% Simulate single point data
[rf lat ax acq_params bf_params]=simscan(0);
SNR=60; %dB
rf_noise=rf+10^(-SNR/20)*max(rf(:))*randn(size(rf));
rf = rf_noise;
save('ref_data.mat','acq_params','bf_params','rf');

%% Show B-mode
env=abs(hilbert(squeeze(sum(rf_noise,2))));
env=db(env/max(env(:)));
figure(1)
imagesc(lat,ax,env,[-60 0]);colormap gray; axis image
xlabel('Lateral (mm)')
ylabel('Axial (mm)')
title('DS PSF')

%% Run MV beamforming
[mvrf W]=MVbeamform(rf_noise,32);

%% Show MV and comparisons
% Calculate envelope
env2=abs(hilbert(mvrf));
env2=db(env2/max(env2(:)));
figure(2)
imagesc(lat,ax,env2,[-60 0]);colormap gray; axis image
xlabel('Lateral (mm)')
ylabel('Axial (mm)')
title('MV PSF')

% Compare with DS
figure(3)
plot(lat,env(117,:),lat,env2(117,:))
xlabel('Lateral (mm)')
ylabel('Amplitude (dB)')
title('PSF')
legend('DS','MV',0)

% Display weighting vector
figure(4)
tr=63:65;
freq=23;
subplot(2,1,1)
plot(abs(W(:,tr,freq)))
xlabel('Channel')
ylabel('Magnitude')
title('Weighting vector')
legend('-1','0','1',0)

subplot(2,1,2)
plot(angle(W(:,tr,freq)))
xlabel('Channel')
ylabel('Phase')

%% Run MV beamforming array pattern
[mvrfs W]=MVbeamformarray(rf_noise,rf);

%% Show MV array pattern
% Calculate envelope
env3=abs(hilbert(mvrfs));
env3=db(env3/max(env3(:)));

% Compare with DS
figure(5)
plot(lat,env3(117,:),lat,env(117,:),'--')
xlabel('Lateral (mm)')
ylabel('Amplitude (dB)')
title('Minimum Variance array pattern')
legend('MV array pattern','DS',0)

%% Simulate three point data
[rf2 lat2 ax2]=simscan([-2 0 2]*1e-3);
rf2=rf2(1:size(rf,1),:,:);
SNR=60; %dB
rf_noise2=rf2+10^(-SNR/20)*max(rf2(:))*randn(size(rf2));

%% Show B-mode
env=abs(hilbert(squeeze(sum(rf_noise2,2))));
env=db(env/max(env(:)));
figure(6)
imagesc(lat,ax,env,[-60 0]);colormap gray; axis image
xlabel('Lateral (mm)')
ylabel('Axial (mm)')
title('DS PSF')

%% Run MV beamforming
[mvrf W]=MVbeamform(rf_noise2);

%% Show MV and comparisons
% Calculate envelope
env2=abs(hilbert(mvrf));
env2=db(env2/max(env2(:)));
figure(7)
imagesc(lat,ax,env2,[-60 0]);colormap gray; axis image
xlabel('Lateral (mm)')
ylabel('Axial (mm)')
title('MV PSF')

% Compare with DS
figure(8)
plot(lat,env(117,:),lat,env2(117,:))
xlabel('Lateral (mm)')
ylabel('Amplitude (dB)')
title('PSF')
legend('DS','MV',0)

% Display weighting vector
figure(9)
tr=63:65;
freq=23;
subplot(2,1,1)
plot(abs(W(:,tr,freq)))
xlabel('Channel')
ylabel('Magnitude')
title('Weighting vector')

subplot(2,1,2)
plot(angle(W(:,tr,freq)))
xlabel('Channel')
ylabel('Phase')

%% Run MV beamforming array pattern
[mvrfs W]=MVbeamformarray(rf_noise2,rf);

%% Show MV array pattern
% Calculate envelope
env3=abs(hilbert(mvrfs));
env3=db(env3/max(env3(:)));

% Compare with DS
figure(10)
plot(lat,env3(117,:),lat,env(117,:),'--')
xlabel('Lateral (mm)')
ylabel('Amplitude (dB)')
title('Minimum Variance array pattern')
legend('MV array pattern','DS',0)

%% Run MV beamforming with diag loading
[mvrfdiag W]=MVbeamformdiag(rf_noise);

%% Show MV and comparisons
% Calculate envelope
env3=abs(hilbert(mvrfdiag));
env3=db(env3/max(env3(:)));
figure(11)
imagesc(lat,ax,env3,[-60 0]);colormap gray; axis image
xlabel('Lateral (mm)')
ylabel('Axial (mm)')
title('MV-diagonal loading-PSF')

% Compare with DS
figure(12)
plot(lat,env(117,:),lat,env2(117,:),lat,env3(117,:))
xlabel('Lateral (mm)')
ylabel('Amplitude (dB)')
title('PSF')
legend('DS','MV-spatial','MV-diag',0)

%% Run MV beamforming with spatial averaging
[mvrfsp W]=MVbeamform(rf_noise,4);

%% Show MV and comparisons
% Calculate envelope
env3=abs(hilbert(mvrfsp));
env3=db(env3/max(env3(:)));
figure(13)
imagesc(lat,ax,env3,[-60 0]);colormap gray; axis image
xlabel('Lateral (mm)')
ylabel('Axial (mm)')
title('MV-spatial-4')

% Compare with DS
figure(14)
plot(lat,env(117,:),lat,env2(117,:),lat,env3(117,:))
xlabel('Lateral (mm)')
ylabel('Amplitude (dB)')
title('PSF')
legend('DS','MV-spatial-32','MV-spatial-4',0)

%% Run MV beamforming with temporal averaging and comparison spatial averaging
N=100;
cutoff=48;
rf_noisetemp=zeros(size(rf,1),size(rf,2)-2*cutoff,size(rf,3),N);
for i=1:N
    SNR=60; %dB
    rf_noisetemp(:,:,:,i)=rf(:,cutoff+1:end-cutoff,:)+...
        10^(-SNR/20)*max(rf(:))*randn(size(rf_noisetemp(:,:,:,1)));
end
[mvrftemp W]=MVbeamformtemp(rf_noisetemp);
[mvrfsp W]=MVbeamform(rf_noisetemp(:,:,:,1),8);

%% Show MV and comparisons
% Calculate envelopes
env=abs(hilbert(squeeze(sum(rf_noisetemp(:,:,:,1),2))));
env=db(env/max(env(:)));

env2=abs(hilbert(mvrfsp));
env2=db(env2/max(env2(:)));

env3=abs(hilbert(mvrftemp));
env3=db(env3/max(env3(:)));
figure(15)
imagesc(lat,ax,env3,[-60 0]);colormap gray; axis image
xlabel('Lateral (mm)')
ylabel('Axial (mm)')
title('MV-temporal')

% Compare with DS
figure(16)
plot(lat,env(117,:),lat,env2(117,:),lat,env3(117,:))
xlabel('Lateral (mm)')
ylabel('Amplitude (dB)')
title('PSF')
legend('DS','MV-spatial','MV-temporal',0)