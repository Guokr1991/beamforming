clear all; close all; clc

load('tmp_32fft.mat')
figure(1); subplot(311); plot(abs(hilbert(rf_out))), title('32 pt FFT'), axis tight
figure(2); subplot(311); plot(rf_out), title('32 pt FFT'), axis tight
figure(3); subplot(311); 
imagesc(20*log10(abs(hilbert(rf_out))./max(abs(hilbert(rf_out(:))))),[-60 0])
colormap gray, title('32 pt FFT'), axis tight

load('tmp_128fft.mat')
figure(1); subplot(312); plot(abs(hilbert(rf_out))), title('128 pt FFT'), axis tight
figure(2); subplot(312); plot(rf_out), title('128 pt FFT'), axis tight
figure(3); subplot(312); 
imagesc(20*log10(abs(hilbert(rf_out))./max(abs(hilbert(rf_out(:))))),[-60 0])
colormap gray, title('128 pt FFT'), axis tight

load('tmp_256fft.mat')
figure(1); subplot(313); plot(abs(hilbert(rf_out))), title('256 pt FFT'), axis tight
figure(2); subplot(313); plot(rf_out), title('256 pt FFT'), axis tight
figure(3); subplot(313); 
imagesc(20*log10(abs(hilbert(rf_out))./max(abs(hilbert(rf_out(:))))),[-60 0])
colormap gray, title('128 pt FFT'), axis tight