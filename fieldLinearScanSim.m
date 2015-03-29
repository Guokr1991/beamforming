close all; clear all; clc;
addpath('../field_ii')
field_init(-1)
% Transducer parameters
f0=5e6; BW=0.7; N_el=256;
elv_focus=0.04; elv_fnum=6.5; kerf_fraction=0.05;
z_focus = elv_focus
focus = [0 0 z_focus];
nr = 129;

% General parameters
c=1540; fs=100e6; % speed of sound c and ultrasound sampling f (default fs)

% Derived parameters
lambda=c/f0;
pitch=lambda/2;
el_height=elv_focus/elv_fnum; % element height (size in y direction)
el_width=(1-kerf_fraction)*pitch; % element width (x direction)
el_kerf=kerf_fraction*pitch;
n_sub_x=ceil(el_width/(lambda/4));
n_sub_y=ceil(el_height/(lambda/4));

tx=xdc_linear_array(N_el,el_width,el_height,el_kerf,n_sub_x,n_sub_y,focus);
rx=xdc_linear_array(N_el,el_width,el_height,el_kerf,n_sub_x,n_sub_y,focus);

% Set impulse responses
tc=gauspuls('cutoff',f0,BW); % Note: default BWR -6 dB and TPE -60 dB
imp_resp=gauspuls(-tc:1/fs:tc,f0,BW);
xdc_impulse(tx,imp_resp); 
xdc_impulse(rx,imp_resp);

apod = zeros(N_el,1);
apod(1:128) = 1;
el_no = (1:N_el)';

tdr_info = xdc_get(tx,'rect');
full_arr_x = unique(tdr_info(24,:));
rx_pos = full_arr_x((128-63):(128+64));

% z = 0.02:0.01:0.06;
z = [z_focus-0.005];
xpos = (el_width+el_kerf)/2;
position = [xpos' zeros(length(z),1) z'];
amplitude = ones(length(z),1);

figure
for nn = 1:nr
    apod_x = circshift(apod,nn-1);
    
    xdc_apodization(tx,0,apod_x');
    xdc_apodization(rx,0,apod_x');
    
    sub_arr_x = full_arr_x(find(apod_x == 1));
    x(nn) = mean([sub_arr_x(64) sub_arr_x(65)]);
    
    xdc_center_focus(tx,[x(nn) 0 0]);
    xdc_focus(tx,0,[x(nn) 0 z_focus]);
    
    xdc_center_focus(rx,[x(nn) 0 0]);
    xdc_focus(rx,0,[x(nn) 0 100]);
    
    clf
    subplot(211)
    plot(xdc_get(rx,'apo').*xdc_get(rx,'focus'))
    subplot(212)
    plot(xdc_get(rx,'apo'))
    pause(0.1)
    
    [tmp st(nn)] = calc_scat_multi(tx,rx,position,amplitude);
    rf{nn} = tmp(:,find(apod_x == 1));
    
end
[rf, t0] = shift_times_multi(rf,st,fs);
SNR = 60;
rf = rf+10^(-SNR/20)*max(rf(:))*randn(size(rf));
xdc_free(rx); xdc_free(tx);

field_end

acq_params.c = c;
acq_params.fs = fs;
acq_params.t0 = t0;
acq_params.rx_pos = rx_pos;

bf_params.x = x;

save('simdata_pt.mat','acq_params','bf_params','rf')