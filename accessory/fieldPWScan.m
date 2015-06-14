function [rf, rf_raw, bf_params, acq_params, flow_params] = fieldPWScan(pointx,pointz,angles)
% [rf, rf_raw, bf_params, acq_params] = fieldPWScan(0,30,[-5:5:5]);

addpath '~/Documents/MATLAB/field_ii/'
addpath '~/Documents/MATLAB/beamforming/'
field_init(-1)
rng(0);

% Transducer parameters
f0=2.5e6;
fs=120e6;
c=1540;
el_width=0.3e-3;
el_height=6/1000;
el_kerf=0.00e-3;
elv_foc = 20;
N_el=128;
BW=0.6;
set_field('c',c);
set_field('fs',fs);

pitch=el_kerf+el_width;
ap_size=N_el*(el_kerf+el_width);

% Derived parameters
lambda=c/f0;

n_sub_x=5;
n_sub_y=5;

% % Single element transmit and full receive
% % (single element transmit reduces on-axis energy)
% tx=xdc_linear_array(N_el,el_width,el_height,el_kerf,n_sub_x,n_sub_y,focus); 
% rx=xdc_linear_array(N_el,el_width,el_height,el_kerf,n_sub_x,n_sub_y,focus);
% 
% % Set impulse responses
% tc=gauspuls('cutoff',f0,BW); % Note: default BWR -6 dB and TPE -60 dB
% imp_resp=gauspuls(-tc:1/fs:tc,f0,BW);
% xdc_impulse(tx,imp_resp); 
% xdc_impulse(rx,imp_resp);
% 
% excitation=sin(2*pi*f0*(0:1/fs:2/f0));
% xdc_excitation (tx, excitation);

x=linspace(-pitch*N_el*0.5,pitch*N_el*0.5,N_el);

% pointz = 30;
% pointx = 0;

zpos = pointz'./1000;
xpos = pointx'./1000;
% alter x positions so that points are directly on a-line

amplitude = ones(length(zpos),1);

% simulation custom parameters for flow
flow_angle = 60;
PRF = 3e3;
vel = 0; % set flow velocity in [m/s]

dr = vel*1/PRF;
dz = -dr*cosd(flow_angle); % in m
dx = dr*sind(flow_angle);
n_img = 1; % number of images to generate
% angles = [-5:5:5]; % angles to SA

trial = repmat(angles,[1 n_img]);
n_tx = length(trial);


for nn = 1:n_tx
    
    angle = trial(nn);
    
    x_foc = sind(angle)*100000;
    z_foc = cosd(angle)*100000;
    focus=[x_foc elv_foc z_foc]./1000;
    % Single element transmit and full receive
    % (single element transmit reduces on-axis energy)
    tx=xdc_linear_array(N_el,el_width,el_height,el_kerf,n_sub_x,n_sub_y,focus); 
    rx=xdc_linear_array(N_el,el_width,el_height,el_kerf,n_sub_x,n_sub_y,focus);
    
    tdr_info = xdc_get(rx,'rect');
    rx_pos = unique(tdr_info(24,:));
    
    % Set impulse responses
    tc=gauspuls('cutoff',f0,BW); % Note: default BWR -6 dB and TPE -60 dB
    imp_resp=gauspuls(-tc:1/fs:tc,f0,BW);
    xdc_impulse(tx,imp_resp); 
    xdc_impulse(rx,imp_resp);

    excitation=sin(2*pi*f0*(0:1/fs:2/f0));
    xdc_excitation (tx, excitation);

    xdc_apodization (tx, 0, hanning(128)');
    xdc_apodization (rx, 0, hanning(128)');
    
    zpos = zpos+dz;
    xpos = xpos+dx; % move scatterers for flow (omit if targets stationary)
    
    position = [xpos zeros(length(zpos),1) zpos];
    xdc_times_focus(rx,0,zeros(1,N_el));
    [tmp st(nn)] = calc_scat_multi(tx,rx,position,amplitude);
    
    tx_delay=xdc_get(tx,'focus');
    
    zfoc_comp  = max(tx_delay)-min(tx_delay)*c;
    xfoc_comp = (max(rx_pos)-min(rx_pos))/2;
    
    angle_calc = rad2deg(atan2(zfoc_comp,xfoc_comp));
    if angle < 0
        angle_calc = -angle_calc;
    elseif angle == 0
        angle_calc = 0;
    end
    rf_tmp{nn} = tmp;
    
    acq_params.angle(nn) = angle_calc;
    xdc_free(rx); xdc_free(tx);
end

tshift=round(size(conv(conv(excitation,imp_resp),imp_resp),2)/2);
[rf, t0] = shift_times_multi(rf_tmp,st,fs);
t0 = t0-tshift;

% add noise to rf
SNR = 60;
rf_raw = rf;
rf = rf+10^(-SNR/20)*max(rf(:))*randn(size(rf));

field_end

acq_params.c = c;
acq_params.fs = fs;
acq_params.t0 = t0;
acq_params.rx_pos = rx_pos;
acq_params.n_img = n_img;
acq_params.angles = angles;
acq_params.f0 = f0;

flow_params.vel = vel;
flow_params.PRF = PRF;
flow_params.flow_angle = 60;

bf_params.x = x;

