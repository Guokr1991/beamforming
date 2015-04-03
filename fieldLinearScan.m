function [rf, rf_raw, bf_params, acq_params] = fieldLinearScan(xpos,zpos,sec)
addpath('../field_ii')
field_init(-1)

% Transducer parameters
f0=7e6;
fs=100e6;
c=1540;
el_width=0.075e-3;
el_height=6/1000;
el_kerf=0.035e-3;
focus=[0 0 45]./1000;
N_el=129;
BW=0.6;
set_field('c',c);
set_field('fs',fs);

pitch=el_kerf+el_width;
ap_size=N_el*(el_kerf+el_width);

% General parameters
c=1540; fs=100e6; % speed of sound c and ultrasound sampling f (default fs)

% Derived parameters
lambda=c/f0;

n_sub_x=5;
n_sub_y=5;

% Single element transmit and full receive
% (single element transmit reduces on-axis energy)
tx=xdc_linear_array(N_el,el_width,el_height,el_kerf,n_sub_x,n_sub_y,focus); 
rx=xdc_linear_array(N_el,el_width,el_height,el_kerf,n_sub_x,n_sub_y,focus);

% Set impulse responses
tc=gauspuls('cutoff',f0,BW); % Note: default BWR -6 dB and TPE -60 dB
imp_resp=gauspuls(-tc:1/fs:tc,f0,BW);
xdc_impulse(tx,imp_resp); 
xdc_impulse(rx,imp_resp);

% Define excitation
excitation=sin(2*pi*f0*(0:1/fs:2/f0));
xdc_excitation(tx,excitation);

tshift=size(conv(conv(excitation,imp_resp),imp_resp),2)/2;

% z = 0.02:0.01:0.06;
% zpos = [35:10:45 30:5:60 42.5:2.5:47.5]./1000;
% xpos = [-5*ones(1,3) zeros(1,6)  5*ones(1,3)]./1000;

zpos = zpos./1000;
xpos = xpos./1000;
position = [xpos' zeros(length(zpos),1) zpos'];
amplitude = ones(length(zpos),1);

sector=sec/1000;
zfoc=focus(3);
lines=8*ceil(sector/(lambda*zfoc/ap_size))+1; % scan lines odd to ensure 
x=linspace(-sector/2,sector/2,lines);

tdr_info = xdc_get(rx,'rect');
rx_pos = unique(tdr_info(24,:));

for nn = 1:lines
    xpos_shift = xpos+x(nn);
    position = [xpos_shift' zeros(length(zpos),1) zpos'];
    fprintf('Generating A-line %d/%d at %.1f mm \n',nn,lines,...
        1000*xpos_shift);
    xdc_times_focus(rx,0,zeros(1,N_el));
    [tmp st(nn)] = calc_scat_multi(tx,rx,position,amplitude);
    rf_tmp{nn} = tmp;
end

[rf, t0] = shift_times_multi(rf_tmp,st,fs);
t0 = t0-tshift;
SNR = 60;

rf_raw = rf;
rf = rf+10^(-SNR/20)*max(rf(:))*randn(size(rf));
xdc_free(rx); xdc_free(tx);

field_end

acq_params.c = c;
acq_params.fs = fs;
acq_params.t0 = t0;
acq_params.rx_pos = rx_pos;

bf_params.x = x;

