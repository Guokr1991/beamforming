function [rf, bf_params, acq_params] = fieldWSRx(pointx,pointz,sec)

% NOT FINISHED... beamformer is different from linearScanDR
if nargin < 3
    pointx = [];
    pointz = [];
    sec = 14;
end

addpath('../field_ii')
field_init(-1)
rng(0);

% Transducer parameters
f0=7e6;
fs=100e6;
c=1540;
el_width=0.075e-3;
el_height=6/1000;
el_kerf=0.035e-3;
focus=[0 0 45]./1000;
N_el=128;
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

% Scan sector
tdr_info = xdc_get(rx,'rect');
zfoc=focus(3);

rx_pos = unique(tdr_info(24,:));
xposlim = [min(rx_pos)-pitch/2, max(rx_pos)+pitch/2];
zposlim = [zfoc-20/1000, zfoc+20/1000];


lat_res = lambda*zfoc/ap_size;
ax_res = lambda;
aResCell = lat_res*ax_res;

scatPerResCell = 5; % minimum for ~ fully developed speckle
scatDensity = scatPerResCell/aResCell;
scatN = round((diff(xposlim)*diff(zposlim))*scatDensity);

xpos = xposlim(1)+(xposlim(2)-xposlim(1)).*rand(scatN,1);
zpos = zposlim(1)+(zposlim(2)-zposlim(1)).*rand(scatN,1);

amplitude = ones(scatN,1);

N_sub = 81; % should be an odd number
lines = N_el-N_sub+1;
dx = pitch; 

for nn = 1:lines
    apo = [zeros(1,nn-1) ones(1,N_sub) zeros(1, N_el-N_sub-nn+1)];
    x(nn)=(nn-1-lines/2)*dx;
%     (nn:nn+N_subarray-1)
    position = [xpos zeros(length(zpos),1) zpos];
    
    xdc_center_focus(tx,[0 0 0]);
    xdc_focus(tx,0,[0 0 zfoc]);

    xdc_times_focus(rx,0,zeros(1,N_el));
    xdc_apodization(rx, 0, apo);    
    
    [tmp st(nn)] = calc_scat_multi(tx,rx,position,amplitude);
    rf_tmp{nn} = tmp;
    fprintf('Line %d/%d simulated.\n',nn,lines);
end

[rf, t0] = shift_times_multi(rf_tmp,st,fs);
t0 = t0-tshift;

field_end

acq_params.c = c;
acq_params.fs = fs;
acq_params.t0 = t0;
acq_params.rx_pos = rx_pos;

bf_params.x = x;
