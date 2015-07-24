clear all
close all
clc
addpath ../accessory/
% load /getlab/wjl11/scratch/verasonics_data/canine_pw_072315/complete_pwAcq_0deg_1n_100depth_20150723_114744_1.mat

load ~/Documents/MATLAB/data_files/beamforming/complete_pwAcq_5deg_5n_100depth_20150723_115306_0.mat
%%
tic
angle = rfdata.steerAngles;
fs = rfdata.samplingRateMHz*1e6;
c = rfdata.c;
t0 = rfdata.timeZero; 
elspacing = rfdata.elementSpacingMM*1e-3;
rf = rfdata.data;

rref = (t0:t0+size(rf,1)-1).*c/(2*fs);

xrange = [-0.01 0.01];
zrange = [0 0.11];

nxgrid = 500;
nzgrid = 1200;

xpts = linspace(xrange(1),xrange(2),nxgrid);
zpts = linspace(zrange(1),zrange(2),nzgrid);

xel = (-(size(rf,2)-1)/2:(size(rf,2)-1)/2)*elspacing;
nel = length(xel);
txRcv = [1 600 700];
nTxRcv = length(txRcv);

rf_out = zeros(nzgrid,nxgrid,nTxRcv);
idx = 0;
for k = txRcv
    idx = idx+1;
    r_offset = xel(end)*abs(sind(angle(k)));
    data = squeeze(rf(:,:,k));

    tmp = repmat(xel,length(xpts),1); % element xpos for each img point (first length(xpts) all correspond to channel 1)
    xel_vec = tmp(:)'; clear tmp;
    xel_mat = repmat(xel_vec,length(zpts),1); 

    xpts_vec = repmat(xpts,1,nel);
    xpts_mat = repmat(xpts_vec,length(zpts),1);

    dx_mat = abs(xpts_mat-xel_mat);
    dz_mat = repmat(zpts',1,length(xpts)*nel);

    rrcv = sqrt(dz_mat.^2+dx_mat.^2);
    rtx = dz_mat.*cosd(angle(k));

    rtot = rrcv+rtx;
    ri = rtot./2; % accommodate for rref which is distance given 2 way propagation
    
    
    tmp = repmat(data,1,length(xpts));
    tmp = reshape(tmp,size(data,1),nel,length(xpts));
    tmp = permute(tmp,[1 3 2]);
    data_rep = reshape(tmp,size(data,1),nel*length(xpts));
    interpdat = linearInterp(rref-r_offset/2,data_rep,ri);
    tmp = reshape(interpdat,nzgrid,nxgrid,nel);
%     imagesc(squeeze(tmp(:,140,:)));
%     pause
    rf_out(:,:,idx) = sum(tmp,3);
end
rf_out(find(isnan(rf_out))) = 0;
toc