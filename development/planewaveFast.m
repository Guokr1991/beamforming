clear all
close all
clc
addpath ../accessory/
load ~/Google' Drive'/Trahey' Lab'/DTU_summer_school/flow_sa_project/point_target.mat
angle = 0;

tic

fs = acq_params.fs;
c = acq_params.c;
t0 = acq_params.t0; 

rref = (t0:t0+size(rf,1)-1).*c/(2*fs);

xrange = [-0.01 0.01];
zrange = [0.029 0.032];

nxgrid = 190;
nzgrid = 200;

xpts = linspace(xrange(1),xrange(2),nxgrid);
zpts = linspace(zrange(1),zrange(2),nzgrid);

xel = acq_params.rx_pos;
nel = length(xel);
nTxRcv = size(rf,3);

rf_out = zeros(nzgrid,nxgrid,nTxRcv);

for k = 1:nTxRcv
    data = squeeze(rf(:,:,k));
    angle = acq_params.angle(k);

    tmp = repmat(xel,length(xpts),1); % element xpos for each img point (first length(xpts) all correspond to channel 1)
    xel_vec = tmp(:)'; clear tmp;
    xel_mat = repmat(xel_vec,length(zpts),1); 

    xpts_vec = repmat(xpts,1,nel);
    xpts_mat = repmat(xpts_vec,length(zpts),1);

    dx_mat = abs(xpts_mat-xel_mat);
    dz_mat = repmat(zpts',1,length(xpts)*nel);

    rrcv = sqrt(dz_mat.^2+dx_mat.^2);
    rtx = dz_mat.*cosd(angle);

    rtot = rrcv+rtx;
    ri = rtot./2; % accommodate for rref which is distance given 2 way propagation
    
    
    tmp = repmat(data,1,length(xpts));
    tmp = reshape(tmp,size(data,1),nel,length(xpts));
    tmp = permute(tmp,[1 3 2]);
    data_rep = reshape(tmp,size(data,1),nel*length(xpts));
    interpdat = linearInterp(rref,data_rep,ri);
    
    tmp = reshape(interpdat,nzgrid,nxgrid,nel);
    rf_out(:,:,k) = sum(tmp,3);
end
toc