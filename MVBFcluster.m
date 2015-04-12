function MVBFcluster(IDX,label)

load([label '.mat'])

[rf_out, x, z] = linearScanMVfast(rf,acq_params,bf_params,IDX,1);
rf_out = squeeze(rf_out);
save(sprintf('./data/%s_line%d.mat',label,IDX),'rf_out','x','z');