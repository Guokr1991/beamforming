function [rf_out, rf_pre, x, z] = linearScanDR(rf_in,acq_params,bf_params,type)
% [rf_out, ~, x, z] = linearScanDR(rf,acq_params,bf_params)
%
% Dynamic receive code - Will Long. Latest revision: 10/1/14
% Inputs: 
% rf_in - raw rf data organized [rf_line,rx_chan, tx_event]
% acq_params - parameters include rx_pos, c, t0
% bf_params - parameters include x (tx_pos or A line lateral location)
% type - window function for rcv element weights (default: boxcar)
% 
% NOTE: for the case when all rx positions are relative and identical for
% each tx location (i.e. image using 128 subaperture of 256 array with tx
% focus at center of subaperture. Walk subaperture to generate multiple
% scan lines. In this case, lateral distance from focus is only the rx pos.


%all focusing performed on axis

% acq_params.rx_pos; %lateral positions of each rx element relative to tx
% bf_params.x; %lateral positions of each tx focus (per tx event)
% acq_params.c;
% acq_params.t0; start time of rf data reference in terms of sample number
tic

if nargin < 4
    type = 'boxcar';
end

x = bf_params.x;
z_ref = ((acq_params.t0+1:acq_params.t0+size(rf_in,1))/acq_params.fs)*acq_params.c;
z = z_ref/2;

[b,a]=butter(2,[.05 .95]);

% intialize matrices and arrays for speed
n_tx = length(bf_params.x);
n_rcv_chn = length(acq_params.rx_pos);
n_depth = length(z);
rf_out = zeros(length(z),n_tx);

if n_rcv_chn ~= size(rf_in,2) || n_tx ~= size(rf_in,3)
    disp('Mismatch in RF data.')
end
    
dz = repmat(z',1,n_rcv_chn);
dx = repmat(acq_params.rx_pos,n_depth,1);
dr = sqrt(dz.^2+dx.^2);
t_samp = (dr+dz)./acq_params.c;
switch type
    case 'boxcar'
        w = ones(1,n_rcv_chn); %weighting function defined
    case 'hann'
        w = hann(n_rcv_chn)';
end      
    
w_mat = repmat(w,size(rf_in,1),1);
for j = 1:n_tx % iterate for every tx event
    rf_pre(:,:,j) = w_mat.*linearInterp(z_ref'/acq_params.c,squeeze(rf_in(:,:,j)),t_samp);
    rf_out(:,j) = sum(rf_pre(:,:,j),2);
end

rf_out(find(isnan(rf_out))) = 0; 
rf_out = filter(b,a,rf_out);
toc

