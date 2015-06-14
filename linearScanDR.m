function [rf_out, rf_steer, x, z] = linearScanDR(rf_in,acq_params,bf_params,type)
% [rf_out, ~, x, z] = linearScanDR(rf,acq_params,bf_params)
%
% Linear Scan DS beamforming code - Will Long. Latest revision: 4/2/15
% Inputs: 
% rf_in - element array data [array signal x array element x Tx event]
% acq_params - parameters include rx_pos, c, t0
% bf_params - parameters include x (tx_pos or A line lateral location)
% type - fixed window function for rcv element weights (default: boxcar)

tic

if nargin < 4
    type = 'boxcar';
end

x = bf_params.x;
z_ref = ((acq_params.t0+1:acq_params.t0+size(rf_in,1))/acq_params.fs)*acq_params.c;
z = z_ref/2;

[b,a]=butter(2,[.05 .95]);

% intialize matrices and arrays for speed
n_tx = size(rf_in,3);
n_rcv_chn = length(acq_params.rx_pos);
n_depth = length(z);
rf_out = zeros(length(z),n_tx);

if n_rcv_chn ~= size(rf_in,2) || n_tx ~= size(rf_in,3)
    disp('Mismatch in RF data.')
end
    

% define incremental distances for estimation of echo propagation time
dz = repmat(z',1,n_rcv_chn);
dx = repmat(acq_params.rx_pos,n_depth,1);
dr = sqrt(dz.^2+dx.^2);
t_samp = (dr+dz)./acq_params.c;

% define fixed apodizations
switch type
    case 'boxcar'
        w = ones(1,n_rcv_chn); %weighting function defined
    case 'hann'
        w = hann(n_rcv_chn)';
end      
    
w_mat = repmat(w,size(rf_in,1),1);

% perform DS beamforming with DR focusing
for j = 1:n_tx % iterate for every tx event
    % linear interpolation to extract response at focal point
    % (apply focal delay)
    rf_steer(:,:,j) = w_mat.*linearInterp(z_ref'/acq_params.c,squeeze(rf_in(:,:,j)),t_samp);
    
    % sum focal delayed signals
    rf_out(:,j) = sum(rf_steer(:,:,j),2);
end

rf_out(find(isnan(rf_out))) = 0; 
rf_out = filter(b,a,rf_out);
toc

