function [rf_out, x, z] = linearScanMVtest(rf_in,acq_params,bf_params)
% [rf_out, x, z] = linearScanMVtest(rf,acq_params,bf_params)
%
% Dynamic receive code - Will Long. Latest revision: 10/1/14
% Inputs: 
% rf_in - raw rf data organized [rf_line,rx_chan, tx_event]
% acq_params - parameters include rx_pos, c, t0
% bf_params - parameters include x (tx_pos or A line lateral location)
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

x = bf_params.x;
z_ref = ((acq_params.t0+1:acq_params.t0+size(rf_in,1))/acq_params.fs)*acq_params.c;
z = z_ref/2;

[b,a]=butter(2,[.05 .95]);

% intialize matrices and arrays for speed
n_tx = length(bf_params.x);         % # of transmit events
M = length(acq_params.rx_pos);      % # of receive channels per transmit
n_depth = length(z);
% rf_out = zeros(length(z),n_tx);
rf_bf = zeros(length(z),M,n_tx);

if M ~= size(rf_in,2) || n_tx ~= size(rf_in,3)
    disp('Mismatch in RF data.')
end
    
dz = repmat(z',1,M);
dx = repmat(acq_params.rx_pos,n_depth,1);
dr = sqrt(dz.^2+dx.^2);
t_samp = (dr+dz)./acq_params.c;

nZ = length(z); 
% nZ = 2^(nextpow2(nZ)-1);            % # of data points per frequency window
Mp = floor(M/4);                % # of subarray elements for subarray avg (Mp <= M/2)
e = ones(Mp,1);                  % steering vector for pre-steered data

for l = 1:n_tx 
    fprintf('%d/%d tx line processed... \n',l,n_tx)
    rf_bf(:,:,l) = linearInterp(z_ref'/acq_params.c,squeeze(rf_in(:,:,l)),t_samp);
    rf_bf(find(isnan(rf_bf))) = 0; 
    Yl = fft(rf_bf(:,:,l),nZ,1)';
        
    for k = 1:nZ
        Rl = zeros(Mp,Mp);
        Gav = zeros(Mp,1);
        for p = 1:M-Mp+1
            G = Yl(p:p+Mp-1,k);
            Rl = G*G'+Rl;
            Gav = G+Gav;
        end
        Rl = 1/Mp*Rl;
        Gav = 1/Mp*Gav;
        wl = inv(Rl)*e/(e'*inv(Rl)*e);
        Bl(k) = wl'*Gav; % beamform in frequency domain
    end
    rf_out(:,l) = ifft(Bl); 
end
rf_out(find(isnan(rf_out))) = 0; 
rf_out = flipud(rf_out);

