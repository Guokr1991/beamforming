function [rf_out, x, z] = linearScanMV_mod(rf_in,acq_params,bf_params,lines,flag)
% [rf_out, x, z] = linearScanMV_mod(rf,acq_params,bf_params)
%
% Dynamic receive code - Will Long. Latest revision: 10/1/14
% Inputs: 
% rf_in - raw rf data organized [rf_line,rx_chan, tx_event]
% acq_params - parameters include rx_pos, c, t0
% bf_params - parameters include x (tx_pos or A line lateral location)
% flag - '1' for beamformed rf input and '0' raw rf input
% 
% NOTE: for the case when all rx positions are relative and identical for
% each tx location (i.e. image using 128 subaperture of 256 array with tx
% focus at center of subaperture. Walk subaperture to generate multiple
% scan lines. In this case, lateral distance from focus is only the rx pos.
% 
% all focusing performed on axis
%
% acq_params.rx_pos; %lateral positions of each rx element relative to tx
% bf_params.x; %lateral positions of each tx focus (per tx event)
% acq_params.c;
% acq_params.t0; start time of rf data reference in terms of sample number

if nargin < 4 || isempty(lines)
    lines = 1:length(bf_params.x);
end
if nargin < 5
    flag = 0;
end

x = bf_params.x(lines);
z_ref = ((acq_params.t0+1:acq_params.t0+size(rf_in,1))/acq_params.fs)*acq_params.c;
z = z_ref/2;

[b,a]=butter(2,[.05 .95]);

% intialize matrices and arrays for speed
n_tx = length(x);         % # of transmit events
M = length(acq_params.rx_pos);      % # of receive channels per transmit
n_depth = length(z);
% rf_out = zeros(length(z),n_tx);

dz = repmat(z',1,M);
dx = repmat(acq_params.rx_pos,n_depth,1);
dr = sqrt(dz.^2+dx.^2);
t_samp = (dr+dz)./acq_params.c;

nZ = length(z); 
% nZ = 2^(nextpow2(nZ)-1);            % # of data points per frequency window
Mp = floor(M/4);                % # of subarray elements for subarray avg (Mp <= M/2)
fprintf('# elements in subarray: %d \n',Mp)
e = ones(Mp,1);                  % steering vector for pre-beamformed data

rf_bf = zeros(length(z),M,length(lines));
switch flag
    case 0 
        idx = 0;
        fprintf('Dynamic receive focusing... /n');
        for l = lines
            idx = idx+1;
            rf_bf(:,:,idx) = linearInterp(z_ref'/acq_params.c,squeeze(rf_in(:,:,l)),t_samp);
        end
    case 1
        rf_bf = rf_in;
end

idx = 0;
for l = lines
    fprintf('%d/%d A-line processed... \n',l,length(bf_params.x))
    idx = idx+1;
    rf_bf(find(isnan(rf_bf))) = 0; 
    Yl = fft(rf_bf(:,:,idx),nZ,1)';
        
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
        Bl(k) = wl'*Gav;                % beamform operation in k-domain
    end
    rf_out(:,idx) = ifft(Bl);             % inverse fft to retreive beamformed rf
end

% rf_out = filter(b,a,rf_out);
rf_out = flipud(rf_out);

