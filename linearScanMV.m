function [rf_out, x, z] = linearScanMV(rf_in,acq_params,bf_params)
% [rf_out, x, z] = linearScanMV(rf,acq_params,bf_params)
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

nZ = 64; 
nZ = 2^nextpow2(nZ);            % # of data points per frequency window
nz0 = nZ/2:length(z)-nZ/2;
Mp = floor(M/4);                % # of subarray elements for subarray avg (Mp <= M/2)
e = ones(Mp,1);                  % steering vector for pre-steered data

% initialize matrices
Yl = zeros(M,nZ,length(nz0));
rf_mv = zeros(length(nz0),M,n_tx);
% covG = zeros(Mp,nZ,M-Mp); 

for l = 1:n_tx 
    fprintf('%d/%d tx line processed... \n',l,n_tx)
    rf_bf(:,:,l) = linearInterp(z_ref'/acq_params.c,squeeze(rf_in(:,:,l)),t_samp);
    for z0_i = 1:length(nz0)
        z0_win = nz0(z0_i)-nZ/2+1:nz0(z0_i)+nZ/2;
        Yl(:,:,z0_i) = fft(rf_bf(z0_win,:,l),[],1)';
        for p = 1:M-Mp+1
            G = Yl(p:p+Mp-1,:,z0_i);
            covG(:,:,p) = G*G';
        end
        Rl(:,:,z0_i) = sum(covG,3);
        tmp = repmat(real(inv(Rl(:,:,z0_i))*e/(e'*inv(Rl(:,:,z0_i))*e))',M/Mp,1);
        Wl(z0_i,:) = tmp(:);
        wl = 
        rf_mv(z0_i,:,l) = rf_bf(nz0(z0_i),:,l).* wl(z0_i,:); 
    end
    rf_out(:,l) = sum(rf_mv(:,:,l),2);
end

%         Y{z0_i,l} = fft(rf_bf(z0_win,:,l),[],1); % fft at each rcv elem 
% 
%         covar(Y_l(:,:,z0_i))

% for j = 1:n_tx % iterate for every tx event
% %     w = ones(1,M); %weighting function defined 
% %     w = hann(M)';
% %     w_mat = repmat(w,size(rf_in,1),1);
%     rf_out(:,j) = sum(w_mat.*linearInterp(z_ref'/acq_params.c,squeeze(rf_in(:,:,j)),t_samp),2);
% end
% 

% rf_out = filter(b,a,rf_out);
z = z(nz0);
rf_out(find(isnan(rf_out))) = 0; 


