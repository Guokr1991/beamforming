function [rf_out, x, z, rf_steer] = sweptScanMV(rf_steer,x_in,z_in,lines,Mp,type)
% [rf_out, x_out, z_out] = sweptScanMV(rf,x,z);
%
% SSA MV beamforming code - Will Long. Latest revision: 4/24/15
% Inputs: 
% rf_steer - pre-focused SSA rf data [array signal x array element x Tx event]
% x_in - lateral positions of beamformed ROI
% z_in - axial positions of beamformed ROI
% lines - A-lines to beamform (corresponding to x_in)
% Mp - subarray length for covariance matrix estimation
% type - flag 'bf' to beamform data, 'p' to skip beamforming and output
% steered data only


if nargin < 4 || isempty(lines)
    lines = 1:size(rf_steer,3); % beamform all a-lines (all tx events)
end
if nargin < 5 || isempty(Mp)
    Mp = floor(size(rf_steer,2)/4); % # of subarray elements for subarray avg (Mp <= M/2)
end  
if nargin < 6 || isempty(type)
    type = 'bf';
end

x = x_in;
z = z_in; 

M = size(rf_steer,2); % total number of array elements
rf_steer(isnan(rf_steer)) = 0;

% initialize MV beamformer parameters
e = ones(Mp,1); % steering vector (all ones for planar wavefront)
fprintf('# elements for subarray avg: %d \n',Mp)

Nwin = 128; % fft window size (should be larger than 2 way conv of pulse)
if Nwin > length(z)
    error('Specified fft window size greater than available data'); 
end

% matrix of overlapping short time windows pre-defined for each desired focal depth
tmp = buffer(1:length(z),Nwin,Nwin-1); 
izwin = tmp(:,Nwin:end); clear tmp;
iz0 = izwin(Nwin/2,:);
switch type
    case 'bf'
        % memory pre-allocation for speed
        Bl = zeros(1,Nwin);
        rf_out = zeros(length(iz0),length(lines));

        % minimum variance beamform at each depth and for each subband
        idx = 0;
        for l = lines
            fprintf('Beamforming %d/%d A-line... \n',l,size(rf_steer,3))
            tic
            idx = idx+1;

            % calculate STDFT for subband beamforming
            Yl_mat = zeros(Nwin,length(iz0),size(rf_steer,2));
            for rcv = 1:size(rf_steer,2)
                tmp = rf_steer(izwin(:),rcv,l);
                rcv_mat = reshape(tmp,size(izwin)); 
                % take DFT of each window representing a different depth
                % (window values x depths x array element number)
                Yl_mat(:,:,rcv) = fft(rcv_mat,Nwin,1);
            end

            for zi = 1:length(iz0)
                % extract DFT of window at specfic depth for all array signals
                Yl = squeeze(Yl_mat(:,zi,:)).';

                % solve MV problem for each frequency band k
                for k = 1:Nwin
                    Rl = zeros(Mp,Mp);
                    Gav = zeros(Mp,1);
                    for p = 1:M-Mp+1
                        G = Yl(p:p+Mp-1,k); % subband values in subarray p
                        Rl = G*G'+Rl; % calculate subarray covariance matrix
                        Gav = G+Gav;
                    end
                    Rl = 1/Mp*Rl; % spatial smoothing of covariance matrix
                    Gav = 1/Mp*Gav; % subarray average of subband values
                    wl = (Rl\e)/(e'*(Rl\e)); % solution of MV problem
                    Bl(k) = wl'*Gav; % beamform in Fourier domain
                end
                % IDFT to extract beamformed signal
                bl = ifft(Bl);
                % extract beamformed signal at focal point
                rf_out(zi,idx) = bl(floor(Nwin/2));
            end
            toc
        end
    case 'p'
        fprintf('Skipping beamforming \n');
        rf_out = [];
end
z = z(iz0);