% SLSC - Short-Lag Spatial Coherence Imaging
%
% [CC] = SLSC(X,K) computes the van Cittert-Zernike theorem data in order to
% create an SLSC image.  X is the PxNxQ channel data and K is the size of
% the correlation kernel used to compute the spatial coherence.  P is the
% number of samples in the time/depth dimension, N is the number of
% channels, and Q is the number of image beams.  The output CC is the
% cross-correlation values as a function of element lag.
%
% [CC] = SLSC(X,K,MAXLAG) computes the same SLSC data above, where MAXLAG is
% the maximum spatial lag to compute the cross-correlations over (i.e. the
% largest spacing between elements to compute the cross-correlation).
%
% [CC] = SLSC(X,K,MAXLAG,DS) computes the same as above, where DS is the
% downsample factor.  DS forces the cross-correlations to be computed every
% DS depth sample.  If MAXLAG and DS are not supplied, they will default to
% 6 and 12, respectively.
%
% [CC] = SLSC(X,K,MAXLAG,DS,A) computes the same as above, where A is the
% maximum aperture size used to compute the SLSC data.  SLSC data will be
% computed for N-A apertures to be averaged together (default).
%
% [CC] = SLSC(X,K,MAXLAG,DS,A,MODEA) computes the same as above, where MODEA
% specifies the number of apertures.  For MODEA 'max', the maximum number of
% apertures (N-A), as defined previously, is used.  For MODEA 'min', only 1
% aperture is used, centered about the full aperture.  MODEA defaults to
% 'min' when A is undefined, and defaults to 'max' when A is defined.

% %%Common configuration Example:
% x = data; % Time-delayed channel signal data (usually a 3D matrix)
% k = numsamples; % Set to number of samples equal to a wavelength
% maxlag = 25; % Maximum aperture lag to compute correlations
% ds = 12; % For fast computation with a sampling frequency of 40 MHz
% a = numelements; % Set to maximum number of elements in array
% modea = 'max';
% cc = slsc(x,k,maxlag,ds,a,modea);
%
% %%To make an image
% M = 12; % Or whatever short-lag integration value you want
% c = squeeze(sum(cc(:,2:M,:),2)); % Integrate
% c = c./max(c(:)); % Normalize image
% imagesc(c,[0 1]) % Display, restrict dynamic range to [0,1])

% Change line with "for snapshot = 1:numsnapshots" to 
% "parfor snapshot = 1:numsnapshots" to take advantage
% of the multiple cpu/parallel processing code of matlab 

function [cc] = slsc(x,k,varargin)
    snapshot = 1; % to bypass matlab function

    [numsamp,num_el,numsnapshots] = size(x);

    if (nargin == 2) 
        maxlag = 6;
        ds = 12;
        a = num_el;
        modea = 'min';
    elseif (nargin == 3)
        maxlag = varargin{1};
        ds = 12;
        a = num_el;
        modea = 'min';
    elseif (nargin == 4)
        maxlag = varargin{1};
        ds = varargin{2};
        a = num_el;
        modea = 'min';
    elseif (nargin == 5)
        maxlag = varargin{1};
        ds = varargin{2};
        a = varargin{3};  
        modea = 'max';
    elseif (nargin == 6)
        maxlag = varargin{1};
        ds = varargin{2};
        a = varargin{3};  
        modea = varargin{4};
    else
        error('Too many/few input parameters!')
    end

    % Out samples in Z
    outsz = floor(numsamp/ds);

    % Check Kenel size
    if (mod(k,2)==0)
        k = k + 1;
        disp(['Kernel size is even, rounding up for odd sized kernel of ' ...
            num2str(k) ' samples'])
    end

    % Correlation indices
    cc = ones(outsz,maxlag,numsnapshots);

    % Half kernel size
    hk = (k-1)/2;

    for snapshot = 1:numsnapshots
        % Reserve memory space for output
        cctmp = ones(outsz,maxlag);
        for l = ds:ds:numsamp
            % Determine kernel for coherence computation
            samp = round(l/ds);
            g = [max(1,l-hk) min(l+hk,numsamp)];

            % Create sub-apertures
            if (strcmp(modea,'max') == 1)
                % Use maximum number of subapertures
                numap = num_el-a+1;
                rtmp = zeros(a,a,numap);
                for y = 1:numap
                    tmp = x(g(1):g(2),y:y+a-1,snapshot);
                    % Gated data
                    bt = tmp - repmat(mean(tmp,1),length(g(1):g(2)),1);
                    % Standard deviation
                    st = std(bt,0,1);
                    % Normalized cross correlation
                    rtmp(:,:,y) = (1/(k-1))*bt'*bt./(st'*st);
                end
                r = mean(rtmp,3);
            elseif (strcmp(modea,'min') == 1)
                % Use one subaperture
                ap = max(1,round(num_el/2)-round(a/2)):min(num_el,round(num_el/2)+round(a/2));
                tmp = x(g(1):g(2),ap,snapshot);
                % Gated data
                bt = tmp - repmat(mean(tmp,1),length(g(1):g(2)),1);
                % Standard deviation
                st = std(bt,0,1);
                % Normalized cross correlation
                r = (1/(k-1))*bt'*bt./(st'*st);
            else
                error('Aperture mode unknown')
            end
            r(isnan(r))=0;
            for d = 1:maxlag
                % Get mean cross correlations as a function of lag
                c = diag(r,d-1);
                cctmp(samp,d) = mean(c);
            end
        end
        cc(:,:,snapshot) = cctmp;
    end
end