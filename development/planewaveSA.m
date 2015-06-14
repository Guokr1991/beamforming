clear all; close all; clc;
load point_target.mat

num_tx = size(rf,3);
num_el = length(acq_params.rx_pos);

fs = acq_params.fs;
c = acq_params.c;
t0 = acq_params.t0; % start sample

rref = (t0:t0+size(rf,1)-1).*c/(2*fs);

z = rref; % reference z values use first tx 
x = acq_params.rx_pos;

sample_origin = zeros(num_el,3); % matrix containing all rcv element positions [x,y,z];
sample_origin(:,1) = acq_params.rx_pos;

z_origin = unique(sample_origin(:,3));

nxgrid = 100; % grid sampling
nzgrid = 100;

xrange = [-0.01 0.01]; % min and max of x for full image
zrange = [0.029 0.032];

if nxgrid == 1
    xrange = [0 0];
end

xaxis = linspace(xrange(1),xrange(2),nxgrid);
zaxis = linspace(zrange(1),zrange(2),nzgrid);

[x_img, z_img] = meshgrid(xaxis,zaxis);

% IQdata = [t,n_el,angles]
h = zeros(size(z_img,1),num_tx,size(x_img,2));

for k = 1:num_tx
    RFdata = squeeze(rf(:,:,k));
    angle = acq_params.angle(k);
    dz = abs(z_img-repmat(z_origin,size(z_img)));
    rtx = dz.*cosd(angle);
    
    for x_id = 1:size(x_img,2)
        for z_id = 1:size(z_img,1)
            dz = abs(z_img(z_id,x_id)-z_origin); % dz to receive element
            rtx = dz.*cosd(angle);

            for el_id = 1:num_el
                dx = abs(x_img(z_id,x_id)-sample_origin(el_id,1)); % dx to receive element
                rrcv = sqrt(dz^2+dx^2);
                ri = (rtx+rrcv)/2; % half the distance for 2 way propagation

                yi(el_id) = interp1(rref, RFdata(:,el_id),ri);
                if isnan(yi(el_id))
                    yi(el_id) = 0;
                end
            end
            h(z_id,k,x_id) = sum(yi); clear yi;
        end
        disp('beamforming...')
    end
    fprintf('emission %.d processed. \n',k)
end

rf_out = squeeze(h);
% figure
% imagesc(1:num_tx, zaxis,h)


% plot(zaxis,rf_out)

% rf_sa = sum

% sum every three to generate high resolution images


