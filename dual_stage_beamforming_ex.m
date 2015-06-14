clear all; close all; clc
cd ~/Google' Drive'/Trahey' Lab'/DTU_summer_school/exercises/exercise_4/

load Ex1.mat

x = squeeze(xmt.sample_origin(:,1));
z0 = unique(xmt.sample_origin(:,3));
zf = unique(xmt.focus_point(:,3));
zr = unique(xmt.scanline_ref_point(:,3));
z = z0:xmt.dr:z0+xmt.dr*(size(RFdata,1)-1);
env = abs(hilbert(RFdata));
imagesc(x,z,20.*log10(env/max(env(:))),[-60 0]); colormap gray; axis image

nxgrid = 5;
xaxis = linspace(min(x),max(x),nxgrid);
nzgrid = 401;
zaxis = linspace(min(z),max(z),nzgrid);

[x_img,z_img] = meshgrid(xaxis, zaxis);


IQdata = hilbert(RFdata);
f_num = 2;
c = 1540;

h = zeros(size(z_img,1),xmt.no_lines,size(x_img,2));
% Loop over all scan lines
for k = 1:xmt.no_lines
    for x_id = 1:size(x_img,2) % loop in lateral direction
        for z_id = 1:size(z_img,1) % loop in axial direction
            
            % determine axial distance from image point to virtual source
            dz = abs(z_img(z_id,x_id)-xmt.focus_point(k,3));
            % determine lateral distance from image point to scan line
            dx = abs(x_img(z_id,x_id)-xmt.focus_point(k,1));
            % determine apodization value
            delta = dz/f_num;
            n = dx/delta;
            if(n <= 0.5)
                W = 0.54+0.46*cos(2*n*pi);
            else
                W = 0;
            end
            
            % determine time of flight (tof)
            if (z_img(z_id,x_id) >= xmt.focus_point(k,3)) 
                tof = 2/c*(abs(zf-zr)+sqrt(dz^2+dx^2));
            else
                tof = 2/c*(abs(zf-zr)-sqrt(dz^2+dx^2));
            end
            zi = tof*c/2;
            % interpolate with tof into the IQdata1
            yi = interp1(z, IQdata(:,k),zi);
            if(isnan(yi))
                yi = 0;
            end
            
            h(z_id,k,x_id) = yi*W;
            
        end
    end
    fprintf('Line %.d processed. \n',k)
end

figure
env = abs(squeeze(sum(h,2)));
imagesc(xaxis, zaxis, 20.*log10(env/max(env(:))),[-60 0]); colormap gray; axis image

% for i = 1:length(X)
%     for j = 1:length(Z)
%         if Z(i) > 
%         elseif 
%     end
% end
% Zf = repmat(zf,size(Z));
% Zr = repmat(zr,size(Z));
% 
% abs(Zf-Zr)