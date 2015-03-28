% simdata.m - Simulate single point target at focus for SA data
% BME 265 - Lab 6
% Nick Bottenus
% 3/15/12
% 
% Parameters:
%   x_pts   - x-location of scatterers
%
% Returns:
%   rf      - Axial x Channel x Transmit
%   ax, lat - Axial and lateral vectors in meters

function [rf lat ax acq_params bf_params] = simscan(x_pts)

%Trandsucer parameters
f0=5e6;                     % Center frequency [Hz]
num_elements=128;           % Number of elements
c=1540;                     % Speed of sound [m/s]
BW=4e6/f0;                  % Fractional BW
fs=120e6;                   % Sampling frequency [Hz]
lambda=c/f0;                % Wavelength [m]
kerf=.01/1000;              % Kerf [m]
element_width=lambda/2-kerf;% Element width [m]
element_height=5/1000;      % Element height [m]
ax_focus=40/1000;           % Axial focus [m]
sub_x=1;                    % Number of subelements in x and y
sub_y=round(element_height/element_width);

%Set field parameters
set_field('c',c);
set_field('fs',fs);

%Set up impulse response
tc = gauspuls('cutoff',f0,BW,-6,-40);
t = -tc:1/fs:tc;
impulse_response = gauspuls(t,f0,BW);

%Set up excitation
excitation=1;

%Generate transmit aperture (single element)
emit_aperture = xdc_linear_array (num_elements, element_width, ...
    element_height, kerf, sub_x, sub_y, [0 0 ax_focus]);
xdc_impulse (emit_aperture, impulse_response);
xdc_excitation (emit_aperture, excitation);

%Generate receive aperture
receive_aperture = xdc_linear_array (num_elements, element_width,...
    element_height, kerf, sub_x, sub_y, [0 0 ax_focus]);
xdc_impulse (receive_aperture, impulse_response);

%Create phantom
y_pts = zeros(1,length(x_pts));
z_pts = ones(1,length(x_pts))*ax_focus;
amp = ones(length(x_pts),1);

%Set up scan parameters
num_lines=num_elements-1;
scan_radius=(num_elements/2-1)*element_width; %Radius of scan [m]
lat=linspace(-scan_radius,scan_radius,num_lines);


for i=1:num_lines
    
    %Simulate
    pts = [x_pts-lat(i);y_pts;z_pts]';
    [v1, t1]=calc_scat_multi(emit_aperture,receive_aperture,pts,amp);
    
    %Store the result
    image_data(1:size(v1,1),:,i)=v1;
    times(i) = t1;
    
end

%Adjust the data in time
min_sample=min(times)*fs;
max_sample=max(times)*fs;
n=size(image_data,1);
n=n+(max_sample-min_sample);
rf=zeros(n,num_elements,num_lines);
for i=1:num_lines
    for j=1:num_elements
        rf_temp=[zeros(times(i)*fs-min_sample,1); image_data(:,j,i)];
        rf(1:size(rf_temp,1),j,i)=rf_temp;
    end
end

%Create axis vectors (mm)
lat=lat*1000;
offset=length(conv(impulse_response,conv(impulse_response,excitation)))/(2*fs);
ax=(min(times):1/fs:(min(times)+(size(rf,1)-1)/fs))-offset;
ax=ax*c/2*1000;

%Free space for aperture


tdr_info = xdc_get(receive_aperture,'rect');
rx_pos = unique(tdr_info(24,:));
xdc_free (emit_aperture)
xdc_free (receive_aperture)

acq_params.c = c;
acq_params.fs = fs;
acq_params.t0 = min_sample;
acq_params.rx_pos = rx_pos;

bf_params.x = lat/1000;

