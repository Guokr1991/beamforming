%% MV simulations in FIELD II
%
% by Marko Jakovljevic
% 
%last updated on: 11/28/2011
%
clear all
addpath(genpath('/data/mj66/code/Field_II_linux.tar.gz_FILES/'))
cd /home/mj66/Classes/Intro_to_Ultrasound/
field_init(-1)


%%
% set parameters from Gran's paper
f0=7e6;
fs=100e6;
c=1540;
width=0.075e-3;
element_height=6/1000;
kerf=0.035e-3;
focal=[0 0 45]/1000;
N_elements=129;
fbwd = 0.6;
set_field('c',c);
set_field('fs',fs);

ap_size = N_elements*(kerf+width);
pitch = kerf+width;

%%   Generate transmit aperture
%transmit_aperture = xdc_linear_array (N_elements, width, element_height, kerf, 5, 5,focal);
transmit_aperture = xdc_linear_array (1, width, element_height, kerf, 5, 5,[0 0 100]); % single element transmission
receive_aperture_sch = xdc_linear_array (N_elements, width, element_height, kerf, 5, 5,[0 0 40/1000]); % for sch

impulse_response = makeImpulseResponse(fbwd,f0,fs);
xdc_impulse (transmit_aperture, impulse_response); % define impusle response of the emit aperture
xdc_impulse (receive_aperture_sch, impulse_response); % define impusle response of the receive aperture

excitation=sin(2*pi*f0*(0:1/fs:2/f0));
xdc_excitation (transmit_aperture, excitation);

% calc time shift to account for difference between FIELD and usp_array
tshift = size(conv(conv(excitation,impulse_response),impulse_response),2)/2/fs;
%% define field size

imgszax = 2/100; % image size in axial direction
imgszlat = ap_size + 2/1000; % image size in lateral direction
axspacing=c/2/fs;
latspacing=axspacing;

aximg = (-imgszax/2:axspacing:imgszax/2)+focal(3);
latimg=-imgszlat/2:latspacing:imgszlat/2;


%%  create point scatterers

% r = 0.035/1000; % radius of a single target is 0.035 mm
% t_spacing = 1e-4; % target spacing of 0.1 mm
% num_target = 3; % number of circle targets (make it odd)

%ax_pos = focal(3); % axial positions of targets (center)
%lat_pos = linspace(-(num_target-1)/2,(num_target-1)/2,num_target)*t_spacing;
ax_pos = 40/1000;
lat_pos = 0;
% test = zeros(size(aximg,2),size(latimg,2));
% x_points = [];
% z_points = [];

x_points = [-5 0 5]/1000;
z_points = ax_pos*ones(size(x_points));
%for ii=1:num_target
% for ii=1    
%     for jj=1000:1500 % for axial; hardcode to cut down the time
%         for kk=700:1500 % for lateral; hardcode to cut down the time
%             if (aximg(jj)-ax_pos)^2+(latimg(kk)-lat_pos(ii))^2 <= r^2
%                 z_points = cat(2,z_points,aximg(jj));
%                 x_points = cat(2,x_points,latimg(kk));
%                 test(jj,kk) = 1;
%             end
%         end
%     end
% end

y_points = zeros(size(z_points));
points = [x_points' y_points' z_points'];
amplitudes = ones(size(points,1),1);
%%

%no_lines=40;            %  Number of A-lines in image; INCREASE LATER
sector=16/1000;        %  Size of image sector; INCREASE LATER
z_focus=45/1000; % ASK FOR THE DIFFERENCE BETWEEN THIS AND FOCUS ABOVE
lambda = c/f0;
no_lines = 8*ceil(sector/(lambda*z_focus/ap_size))+1; % ensures its odd number

%  Pre-allocate some storage
image_data=zeros(800,no_lines);
sch_data = zeros(800,N_elements,no_lines);
start_times = zeros(1,no_lines);

x = linspace(-sector/2,sector/2,no_lines);

for ii=1:no_lines
% %for ii=20
%     xdc_center_focus (transmit_aperture, [x(ii) 0 0]);
%     xdc_focus (transmit_aperture, 0, [x(ii) 0 z_focus]);
    
    x_points_shift = x_points + x(ii); % shift the scatterers
    points = [x_points_shift' y_points' z_points'];
    %xdc_center_focus (receive_aperture_sch, [0 0 0]);
    % xdc_focus (receive_aperture_sch, 0, [x(ii) 0 z_focus]);
    %xdc_dynamic_focus (receive_aperture_sch,0,0,0);
    
    % getting single-channel data
    xdc_times_focus(receive_aperture_sch,0,zeros(1,N_elements)) % set delays to 0 to get unfocused data
    [temp t2]=calc_scat_multi(transmit_aperture, receive_aperture_sch, points, amplitudes);
    sch_data(1:size(temp,1),:,ii) = temp;
    start_times(ii) = t2;
    

end

%% re-focus the data
numSamples = size(sch_data,1);

% Create defocusing object
% defoc = usp_focus;
% defoc.type = 'fixed';
% defoc.origin = [0; 0; 0]';
% defoc.direction = [0; 0; 1]';
% defoc.range = 40/1000;

% Create focus object
foc = usp_focus;
foc.type = 'dynamic';
%foc.type = 'fixed';
foc.origin = [zeros(1,no_lines); ...
    zeros(1,no_lines); zeros(1,no_lines)]';
foc.direction = [zeros(1,no_lines); zeros(1,no_lines); ...
    ones(1,no_lines)]';
% foc.range = z_focus; % Tx focus
foc.range = 40/1000; % Tx focus
foc.sound_speed = c;

% Sort the data; bmode actually contains focused single channel data
bmode = usp_array;
bmode.signal = sch_data - repmat(mean(sch_data,1),[numSamples,1,1]);
start_sample = round((min(start_times) - tshift)*fs);
bmode.start_sample = start_sample;
bmode.sampling_frequency = fs;
% bmode.position = [-pitch*(N_elements/2-1):pitch:pitch*(N_elements/2); ...
%     zeros(1,N_elements); zeros(1,N_elements)]';
bmode.position = [-pitch*64:pitch:pitch*64; ... % HARD-CODED
    zeros(1,N_elements); zeros(1,N_elements)]';
bmode.index = [1:N_elements]';
bmode = focus(bmode,foc);

%% make the B-mode
sd_rf_data = squeeze(sum(bmode.signal,2)); % sum-and-delay rf data
s_time = start_times-tshift;
[sd_rf_data st] = ShiftStartTimes(s_time, sd_rf_data, fs);

env = abs(hilbert(sd_rf_data));
env = env/max(env(:));

figure;
lat = x;
t = st:1/fs:(st+size(env,1)*(1/fs));
axial = t*c/2;
imagesc(lat,axial,db(env),[-50 0])
colormap(gray)
figure; plot(lat,db(env(73,:)))
ylim([-60 0])
xlim([min(lat) max(lat)])
%% DFT, estimating R, and IDFT
% add white noise to the data
noise_power = 1/(10^6)*sum(bmode.signal.^2,1); % hardcoded for 60 dB
white_noise = randn(size(bmode.signal)).*repmat(noise_power,[size(bmode.signal,1),1,1]);
sch_data_noise = bmode.signal + white_noise;

num_samples = numSamples; % number of time samples to take
num_freq = 256; % number of descrete frequencies
L = floor(N_elements/4); % number of elements in the subaperture

steer = ones(L,1); % steering vector
discrete_times = [0:1:num_samples-1];
discrete_times = repmat(discrete_times,num_freq,1);
%freqs = repmat([linspace(0,1-1/num_freq,num_freq)]',1,num_samples);
freqs = repmat([linspace(0,0.22,num_freq)]',1,num_samples);
exponents = freqs .* discrete_times;

dft_matrix = exp(-1i*2*pi*exponents);
idft_matrix = exp(1i*2*pi*exponents');

mv_data_dft = zeros(num_freq,no_lines);
for ii=1
%for ii=1:floor(size(sch_data,1)/num_samples)
    %for jj = 20
    for jj=1:no_lines
        disp(['Processing A-line ' num2str(jj) ' ...'])
        
        temp2 = dft_matrix * sch_data_noise(((ii-1)*num_samples)+1:ii*num_samples,:,jj); % take a DFT
        
        for omega = 1:num_freq
            R = zeros(L,L); % reset the covariance matrix
            G_av = zeros(L,1);
            for p = 1:N_elements-L+1
                G = temp2(omega,p:p+L-1).';
                R = G*G' + R; % estimate R
                G_av = G+G_av; % average G
            end
            R = 1/p*R;
            G_av = 1/p*G_av;
            weights = (inv(R)*steer)/(steer'*inv(R)*steer); % calculate the weights
            out_temp = weights' * G_av; % beamform in the freq domain    
            mv_data_dft(omega,jj) = out_temp; % store result
        end
        
    end
    
    % inverse DFT to give broadband output
    mv_data_idft = 1/num_freq * idft_matrix * mv_data_dft;
    
end

%% make MVDR image

% [mv_data_shifted st] = ShiftStartTimes(start_times, mv_data_idft, fs);
st = min(start_times);
env2 = abs(mv_data_idft);
env2 = env2/max(env2(:));

figure;
lat = x;
t = st:1/fs:(st+size(env2,1)*(1/fs));
axial = t*c/2;
imagesc(lat,axial,db(env2),[-50 0])
colormap(gray)

figure; plot(lat,db(env2(73,:)));










