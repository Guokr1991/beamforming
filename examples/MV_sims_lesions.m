%% MV simulations in FIELD II
%
% by Marko Jakovljevic
% 
%last updated on: 11/28/2011
%
clear all
addpath(genpath('/data/mj66/code/Field_II_linux.tar.gz_FILES/'));
addpath(genpath('/data/mj66/code/usp/'));
addpath(genpath('/nefs/jjd/simulation/VF7-3/contrast/simdata/'));
cd /home/mj66/Classes/Intro_to_Ultrasound/
field_init(-1)


%%
% set parameters from Gran's paper
% f0=7e6;
run /nefs/jjd/simulation/VF7-3/contrast/code/param
%fsample=160e6;
%c=1540;
% width=0.075e-3;
% element_height=6/1000;
% kerf=0.035e-3;
% focal=[0 0 45]/1000;
N_elements=96;
% fbwd = 0.6;
% set_field('c',c);
% set_field('fs',fs);
% 
% ap_size = N_elements*(kerf+width);
% pitch = kerf+width;

%% define field size

% imgszax = 2/100; % image size in axial direction
% imgszlat = ap_size + 2/1000; % image size in lateral direction
% axspacing=c/2/fs;
% latspacing=axspacing;
% 
% aximg = (-imgszax/2:axspacing:imgszax/2)+focal(3);
% latimg=-imgszlat/2:latspacing:imgszlat/2;

%%

% sector=8/1000;        %  Size of image sector; INCREASE LATER
% z_focus=45/1000; % ASK FOR THE DIFFERENCE BETWEEN THIS AND FOCUS ABOVE
% lambda = c/f0;
% numbeams = 80;

%  Pre-allocate some storage
%image_data=zeros(800,no_lines);
%sch_data = zeros(2500,N_elements,no_lines);
sch_data = [];
%start_samp = zeros(1,no_lines);
%lat = zeros(1,no_lines);

for ii=1:8

    load(sprintf('/nefs/jjd/simulation/VF7-3/contrast/data/rfdata1_2_%d.mat',ii));
    %sch_data(1:size(signal,1),:,ii) = signal;
    sch_data = cat(3,sch_data,bmode.signal);
%     start_samp(ii) = start_sample;
%     lat(ii) = latshift;
    disp(['Loading sch data block ' num2str(ii) ' ...'])
    
end

numSamples = size(sch_data,1);

%% make the B-mode
sd_rf_data = squeeze(sum(sch_data,2)); % sum-and-delay rf data

% min_sample=min(start_samp);
% max_sample=max(start_samp);
% [n,m]=size(sd_rf_data);
% store=zeros(n+max_sample-min_sample,m);
% for k=1:size(sd_rf_data,2),
%     v=[zeros(start_samp(k)-min_sample,1); sd_rf_data(1:n,k)];
%     store(1:max(size(v)),k)=v;
% end;
% sd_rf_data = store;
% st=min(start_samp)*(1/fs);

%[sd_rf_data st] = ShiftStartTimes(s_time, sd_rf_data, fs);

st = bmode.start_sample*(1/fsample);
lat = ((1:1:numbeams) - numbeams/2)* beamspacing;

env = abs(hilbert(sd_rf_data));
env = env/max(env(:));

figure;
s1 = subplot(1,1,1);
t = st:1/fsample:(st+size(env,1)*(1/fsample));
axial = t*c/2;
imagesc(lat,axial,db(env),[-50 0])
colormap(gray)
colorbar
%axis normal
xlabel({'Lateral Distance (mm)'},'FontSize',11)
ylabel('Axial Distance (mm)','FontSize',11)
ylim([35.5 41]*10^-3)
%title('B-mode','FontSize',11)

set(s1, 'FontSize', 11)
set(s1, 'LineWidth', 2)
set(s1, 'XColor', [0.25 0.25 0.25])
set(s1, 'YColor', [0.25 0.25 0.25])
set(s1, 'ZColor', [0.25 0.25 0.25])

% set(s1,'XTick',[-8 -4 0 4 8]*10^-3)
% set(s1,'XTickLabel',[-8 -4 0 4 8])
set(s1,'XTick',[-4 -2 0 2 4]*10^-3)
set(s1,'XTickLabel',[-4 -2 0 2 4])
set(s1,'YTick',[36 37 38 39 40 41]*10^-3)
set(s1,'YTickLabel',[36 37 38 39 40 41])

% figure;
% s2 = subplot(1,1,1);
% plot(lat,db(env(640,:)),'LineWidth',1.5)
% ylim([-80 0])
% xlim([min(lat) max(lat)])
% xlabel({'Lateral Distance (mm)'},'FontSize',11)
% ylabel('Power (dB)','FontSize',11)
% % set(s2,'XTick',[-8 -4 0 4 8]*10^-3)
% % set(s2,'XTickLabel',[-8 -4 0 4 8])
% set(s2,'XTick',[-4 -2 0 2 4]*10^-3)
% set(s2,'XTickLabel',[-4 -2 0 2 4])

%% DFT, estimating R, and IDFT
% add white noise to the data
sch_data = sch_data(400:1700,:,:); % cut off the top and bottom
noise_power = 1/(10^6)*sum(sch_data.^2,1); % hardcoded for 60 dB
white_noise = randn(size(sch_data)).*repmat(noise_power,[size(sch_data,1),1,1]);
sch_data_noise = sch_data + white_noise;

num_samples = 50; % number of time samples to take
num_freq = 256; % number of descrete frequencies
L = floor(N_elements/4); % number of elements in the subaperture

steer = ones(L,1); % steering vector
discrete_times = [0:1:num_samples-1];
discrete_times = repmat(discrete_times,num_freq,1);
%freqs = repmat([linspace(0,1-1/num_freq,num_freq)]',1,num_samples);
freqs = repmat([linspace(0,0.2,num_freq)]',1,num_samples);
exponents = freqs .* discrete_times;

dft_matrix = exp(-1i*2*pi*exponents);
idft_matrix = exp(1i*2*pi*exponents');

mv_data_dft = zeros(num_freq,numbeams);
mv_data_idft = zeros(size(sch_data,1),numbeams);

%for ii=1
for ii=1:floor(size(sch_data,1)/num_samples)
    %for jj = 93
    for jj=1:numbeams
        disp(['Processing A-line ' num2str(jj) ' ...'])
        
        temp2 = dft_matrix * sch_data_noise(((ii-1)*num_samples)+1:ii*num_samples,:,jj); % take a DFT
        %temp2 = dft_matrix * sch_data_noise(900:1099,:,jj); % take a DFT
        
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
    mv_data_idft(((ii-1)*num_samples)+1:ii*num_samples,:) = 1/num_freq * idft_matrix * mv_data_dft;
    
end

%% make MVDR image

% [mv_data_shifted st] = ShiftStartTimes(start_times, mv_data_idft, fs);
% min_sample=min(start_samp);
% max_sample=max(start_samp);
% [n,m]=size(mv_data_idft);
% store=zeros(n+max_sample-min_sample,m);
% for k=1:size(mv_data_idft,2),
%     v=[zeros(start_samp(k)-min_sample,1); mv_data_idft(1:n,k)];
%     store(1:max(size(v)),k)=v;
% end;
% mv_data_idft = store;
%st=min(start_samp)*(1/fs);
st=min(bmode.start_samp+400)*(1/fsample);

env2 = abs(mv_data_idft);
env2 = env2/max(env2(:));

figure;
s1=subplot(1,1,1);
t2 = st:1/fsample:(st+size(env2,1)*(1/fsample));
axial2 = t2*c/2;
imagesc(lat,axial2,db(env2),[-50 0])
colormap(gray)
colorbar
xlabel({'Lateral Distance (mm)'},'FontSize',11)
ylabel('Axial Distance (mm)','FontSize',11)
ylim([35.5 41]*10^-3)

set(s1, 'FontSize', 11)
set(s1, 'LineWidth', 2)
set(s1, 'XColor', [0.25 0.25 0.25])
set(s1, 'YColor', [0.25 0.25 0.25])
set(s1, 'ZColor', [0.25 0.25 0.25])

% set(s1,'XTick',[-8 -4 0 4 8]*10^-3)
% set(s1,'XTickLabel',[-8 -4 0 4 8])
set(s1,'XTick',[-4 -2 0 2 4]*10^-3)
set(s1,'XTickLabel',[-4 -2 0 2 4])
set(s1,'YTick',[36 37 38 39 40 41]*10^-3)
set(s1,'YTickLabel',[36 37 38 39 40 41])

% figure;
% s2 = subplot(1,1,1);
% plot(lat,db(env2(73,:)),'b','LineWidth',1.5);
% ylim([-60 0])
% xlim([min(lat) max(lat)])
% xlabel({'Lateral Distance (mm)'},'FontSize',11)
% ylabel('Power (dB)','FontSize',11)
% set(s2,'XTick',[-4 -2 0 2 4]*10^-3)
% set(s2,'XTickLabel',[-4 -2 0 2 4])

%% raw rf data figures (just for the project)

figure;
s1=subplot(1,1,1);
% t = [st:1/fs:(st+size(env2,1)*(1/fs))]-40/1000/c*2;
t = st:1/fsample:(st+size(env2,1)*(1/fsample));
imagesc([1:1:N_elements],t,sch_data(:,:,1))
%imagesc([1:1:N_elements],t,bmode.signal(:,:,1))
colormap(gray)
xlabel({'Element Number'},'FontSize',14)
ylabel('Time (s)','FontSize',14)
set(s1, 'FontSize', 14)
set(s1, 'LineWidth', 2)
set(s1, 'XColor', [0.25 0.25 0.25])
set(s1, 'YColor', [0.25 0.25 0.25])
set(s1, 'ZColor', [0.25 0.25 0.25])





