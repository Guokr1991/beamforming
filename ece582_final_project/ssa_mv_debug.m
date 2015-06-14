%% data for SSA
close all; clear all; clc

figure(1)
load ./ece582_dataset/SSA_points.mat
rf_out = squeeze(sum(rf,2));
fprintf('Data loaded. \n');
for l = 1:size(rf,3)
    subplot(121)
    imagesc(squeeze(rf(:,:,l))), colormap gray; 
    title(sprintf('A-line %d/%d',l,size(rf,3)))
    
    subplot(122)
    imagesc(rf_out(:,1:l))
    drawnow
    pause(0.1)
end

% %%
% figure(1)
% load ./ece582_dataset/SSA_points.mat
% fprintf('Data loaded. \n');
% 
% 
% for l = 1:size(rf,2)
% 
%     imagesc(squeeze(rf(:,l,:))), colormap gray; colorbar; 
%     title(sprintf('Sub-image %d/%d',l,size(rf,2)))  
%     drawnow
% end
    
%% data for simulation
close all; clear all; clc

figure(1)
load ./ece582_dataset/point_grid_MVdebug.mat
rf_out = squeeze(sum(rf_steer,2));
fprintf('Data loaded. \n');
for l = 1:size(rf_steer,3)
    subplot(121)
    imagesc(squeeze(rf_steer(:,:,l))), colormap gray;
    title(sprintf('A-line %d/%d',l,size(rf_steer,3)))
    
    subplot(122)
    imagesc(rf_out(:,1:l))
    drawnow
    pause(0.1)
end

% %%
% figure(1)
% load ./ece582_dataset/SSA_points.mat
% fprintf('Data loaded. \n');
% 
% 
% for l = 1:size(rf,2)
% 
%     imagesc(squeeze(rf(:,l,:))), colormap gray; colorbar; 
%     title(sprintf('Sub-image %d/%d',l,size(rf,2)))  
%     drawnow
% end
    
%%
clear all; close all;

load ./ece582_dataset/SSA_points.mat
rf_con = squeeze(sum(rf,2));
figure(1)
subplot(211); rf2bmode(rf_con,50,x,z);
figure(2)
subplot(211); imagesc(x,z,rf_con); colorbar, colormap gray

load ssapoints_MV128fft.mat
rf_MV = rf_out;
figure(1)
subplot(212); rf2bmode(rf_out,50,x,z);
figure(2)
subplot(212); imagesc(x,z,rf_out); colorbar, colormap gray