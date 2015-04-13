close all; clear all; clc

load ./ece582_dataset/SSA_points.mat
fprintf('Data loaded. \n');
for l = 1:size(rf,3)
    imagesc(squeeze(rf(:,:,l))), colormap gray; colorbar; drawnow
    title(sprintf('A-line %d/%d',l,size(rf,3)))
end
    