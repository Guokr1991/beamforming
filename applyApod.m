function [rf_out, rf_apo] = applyApod(rf,type,param)
% [rf_out, rf_apo] = applyApod(rf,'gauss',2.5)
% 
% rf - steered data pre-summation
% type - type of apodization to apply
% param - parameters for apodization (type specific)
% 
% rf_out - summed apodized rf data
% rf_apo - pre-summed apodized rf data
    
L = size(rf,2); 

switch type
    case 'gauss'
        win = gausswin(L,param)'; 
    case 'hann'

end

winmat = repmat(win,size(rf,1),1);

figure
plot(win)
title('Apodization weights');

rf_apo = zeros(rf);
for l = 1:size(rf,3); 
    rf_apo(:,:,l) = squeeze(rf(:,:,l)).*winmat; 
end
rf_out = squeeze(sum(rf_apo,2));
