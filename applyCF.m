function [rfcf cfmap] = applyCF(rfsteer, rfsum, zsteer, zsum)
% rf_cf = applyCF(rf, rf_out, z, zsum)
% 
% inputs:
% rfsteer - steered data pre-summation
% rfsum - focused and summed rf data
% zsteer - depths corresponding to pre-beamformed data
% zsum - depths corresponding to focused and summed rf data
% 
% outputs:
% rfcf - delayed and summed rf data with applied coherence factor (CF)

M = size(rfsteer,2);
rfcf = zeros(size(rfsum));
cfmap = zeros(size(rfsum));
for l = 1:size(rfsum,2) 
    for ii = 1:length(zsum)
        zi = find(zsum(ii) == zsteer);
        CF(ii) = abs(sum(rfsteer(zi,:,l),2)).^2./(M.*sum(abs(rfsteer(zi,:,l)).^2,2));
    end
    rfcf(:,l) = CF'.*rfsum(:,l);
    cfmap(:,l) = CF';
    clear CF
end
