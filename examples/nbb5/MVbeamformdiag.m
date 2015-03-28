% MVbeamformdiag.m - Minimum variance beamform - diagonal loading
%
% Nick Bottenus - 3/20/2012
% BME 265 - Lab 6
%
% Parameters:
%   rf_in - RF channel data (focused)
%
% Returns: 
%   rf - Summed RF data
%   W  - Weighting vector used (ch x tr x freq)
function [rf W] = MVbeamformdiag(rf_in)

%Preallocate
rf=zeros(size(rf_in,1),size(rf_in,3));
W=zeros(size(rf_in,2),size(rf_in,3),size(rf_in,1));
%Loop through transmits (A-lines)
for tr=1:size(rf_in,3)
    %Get FT of single transmit data
    rf_ft=fft(rf_in(:,:,tr));
    rf_line=zeros(size(rf_ft,1),1);
    
    %Loop through frequency bands
    for freq=1:size(rf_ft,1)
        %Calculate cov matrix
        Y=rf_ft(freq,:);
        Y=Y(:);
        R=Y*Y';
        
        %Diagonal loading
        loading=mean(R(:));
        R=R+eye(size(R))*loading;
        
        %Calculate W
        e=ones(length(R),1);
        W(:,tr,freq)=(R\e)/(e'*(R\e));
        
        %Apply to data
        rf_line(freq)=W(:,tr,freq)'*Y;
    end
    
    %IFT to get back A-line
    rf(:,tr)=ifft(rf_line);
end