% MVbeamformtemp.m - Minimum variance beamform - temporal averaging
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

function [rf W] = MVbeamformtemp(rf_in)

%Preallocate
rf=zeros(size(rf_in,1),size(rf_in,3));
W=zeros(size(rf_in,2),size(rf_in,3),size(rf_in,1));
rf_ft=zeros(size(rf_in,1),size(rf_in,2),size(rf_in,4));
%Loop through transmits (A-lines)
for tr=1:size(rf_in,3)
    for n=1:size(rf_in,4)
        %Get FT of single transmit data
        rf_ft(:,:,n)=fft(rf_in(:,:,tr,n));
    end
    rf_line=zeros(size(rf_in,1),1);
    
    %Loop through frequency bands
    for freq=1:size(rf_ft,1)
        %Calculate temporal-averaged cov matrix
        R=zeros(size(rf_ft,2));
        Yavg=zeros(size(rf_ft,2),1);
        for n=1:size(rf_ft,3)
            Y=rf_ft(freq,:,n);
            Y=Y(:);
            R=R+Y*Y';
            Yavg=Y+Yavg;
        end
        R=R/size(rf_ft,3);
        Yavg=Yavg/size(rf_ft,3);
        
        %Calculate W
        e=ones(length(R),1);
        W(:,tr,freq)=(R\e)/(e'*(R\e));
        
        %Apply to data
        rf_line(freq)=W(:,tr,freq)'*Yavg;
    end
    
    %IFT to get back A-line
    rf(:,tr)=ifft(rf_line);
end