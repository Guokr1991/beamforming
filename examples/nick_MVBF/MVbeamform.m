% MVbeamform.m - Minimum variance beamform
%
% Nick Bottenus - 3/20/2012
% BME 265 - Lab 6
%
% Parameters:
%   rf_in - RF channel data (focused)
%   subsize - Spatial averaging data subsize
%
% Returns: 
%   rf - Summed RF data
%   W  - Weighting vector used (ch x tr x freq)

function [rf W] = MVbeamform(rf_in,subsize)

%Preallocate
rf=zeros(size(rf_in,1),size(rf_in,3));
W=zeros(subsize,size(rf_in,3),size(rf_in,1));
%Loop through transmits (A-lines)
for tr=1:size(rf_in,3)
    %Get FT of single transmit data
    rf_ft=fft(rf_in(:,:,tr));
    rf_line=zeros(size(rf_ft,1),1);
    
    %Loop through frequency bands
    for freq=1:size(rf_ft,1)
        %Calculate spatial-averaged cov matrix
        num=size(rf_ft,2)-subsize+1;
        R=zeros(subsize,subsize);
        Yavg=zeros(subsize,1);
        for ch=1:num
            Y=rf_ft(freq,ch:ch+subsize-1);
            Y=Y(:);
            R=R+Y*Y';
            Yavg=Y+Yavg;
        end
        R=R/num;
        Yavg=Yavg/num;
        
        %Calculate W
        e=ones(length(R),1);
        W(:,tr,freq)=(R\e)/(e'*(R\e));
        
        %Apply to data
        rf_line(freq)=W(:,tr,freq)'*Yavg;
    end
    
    %IFT to get back A-line
    rf(:,tr)=ifft(rf_line);
end