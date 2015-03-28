% MVbeamformarray.m - Minimum variance beamform array response
%
% Nick Bottenus - 3/20/2012
% BME 265 - Lab 6
%
% Parameters:
%   rf_in - RF channel data (focused)
%   rf_pt - RF channel data for single point (focused)
%
% Returns: 
%   rf - Summed RF data
%   W  - Weighting vector used (ch x freq)

function [rf W] = MVbeamformarray(rf_in,rf_pt)

%Preallocate
rf=zeros(size(rf_in,1),size(rf_in,3));
subsize=32;
W=zeros(subsize,size(rf_in,1));

%Calculate weight vector down center line
rf_ft=fft(rf_in(:,:,(size(rf_in,3)+1)/2));
%Loop through frequency bands
for freq=1:size(rf_ft,1)
    %Calculate spatial-averaged cov matrix
    num=size(rf_ft,2)-subsize+1;
    R=zeros(subsize,subsize);
    for ch=1:num
        Y=rf_ft(freq,ch:ch+subsize-1);
        Y=Y(:);
        R=R+Y*Y';
    end
    R=R/num;

    %Calculate W
    e=ones(length(R),1);
    W(:,freq)=(R\e)/(e'*(R\e));
end
    
%Loop through transmits (A-lines)
for tr=1:size(rf_pt,3)
    %Get FT of single transmit data
    rf_ft=fft(rf_pt(:,:,tr));
    rf_line=zeros(size(rf_ft,1),1);
    
    %Loop through frequency bands
    for freq=1:size(rf_ft,1)
        Yavg=zeros(subsize,1);
        for ch=1:num
            Y=rf_ft(freq,ch:ch+subsize-1);
            Yavg=Y(:)+Yavg;
        end
        Yavg=Yavg/num;
        
        %Apply to data
        rf_line(freq)=W(:,freq)'*Yavg;
    end
    
    %IFT to get back A-line
    rf(:,tr)=ifft(rf_line);
end