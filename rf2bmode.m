function [env, rf_win, x_win, z_win] = rf2bmode(rf_out, dbrange, x, z, xlims, zlims, flag);
%[env_win, rf_win, x_win, z_win] = rf2bmode(rf_out, 40, x, z)
%
% rf_out - beamformer output of focused radiofrequency signal array
% dbrange - dynamic range of output image in dB
% x - lateral positions to envelope detect [m]
% z - axial depth positions to envelope detect [m]
% xlims - lateral limits of output image [m]
% zlims - axial limits of output image [m]
% flag - '0' output only, '1' output and show image in current figure handle

xlab = [];
ylab = [];
if nargin < 2 || isempty(dbrange)
    dbrange = 40;
end
if nargin < 3 || isempty(x)
    xlims = [1 size(rf_out,2)];
    x = 1:size(rf_out,2);
    xlab = 'jth index';
end
if nargin < 4 || isempty(z)
    zlims = [1 size(rf_out,1)];
    z = 1:size(rf_out,1);
    ylab = 'ith index';
end
if nargin < 5 || isempty(xlims)
    xlims = [min(x) max(x)];
end
if nargin < 6 || isempty(zlims)
    zlims = [min(z) max(z)];
end
if nargin < 7 || isempty(flag)
    flag = 1;
end

z_idx = find(z >= zlims(1) & z <= zlims(2));
x_idx = find(x >= xlims(1) & x <= xlims(2));

z_win = z(z_idx);
x_win = x(x_idx);

rf_win = rf_out(z_idx, x_idx);
env=abs(hilbert(rf_win));
env_win=20*log10(env/max(env(:)));
xtic = x_win;
ytic = z_win;

% display detected ultrasound image in dB scale given beamformed RF data
switch flag
    case 0
    case 1
        if isempty(xlab) && isempty(ylab), atype = 'image'; 
        else, atype = 'square'; end
        if isempty(ylab),ylab = 'Axial Distance [mm]'; ytic = 1000*z_win; end
        if isempty(xlab),xlab = 'Lateral Distance [mm]'; xtic = 1000*x_win; end
        imagesc(xtic, ytic, env_win, [-dbrange 0]); colormap('gray'); 
        xlabel(xlab), ylabel(ylab), axis(atype);
        
end