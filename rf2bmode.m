function [env, rf_win, x_win, z_win] = rf2bmode(rf_out, dbrange, x, z, xlims, zlims, flag);
%[env_win, rf_win, x_win, z_win] = rf2bmode(rf_out, x, z)
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

switch flag
    case 0
    case 1
        imagesc(1000*x_win, 1000*z_win, env_win, [-dbrange 0]); colormap('gray'); 
        xlabel(xlab), ylabel(ylab)
        if isempty(ylab),ylabel('Axial Distance [mm]'); end
        if isempty(xlab),xlabel('Lateral Distance [mm]'); end
        if isempty(xlab) && isempty(ylab), axis image; end
end