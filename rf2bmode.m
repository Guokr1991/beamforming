function [env_win, rf_win, x_win, z_win] = rf2bmode(rf_out, db, x, z, x_range, z_range);
%[env_win, rf_win, x_win, z_win] = rf2bmode(rf_out, x, z)

if nargin == 1 || isempty(db)
    db = 40;
    disp('40 dB scale.')
end
if nargin <= 2
    x_range = [1 size(rf_out,2)];
    z_range = [1 size(rf_out,1)];
    x = 1:size(rf_out,2);
    z = 1:size(rf_out,1);
elseif nargin == 4
    x_range = [min(x) max(x)];
    z_range = [min(z) max(z)];
    disp('Use full data.')
elseif nargin == 5
    z_range = [min(z) max(z)];
    disp('Use full depth data.')
elseif nargin == 6
    if isempty(x_range)
        x_range = [min(x) max(x)];
    end
    if isempty(z_range)
        z_range = [min(z) max(z)];
    end
else
    
end
z_idx = find(z >= z_range(1) & z <= z_range(2));
x_idx = find(x >= x_range(1) & x <= x_range(2));

z_win = z(z_idx);
x_win = x(x_idx);

rf_win = rf_out(z_idx, x_idx);
env=abs(hilbert(rf_win));
env_win=20*log10(env/max(env(:)));

if nargin <= 2
    imagesc(env_win,[-db 0]); colormap('gray');
else
    imagesc(x_win, z_win, env_win, [-db 0]); colormap('gray'); axis image;
end