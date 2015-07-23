function yi = linearInterp(x, y, xi)
% xi - matrix of x values to be interpolated, columns of xi will be 
%   interpolated based on reference data in corresponding column of y
% x - vector of independent data series with known y values
% y - matrix of known dependent data (function of x), datasets are
%   arranged column-wise
% yi - the column by column interpolated values based on x, xi, y inputs 

if size(y,2) == 1
    y = repmat(y,[1 size(xi,2)]);
end
    
tmp = unique(diff(x));
dx = tmp(1);

xi(xi>=max(x)) = NaN;
yi = zeros(size(xi));

if size(x,1) == 1
    x = x';
end

for jj = 1:size(xi,2)
    m = [diff(y(:,jj))./dx; 0];
    nn_i = floor((xi(:,jj)-x(1))./dx);
    yi(:,jj) = m(nn_i).*(xi(:,jj)-x(nn_i))+y(nn_i,jj);
end
