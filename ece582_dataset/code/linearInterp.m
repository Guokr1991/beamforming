function yi = linearInterp(x, y, xi)
% Linear interp for focal steering - Will Long. Latest revision: 4/2/15
% xi - matrix of x values to be interpolated, columns of xi will be 
%   interpolated based on reference data in corresponding column of y
% x - vector of independent data series with known y values
% y - matrix of known dependent data (function of x), datasets are
%   arranged column-wise
% yi - the column by column interpolated values based on x, xi, y inputs 

if size(y,2) == 1
    y = repmat(y,[1 size(xi,2)]);
end
    
dx = diff([x(1) x(2)]);

xi(xi>=max(x)) = NaN;
yi = zeros(size(xi,1),size(xi,2));
for j = 1:size(xi,2)
    m = [diff(y(:,j))./dx; NaN];
    nn_i = floor((xi(:,j)-x(1))./dx)+1;
    nn_i(isnan(nn_i)) = 1;
    yi(:,j) = m(nn_i).*(xi(:,j)-x(nn_i))+y(nn_i,j);
end
