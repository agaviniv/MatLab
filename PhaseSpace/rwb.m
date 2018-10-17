function vba = rwb(m)
%   By R. Sarthour


if nargin < 1, m = size(get(gcf,'colormap'),1); end

n = fix(m / 2);

zu = (0:n-1)'/(n-1); 

r = [ones(n, 1); (1 - zu)];
g = [zu; (1 - zu)];
b = [zu; ones(n, 1)];

vba = [r g b];


