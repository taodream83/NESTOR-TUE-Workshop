function y = d2b(x,m)
% Fast decimal to binary function, supports vector inputs
% Set number of output bits with m

c = floor(log2(x)) + 1;  % Number of divisions necessary ( rounding up the log2(x) )
y = zeros(numel(x),m);   % Initialize output matrix
for i = 1:max(c)
    r = floor(x / 2);
    y(:,m+1-i) = x - 2.*r;
    x = r;
end