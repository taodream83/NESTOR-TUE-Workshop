function y = b2d(x)
% Fast binary to decimal function, supports matrix inputs

y = sum(bsxfun(@times, x, 2.^(size(x,2)-1:-1:0)), 2);
end
