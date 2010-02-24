function y = mediansmooth(x,w)
%y = mediansmooth(x,w)
%   Sliding median filter of width 2*w+1.
%
%   S. H. Muller, 2010/02/23

[N,J] = size(x);
y = zeros(size(x));

for j = 1:J
    ind = mediansmoothmex(x(:,j),int32(w));
    y(:,j) = x(ind+1,j);
end