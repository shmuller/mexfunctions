function y = mediansmooth(x,w,corr_bdry)
%y = mediansmooth(x,w,corr_bdry)
%   Sliding median filter of width 2*w+1.
%
%   S. H. Muller, 2010/02/23

if nargin < 3
    corr_bdry = 1;
end

[N,J] = size(x);
y = zeros(size(x),class(x));

if corr_bdry
    for j = 1:J
        [x0,x1] = sm_extrap(x(:,j),w,'polyfit',1);
        xj = [x0; x(:,j); x1];

        ind = mediansmoothmex(xj,int32(w));
        y(:,j) = xj(ind(w+1:w+N)+1);
    end
else
    for j = 1:J
        ind = mediansmoothmex(x(:,j),int32(w));
        y(:,j) = x(ind+1,j);
    end
end
