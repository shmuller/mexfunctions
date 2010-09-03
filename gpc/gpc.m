function [x,y] = gpc(x1,y1,x2,y2)

if ~isa(x1,'double'); x1 = double(x1); end
if ~isa(y1,'double'); y1 = double(y1); end
if ~isa(x2,'double'); x2 = double(x2); end
if ~isa(y2,'double'); y2 = double(y2); end

[xy,n] = gpcmex([x1(:),y1(:)].',[x2(:),y2(:)].');
xy = xy.';

N = length(n);
x = cell(1,N);
y = cell(1,N);
ii = 0;
for i = 1:N
    x{i} = xy(ii+1:ii+n(i),1);
    y{i} = xy(ii+1:ii+n(i),2);
    ii = ii+n(i);
end
    