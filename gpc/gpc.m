function [xi,yi,x1,y1,x2,y2] = gpc(x1,y1,x2,y2)

if ~isa(x1,'double'); x1 = double(x1); end
if ~isa(y1,'double'); y1 = double(y1); end
if ~isa(x2,'double'); x2 = double(x2); end
if ~isa(y2,'double'); y2 = double(y2); end

if nargout < 3
    [xyi,n] = gpcmex([x1(:),y1(:)].',[x2(:),y2(:)].');
    [xi,yi] = parse_xy(xyi.',n);
else
    [xyi,ni,xy1,n1,xy2,n2] = gpcmex([x1(:),y1(:)].',[x2(:),y2(:)].');
    [xi,yi] = parse_xy(xyi.',ni);
    [x1,y1] = parse_xy(xy1.',n1);
    [x2,y2] = parse_xy(xy2.',n2);
end


%--------------------------------------------------------------------------
function [x,y] = parse_xy(xy,n)
N = length(n);
x = cell(1,N);
y = cell(1,N);
ii = 0;
for i = 1:N
    x{i} = xy(ii+1:ii+n(i),1);
    y{i} = xy(ii+1:ii+n(i),2);
    ii = ii+n(i);
end
    