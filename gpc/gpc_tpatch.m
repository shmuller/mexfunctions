function h = gpc_tpatch(x,y,c,a)

bc = get(gca,'Color');

c1 = alphablend(bc,c(1,:),a(1));
c2 = alphablend(bc,c(2,:),a(2));
ci = alphablend(c1,c2,a(2));

xy1 = [x{1}(:),y{1}(:)].';
xy2 = [x{2}(:),y{2}(:)].';

tic
[xyi,ni,xy1,n1,xy2,n2] = gpcmex(xy1,xy2);
toc

xy = [xy1.'; xy2.'; xyi.'];
n = [n1; n2; ni];
c = [c1(ones(length(n1),1),:); c2(ones(length(n2),1),:); ci(ones(length(ni),1),:)];

h = patchgroup(xy,n,c);


%--------------------------------------------------------------------------
function c = alphablend(cb,ct,at)
c = ct.*at+cb.*(1-at);

%--------------------------------------------------------------------------
function h = patchgroup(xy,n,c)
N = length(n);
h = zeros(1,N);
ii = 0;
for i = 1:N
    h(i) = patch(xy(ii+1:ii+n(i),1),xy(ii+1:ii+n(i),2),0,'FaceColor',c(i,:));
    ii = ii+n(i);
end

