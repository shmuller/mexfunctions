function h = gpc_tpatch(xyca)

bc = get(gca,'Color');

[xy1,N1,c1,a1] = deal(xyca{1}{:});
[xy2,N2,c2,a2] = deal(xyca{2}{:});

xy1 = xy1(:,1:N1(1));
xy2 = xy2(:,1:N2(1));

c1 = alphablend(bc,c1,a1);
c2 = alphablend(bc,c2,a2);
ci = alphablend(c1,c2,a2);

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

