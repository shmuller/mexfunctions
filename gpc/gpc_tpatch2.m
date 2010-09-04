function h = gpc_tpatch2(xy,n,c,a)

%bc = get(gca,'Color');
bc = [1,1,1];

[XY,N,C] = deal([]);
jj = 0;
for j = 1:length(n);
    xyj = xy(:,jj+(1:n(j)));
    jj = jj+n(j);
    cj = alphablend(bc,c(j,:),a(j));
    XY = [XY,xyj];
    N = [N; n(j)];
    C = [C; cj];
    
    ii = 0;
    for i = 1:length(N)-1
        xyi = XY(:,ii+(1:N(i)));
        ii = ii+N(i);
        cji = alphablend(C(i,:),cj,a(j));
        [xyji,nji] = gpcmex(xyj,xyi);
        XY = [XY,xyji];
        N = [N; nji];
        C = [C; cji(ones(length(nji),1),:)];
    end
end

h = patchgroup(XY,N,C);


%--------------------------------------------------------------------------
function c = alphablend(cb,ct,at)
c = ct.*at+cb.*(1-at);

%--------------------------------------------------------------------------
function h = patchgroup(xy,n,c)
xy = xy.';
N = length(n);
h = zeros(1,N);
ii = 0;
for i = 1:N
    z = i(ones(n(i),1));
    h(i) = patch(xy(ii+1:ii+n(i),1),xy(ii+1:ii+n(i),2),z,0,'FaceColor',c(i,:),'EdgeColor','none');
    ii = ii+n(i);
end

