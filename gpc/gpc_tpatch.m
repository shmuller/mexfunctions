function h = gpc_tpatch(xy,n,c,a)

[XY,N,M] = gpc_cut_mex(xy,int32(n));

%bc = get(gca,'Color');
bc = [1,1,1];

m = length(n);
L = pow2(m)-1;
C = zeros(L,3);
p = 1;
for j = 1:length(n)
    C(p,:) = alphablend(bc,c(j,:),a(j));
    pi = p+1;
    for q = 1:p-1
        C(pi,:) = alphablend(C(q,:),C(p,:),a(j));
        pi = pi+1;
    end
    p = pi;
end

tpatch(XY,N,M,C);


%--------------------------------------------------------------------------
function c = alphablend(cb,ct,at)
c = ct.*at+cb.*(1-at);


%--------------------------------------------------------------------------
function tpatch(XY,N,M,C)
jj = 0;
ii = 0;
J = length(M);
for j = 1:J
    nj = N(ii+(1:M(j)));
    ii = ii+M(j);
    cj = C(j,:);
    
    for i = 1:M(j)
        nji = abs(nj(i));
        xyji = XY(:,jj+(1:nji));
        jj = jj+nji;
        
        patch(xyji(1,:),xyji(2,:),j(ones(1,nji)),0,'FaceColor',cj,'EdgeColor','none');
    end
end
