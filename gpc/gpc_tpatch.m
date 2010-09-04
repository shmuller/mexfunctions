function h = gpc_tpatch(xy,n,c,a)

%bc = get(gca,'Color');
bc = [1,1,1];
S0 = struct('xy',cell(1,length(n)),'c',cell(1,length(n)));
jj = 0;
for j = 1:length(n)
    S0(j).xy = xy(:,jj+(1:n(j)));
    S0(j).c = alphablend(bc,c(j,:),a(j));
    jj = jj+n(j);
end

S = struct('xy',{},'c',{});
for j = 1:length(S0)
    sj.xy = S0(j).xy;
    sj.c = S0(j).c;
    for i = 1:length(S)
        si.xy = S(i).xy;
        si.c = S(i).c;
        [sji.xy,nji,si.xy,ni,sj.xy,nj] = gpcmex(si.xy,sj.xy);
        sji.c = alphablend(si.c,sj.c,a(j));
        S(i) = si;
        S = [S,sji];
    end
    S = [S,sj];
end

for i = 1:length(S)
    s = S(i);
    patch(s.xy(1,:),s.xy(2,:),0,'FaceColor',s.c);
end

% [XY,N,C] = deal([]);
% jj = 0;
% for j = 1:length(n);
%     xyj = xy(:,jj+(1:n(j)));
%     jj = jj+n(j);
%     cj = alphablend(bc,c(j,:),a(j));
%     XY = [XY,xyj];
%     N = [N; n(j)];
%     C = [C; cj];
%     
%     ii = 0;
%     for i = 1:length(N)-1
%         xyi = XY(:,ii+(1:N(i)));
%         ii = ii+N(i);
%         cji = alphablend(C(i,:),cj,a(j));
%         [xyji,nji] = gpcmex(xyj,xyi);
%         XY = [XY,xyji];
%         N = [N; nji];
%         C = [C; cji(ones(length(nji),1),:)];
%     end
% end
%
% h = patchgroup(XY,N,C);


%--------------------------------------------------------------------------
function c = alphablend(cb,ct,at)
c = ct.*at+cb.*(1-at);

%--------------------------------------------------------------------------
function S = split(s,n)
jj = 0;
S = struct('xy',{},'c',{});
n = abs(n);
for j = 1:length(n)
    if n(j) > 0
        sj.xy = s.xy(:,jj+(1:n(j)));
        sj.c = s.c;
        S(end+1) = sj;
        jj = jj+n(j);
    else
        jj = jj-n(j);
    end
end


%--------------------------------------------------------------------------
function h = patchgroup(xy,n,c)
xy = xy.';
N = length(n);
h = zeros(1,N);
ii = 0;
for i = 1:N
    h(i) = patch(xy(ii+1:ii+n(i),1),xy(ii+1:ii+n(i),2),0,'FaceColor',c(i,:));
    ii = ii+n(i);
end

