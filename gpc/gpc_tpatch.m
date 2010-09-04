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

E = struct('xy',{},'c',{});
S = E;
for j = 1:length(S0)
    sj = S0(j);
    ih = 0;
    for i = 1:length(S)
        si = S(i);
        [sji.xy,nji,xyi,ni,xyj,nj] = gpcmex(si.xy,sj.xy);
        sji.c = alphablend(si.c,sj.c,a(j));
            
        nji, ni, nj
        if length(ni) == 1
            S(i).xy = xyi;
        end
        if length(nj) == 1
            sj.xy = xyj;
        end
        if length(ni) == 0
            %ih = i;
        end
        S = [S,split(sji,nji)];
    end
    if ih == 0
        S = [S,sj];
    else
        S(ih) = sj;
    end
end

for i = 1:length(S)
    s = S(i);
    z = i(ones(1,size(s.xy,2)));
    patch(s.xy(1,:),s.xy(2,:),z,0,'FaceColor',s.c,'EdgeColor','none');
end


%--------------------------------------------------------------------------
function c = alphablend(cb,ct,at)
c = ct.*at+cb.*(1-at);

%--------------------------------------------------------------------------
function S = split(s,n)
jj = 0;
S = struct('xy',{},'c',{});
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

