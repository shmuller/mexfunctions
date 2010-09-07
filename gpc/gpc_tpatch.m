function h = gpc_tpatch(xy,n,c,a,varargin)
%h = gpc_tpatch(xy,n,c,a,...)
%   Transparent patches.
%
%   S. H. Muller, 2010/09/04

if ~isa(xy,'double')
    xy = double(xy);
end

[XY,N,M] = gpc_cut_mex(xy,int32(n));

bc = get(gca,'Color');

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

h = tpatch(XY,N,M,C,varargin{:});


%--------------------------------------------------------------------------
function c = alphablend(cb,ct,at)
c = ct.*at+cb.*(1-at);


%--------------------------------------------------------------------------
function h = tpatch(XY,N,M,C,varargin)
jj = 0;
ii = 0;
k = 0;
J = length(M);
h = zeros(1,length(N));
% loop over polygon groups
for j = 1:J
    nj = N(ii+(1:M(j)));
    ii = ii+M(j);
    cj = C(j,:);
    
    % loop over polygon parts
    for i = 1:M(j)
        nji = abs(nj(i));
        xyji = XY(:,jj+(1:nji));
        jj = jj+nji;
        
        % don't plot holes
        if nj(i) > 0
            k = k+1;
            h(k) = patch(xyji(1,:),xyji(2,:),0,'FaceColor',cj,varargin{:});
        end
    end
end
h(k+1:end) = [];
