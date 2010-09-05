m = 20;
phi = (0:m-1)*(2*pi/m);
cph = cos(phi); sph = sin(phi);

n1 = 4;
%x1 = 2*rand(n1,1);
%y1 = 2*rand(n1,1);
x1 = [1.2,1,2,2].';
y1 = [1,1,1,1].';
r1 = (1:n1).'./n1;

[xy,n] = deal([]);
for j = 1:n1
    xy = [xy, [x1(j)+r1(j)*cph; y1(j)+r1(j)*sph]];
    n = [n; m];
end
c = [0,0,1; 0,1,0; 1,0,0; 0,0,0];
%c = jet(n1);
a = 0.9*ones(n1,1);

[XY,N,M] = gpc_cut_mex(xy,int32(n));

%gpc_tpatch(xy,n,c,a);

%xlim([0,3]);
%ylim([0,3]);
%axis square;

jj = 0;
ii = 0;
J = length(M);
c = 0;
for j = 1:J
    nc{j} = N(ii+(1:M(j)));
    ii = ii+M(j);
    
    if M(j) ~= 0
        c = c+1;
    end
    
    for i = 1:M(j)
        nji = abs(nc{j}(i));
        xyc{i,j} = XY(:,jj+(1:nji));
        jj = jj+nji;
        
        if nc{j}(i) > 0    
            sty = {};
        else
            sty = {'LineStyle','--'};
        end
        line(xyc{i,j}(1,:),xyc{i,j}(2,:),c(ones(1,nji)),sty{:});
    end
end



