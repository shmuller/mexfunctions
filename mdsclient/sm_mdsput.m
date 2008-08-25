function out = sm_mdsput(sock,node,val)
%out = sm_mdsdir(sock,node,val) - Put value in node
%
%   S. H. Muller, 2008/08/25

out = sm_mdsvalue(sock,'TreePutRecord($,$)',node,val);