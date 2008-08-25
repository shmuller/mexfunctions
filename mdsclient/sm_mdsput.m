function out = sm_mdsput(varargin)
%out = sm_mdsput([sock],node,val) - Put value in node
%
%   S. H. Muller, 2008/08/25

[sock,args] = getsock(varargin);

out = sm_mdsvalue(sock,'TreePutRecord($,$)',args{:});