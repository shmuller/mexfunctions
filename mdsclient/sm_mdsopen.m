function out = sm_mdsopen(varargin)
%out = sm_mdsopen(sock,tree,shot)
%   Wrapper for mdsclient mex file
%
%   S. H. Muller, 2008/02/07

[sock,args] = getsock(varargin);

out = sm_mdsvalue(sock,'TreeOpen($,$)',args{:});