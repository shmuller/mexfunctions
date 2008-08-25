function out = sm_mdsclose(varargin)
%out = sm_mdsclose(sock)
%   Wrapper for mdsclient mex file
%
%   S. H. Muller, 2008/02/07

sock = getsock(varargin);

out = sm_mdsvalue(sock,'TreeClose()');