function out = sm_mdsvalue(sock,expression)
%out = sm_mdsvalue(sock,expression)
%   Wrapper for mdsclient mex file
%
%   S. H. Muller, 2008/02/07

out = mdsclient('mdsvalue',sock,expression);