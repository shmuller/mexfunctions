function sm_mdsclose(sock)
%sm_mdsclose(sock)
%   Wrapper for mdsclient mex file
%
%   S. H. Muller, 2008/02/07

mdsclient('mdsclose',sock);