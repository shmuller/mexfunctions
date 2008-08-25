function sm_mdsdisconnect(varargin)
%sm_mdsdisconnect([sock])
%   Wrapper for mdsclient library
%
%   S. H. Muller, 2008/02/07

sock = getsock(varargin);

mdsclient('mdsdisconnect',sock);

remsock(sock);