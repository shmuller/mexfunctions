function sm_mdsdisconnect(varargin)
%sm_mdsdisconnect([sock])
%   Wrapper for mdsclientmex file
%
%   S. H. Muller, 2008/02/07

sock = getsock(varargin);

mdsclientmex('mdsdisconnect',sock);

remsock(sock);
