function sm_mdsdisconnect(sock)
%sm_mdsdisconnect(sock)
%   Wrapper for mdsclient library
%
%   S. H. Muller, 2008/02/07

mdsclient('mdsdisconnect',sock);
