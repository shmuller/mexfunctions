function sm_mdsopen(sock,tree,shot)
%sm_mdsopen(sock,tree,shot)
%   Wrapper for mdsclient mex file
%
%   S. H. Muller, 2008/02/07

mdsclient('mdsopen',sock,tree,int32(shot));