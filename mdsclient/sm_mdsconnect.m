function sock = sm_mdsconnect(serv)
%sock = sm_mdsconnect(serv)
%   Wrapper for mdsclient mex file
%
%   S. H. Muller, 2008/02/07

sock = mdsclient('mdsconnect',serv);