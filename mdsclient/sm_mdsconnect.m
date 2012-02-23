function sock = sm_mdsconnect(serv)
%sock = sm_mdsconnect(serv)
%   Wrapper for mdsclientmex file
%
%   S. H. Muller, 2008/02/07

sock = mdsclientmex('mdsconnect',serv);

addsock(sock);
