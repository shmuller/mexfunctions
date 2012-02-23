function sock = MdsOpen(serv,tree,shot)
%sock = MdsOpen(serv,tree,shot)
%   Connects to server and opens a tree
%
%   S. H. Muller, 2008/02/07

sock = mdsclientmex('mdsconnect',serv);

switch nargin
    case 2
        mdsclientmex('mdsopen',sock,tree,int32(0));
    case 3
        mdsclientmex('mdsopen',sock,tree,int32(shot));
end
