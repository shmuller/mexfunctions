function addsock(sock)
%addsock(sock) - Add socket to persistent list
%
%   S. H. Muller, 2008/08/15

set(0,'UserData',[sock,get(0,'UserData')]);