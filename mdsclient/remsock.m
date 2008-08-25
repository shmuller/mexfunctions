function remsock(sock)
%remsock(sock) - Remove socket from persistent list
%
%   S. H. Muller, 2008/08/15

socks = get(0,'UserData');
socks(socks == sock) = [];
set(0,'UserData',socks);