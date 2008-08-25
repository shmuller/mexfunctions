function [sock,args] = getsock(args)
%[sock,args] = getsock(args) - Obtain and separate socket
%
%   S. H. Muller, 2008/08/25

socks = get(0,'UserData');
if isempty(socks)
	error('No connection currently open');
end

if isempty(args) || ischar(args{1})
    sock = socks(1);
else
    sock = args{1};
    args(1) = [];
end
k = find(socks == sock);
set(0,'UserData',socks([k,1:k-1,k+1:end]));