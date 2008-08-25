function out = sm_mdsdir(varargin)
%out = sm_mdsdir([sock],d) - List directory in tree
%
%   S. H. Muller, 2008/08/25

[sock,args] = getsock(varargin);

if isempty(args)
	args{1} = '';
else
	args{1} = strrep(args{1},'\','\\');
end
cmd = sprintf('%s ',args{:});

out = sm_mdsvalue(sock,sprintf('Tcl("dir %s",_out); _out',cmd));