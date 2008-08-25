function out = sm_mdsdir(sock,d)
%out = sm_mdsdir(sock,d) - List directory in tree
%
%   S. H. Muller, 2008/08/25

if nargin < 2
	d = '';
end
d = strrep(d,'\','\\');

out = sm_mdsvalue(sock,sprintf('Tcl("dir %s",_out); _out',d));