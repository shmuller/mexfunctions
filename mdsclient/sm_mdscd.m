function out = sm_mdscd(sock,d)
%out = sm_mdsdir(sock,d) - Change directory in tree
%
%   S. H. Muller, 2008/08/25

switch d
	case '..'
		d = '.-';
	case '/'
		d = '\TOP';
end
d = strrep(d,'\','\\');

stat = sm_mdsvalue(sock,sprintf('Tcl("set def %s")',d));
if ~bitand(uint32(stat),1)
	warning('Could not change directory');
end

out = sm_mdsdir(sock);

