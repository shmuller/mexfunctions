function out = sm_mdscd(varargin)
%out = sm_mdsdir([sock],d) - Change directory in tree
%
%   S. H. Muller, 2008/08/25

[sock,args] = getsock(varargin);

if isempty(args)
    args{1} = '\TOP';
end

if strcmp(args{1},'..')
    args{1} = '.-';
else
    args{1} = strrep(args{1},'\','\\');
end
cmd = sprintf('%s ',args{:});

stat = sm_mdsvalue(sock,sprintf('Tcl("set def %s")',cmd));
if ~bitand(uint32(stat),1)
	warning('Could not change directory');
end

out = sm_mdsdir(sock);

