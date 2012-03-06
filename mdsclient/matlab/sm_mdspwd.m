function out = sm_mdspwd(varargin)
%out = sm_mdspwd([sock]) - Show current directory in tree
%
%   S. H. Muller, 2008/08/25

sock = getsock(varargin);

out = sm_mdsvalue(sock,'Tcl("show def",_out); _out');