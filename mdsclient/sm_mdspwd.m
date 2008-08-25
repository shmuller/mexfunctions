function out = sm_mdspwd(sock)
%out = sm_mdspwd(sock) - Show current directory in tree
%
%   S. H. Muller, 2008/08/25

out = sm_mdsvalue(sock,'Tcl("show def",_out); _out');