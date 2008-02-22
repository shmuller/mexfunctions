function MdsClose(sock)
%MdsClose(sock)
%   Closes opened tree and disconnects from server
%
%   S. H. Muller, 2008/02/07

expt = mdsclient('mdsvalue',sock,'$expt');
errstr = '%TREE-W-TreeNOT_OPEN, Tree not currently open';

if ~strcmp(expt,errstr)
    mdsclient('mdsclose',sock);
end

if isunix
    mdsclient('mdsdisconnect',sock);
end