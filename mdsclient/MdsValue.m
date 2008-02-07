function out = MdsValue(sock,expression)
%out = MdsValue(sock,expression)
%   Evaluates TDI expression
%
%   S. H. Muller, 2008/02/07

out = mdsclient('mdsvalue',sock,expression);