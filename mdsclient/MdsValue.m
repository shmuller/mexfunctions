function out = MdsValue(varargin)
%out = MdsValue(sock,expression,var1,var2,...)
%   Evaluates TDI expression
%
%   S. H. Muller, 2008/02/07

out = mdsclientmex('mdsvalue',varargin{:});
