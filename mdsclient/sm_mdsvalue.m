function out = sm_mdsvalue(varargin)
%out = sm_mdsvalue(sock,expression,var1,var2,...)
%   Wrapper for mdsclientmex file
%
%   S. H. Muller, 2008/02/07

[sock,args] = getsock(varargin);

out = mdsclientmex('mdsvalue',sock,args{:});

if ischar(out) && any(strfind(out,'TDI Error'))
    error(out);
end
