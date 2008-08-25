function out = sm_mdsvalue(varargin)
%out = sm_mdsvalue(sock,expression,var1,var2,...)
%   Wrapper for mdsclient mex file
%
%   S. H. Muller, 2008/02/07

[sock,args] = getsock(varargin);

out = mdsclient('mdsvalue',sock,args{:});

if ischar(out) && any(strfind(out,'TDI Error'))
    error(out);
end