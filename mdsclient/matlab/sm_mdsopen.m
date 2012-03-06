function out = sm_mdsopen(varargin)
%out = sm_mdsopen(sock,tree,shot)
%   Wrapper for mdsclientmex file
%
%   S. H. Muller, 2008/02/07

[sock,args] = getsock(varargin);

expt = args{1};
if length(args) > 1
    shot = int32(args{2});
else
    shot = int32(0);
end

out = sm_mdsvalue(sock,'TreeOpen($,$)',expt,shot);
