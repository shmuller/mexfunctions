function [x,t] = MdsGetSig(sock,node)
%MdsGetSig(sock,node)
%   Get signal
%
%   S. H. Muller, 2008/02/07

% expr = ['getsigcal(',node,')'];
% f = mdsclient('mdsvalue',sock,expr);
% 
% expr = ['raw_of(',node,')'];
% x = mdsclient('mdsvalue',sock,expr);
% 
% x = f(1)*cast(x,class(f))+f(2);

x = mdsclient('mdsvalue',sock,node);

if nargout > 1
    t = MdsGetTime(sock,node);
    t = t(1:length(x));
end
