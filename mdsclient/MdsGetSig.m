function [x,t] = MdsGetSig(sock,node,mode)
%MdsGetSig(sock,node)
%   Get signal
%
%   S. H. Muller, 2008/02/07

if nargin < 3
    mode = 'full';
end

switch mode
    case 'full'
        x = mdsclient('mdsvalue',sock,node);
        
    case 'raw'
        expr = ['raw_of(',node,')'];
        x = mdsclient('mdsvalue',sock,expr);
        
    case 'netw'
        expr = ['getsigcal(',node,')'];
        f = mdsclient('mdsvalue',sock,expr);

        expr = ['raw_of(',node,')'];
        x = mdsclient('mdsvalue',sock,expr);

        x = f(1)*cast(x,class(f))+f(2);
end

if nargout > 1
    t = MdsGetTime(sock,node);
    t = t(1:length(x));
end
