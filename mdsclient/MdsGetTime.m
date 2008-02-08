function t = MdsGetTime(sock,node)
%MdsGetTime(sock,node)
%   Get signal time base
%
%   S. H. Muller, 2008/02/07

expr = ['_t=dim_of(',node,'); [slope_of(axis_of(_t)),_t[0]]'];
t = 1000*mdsclient('mdsvalue',sock,expr);

expr = '_i=dim_of(data(_t)); [_i[0],_i[size(_i)-1]]';
i = mdsclient('mdsvalue',sock,expr);

t = cast((i(1):i(2))',class(t))*t(1)+t(2);
