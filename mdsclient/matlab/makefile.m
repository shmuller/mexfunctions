if sm_verlessthan(version,'7.4')
    def = {'-DCOMPAT_V70'};
else
    def = {};
end

switch computer
    case 'PCWIN'
        opt = {'-O'};
        inc = {'-IC:\PROGRA~1\MDSplus\DEVTOOLS\include'};
        lib = {'C:\PROGRA~1\MDSplus\DEVTOOLS\lib\MdsIpShr.lib','C:\MinGW\lib\libws2_32.a'};
    case 'GLNX86'
        opt = {'COPTIMFLAGS=-O3','LDOPTIMFLAGS=-O3'};
        inc = {'-I/usr/local/mdsplus/include'};
        lib = {'/usr/local/mdsplus/lib32/libMdsIpShr.a'};
    case 'GLNXA64'
        opt = {'COPTIMFLAGS=-O3','LDOPTIMFLAGS=-O3'};
        inc = {'-I/usr/local/mdsplus/include'};
        lib = {'/usr/local/mdsplus/lib64/libMdsIpShr.a'};
    case {'MAC','MACI'}
        opt = {'COPTIMFLAGS=-O3','LDOPTIMFLAGS=-O3'};
        inc = {'-I/usr/local/mdsplus/include'};
        lib = {'/usr/local/mdsplus/lib/libMdsIpShr.a'};
end

mex('-v','mdsclientmex.c',def{:},opt{:},inc{:},lib{:});
