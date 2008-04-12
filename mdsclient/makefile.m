if sm_verlessthan(version,'7.4')
    def = {'-DCOMPAT_V70'};
else
    def = {};
end

switch computer
    case 'PCWIN'
        opt = {'-O'};
        mds = {'-IC:\PROGRA~1\MDSplus\DEVTOOLS\include','C:\PROGRA~1\MDSplus\DEVTOOLS\lib\MdsIpShr.lib'};
    case 'GLNX86'
        opt = {'COPTIMFLAGS=-O3','LDOPTIMFLAGS=-O3'};
        mds = {'-I/usr/local/mdsplus/include','/usr/local/mdsplus/lib32/libMdsIpShr.a'};
    case 'GLNXA64'
        opt = {'COPTIMFLAGS=-O3','LDOPTIMFLAGS=-O3'};
        mds = {'-I/usr/local/mdsplus/include','/usr/local/mdsplus/lib64/libMdsIpShr.a'};
end

mex('-v','mdsclient.c',def{:},opt{:},mds{:});
