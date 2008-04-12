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
        mds = {'-I/usr/local/mdsplus/include','-L/usr/local/mdsplus/lib','-lMdsIpShr'};
end

mex('-v','mdsclient.c',def{:},opt{:},mds{:});
