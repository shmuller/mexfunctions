switch computer
    case 'PCWIN'
        opt = {'OPTIMFLAGS=-O3'};
        mds = {'-I"C:\PROGRA~1\MDSplus\DEVTOOLS\include"','-L"C:\PROGRA~1\MDSplus\DEVTOOLS\lib"'};
    case 'GLNX86'
        opt = {'COPTIMFLAGS=-O3','LDOPTIMFLAGS=-O3'};
        mds = {'-I/usr/local/mdsplus/include','-L/usr/local/mdsplus/lib'};
end

mex('-v','mdsclient.c',opt{:},mds{:},'-lMdsIpShr');
