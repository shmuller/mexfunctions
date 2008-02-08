switch computer
    case 'PC_WIN'
        opt = {'OPTIMFLAGS=-O3'};
        mds = {'-I"C:\Program Files\MDSplus\DEVTOOLS\include"','-L"C:\Program Files\MDSplus\DEVTOOLS\lib"'};
    case 'GLNX86'
        opt = {'COPTIMFLAGS=-O3','LDOPTIMFLAGS=-O3'};
        mds = {'-I/usr/local/mdsplus/include','-L/usr/local/mdsplus/lib'};
end

mex('-v','mdsclient.c',opt{:},mds{:},'-lMdsIpShr');
