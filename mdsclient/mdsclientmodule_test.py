import mdsclient

sock = mdsclient.mdsconnect('localhost:8010')
print sock


res = mdsclient.mdsvalue(sock,'[[1,2,3],[4,5,6]]')
print res

res = mdsclient.mdsvalue(sock,'[1D0,2D0,3D0]')
print res

res = mdsclient.mdsvalue(sock,'CMPLX(1.,1.)')
print res

res = mdsclient.mdsvalue(sock,'"Hallo"')
print res

res = mdsclient.mdsvalue(sock,'[]')
print res

res = mdsclient.mdsvalue(sock,'""')
print res

mdsclient.mdsdisconnect(sock)


