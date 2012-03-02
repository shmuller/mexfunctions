import mdsclient

sock = mdsclient.mdsconnect('localhost:8010')
print sock

res = mdsclient.mdsvalue(sock,'[1,2,3]')
print res

mdsclient.mdsdisconnect(sock)


