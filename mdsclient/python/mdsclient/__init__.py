import _mdsclient

def mdsconnect(host):
    return _mdsclient.mdsconnect(host)

def mdsdisconnect(sock):
    return _mdsclient.mdsdisconnect(sock)

def mdsvalue(sock, expr, *args):
    return _mdsclient.mdsvalue(sock, expr, *args)

def mdsopen(sock, tree, shot):
    return mdsvalue(sock, 'TreeOpen($,$)', tree, shot)

def mdsclose(sock):
    return mdsvalue(sock, 'TreeClose()')

