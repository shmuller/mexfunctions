from _mdsclient import mdsconnect, mdsdisconnect, mdsvalue

def mdsopen(sock, tree, shot):
    return mdsvalue(sock, 'TreeOpen($,$)', tree, shot)

def mdsclose(sock):
    return mdsvalue(sock, 'TreeClose()')

