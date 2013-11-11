import _mdsclient

MdsConnectError = _mdsclient.MdsConnectError

class TdiError(Exception):
    def __init__(self, err, mdsfmt, args):
        self.err, self.mdsfmt, self.args = err, mdsfmt, args

    def __str__(self):
        return self.err + '\nmdsfmt:\n' + self.mdsfmt + '\nargs:\n' + pformat(self.args)


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

