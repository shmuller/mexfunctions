import numpy as np
from odepack import odesolve

neq = 2
n = 6
ng = 1

def f(y, t, ydot):
    ydot[0] = -y[0]
    ydot[1] = -y[1]

y = np.ones((n, neq))
t = np.arange(n, dtype=np.float64)

work = np.zeros(neq), np.zeros(1), np.zeros(neq), np.zeros(ng)

points_done = odesolve(f, y, t, work)
t1, y1 = t[:points_done], y[:points_done]

print y1[:,0]
print np.exp(-t1)


def g(y, t, gout):
    gout[0] = y[0] - 0.2

points_done = odesolve(f, y, t, work, g)
t2, y2 = t[:points_done], y[:points_done]

print y2[:,0]
print np.exp(-t2)


ibbox = np.array([0], 'i')
bbox = np.array([0.2])

points_done = odesolve(f, y, t, work, None, ibbox, bbox)
t3, y3 = t[:points_done], y[:points_done]

print y3[:,0]
print np.exp(-t3)


