import numpy as np
import odepack

neq = 2
n = 7
ng = 1

def f(y, t, ydot):
    ydot[0] = -y[0]
    ydot[1] = -y[1]

def g(y, t, gout):
    gout[0] = y[0] - 0.2

y = np.ones((n, neq))

t = np.arange(n, dtype=np.float64)

work = np.zeros(neq), np.zeros(1), np.zeros(neq), np.zeros(ng)

points_done = odepack.odesolve(f, y, t, work, g)

t = t[:points_done]
y = y[:points_done]

print y[:, 0]
print y[:, 1]
print np.exp(-t)

print work

