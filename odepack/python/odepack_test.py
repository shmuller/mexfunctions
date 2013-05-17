import numpy as np
import odepack

neq = 2
n = 7

def dy_dt(y, t, ydot):
    ydot[0] = -y[0]
    ydot[1] = -y[1]

def term(y, t, ydot):
    return t[0] > 3.

#term = None

res = np.ones((n, neq))

time = np.arange(n, dtype=np.float64)

args = np.zeros(neq), np.zeros(1), np.zeros(neq)

points_done = odepack.odesolve(dy_dt, res, time, args, term)

t = time[:points_done]
out = res[:points_done]

print out[:, 0]
print out[:, 1]
print np.exp(-t)

print args

