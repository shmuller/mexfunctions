import numpy as np
import odepack

neq = 2
n = 4

def dy_dt(y, t, ydot):
    ydot[0] = -y[0]
    ydot[1] = -y[1]

res = np.ones((n, neq))

time = np.arange(n, dtype=np.float64)

args = np.zeros(neq), np.zeros(neq), np.zeros(neq)

out = odepack.odesolve(dy_dt, res, time, args)

print out[:, 0]
print out[:, 1]
print np.exp(-time)

print args

