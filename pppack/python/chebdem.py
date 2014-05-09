import numpy as np
import splinetoolbox

diff, cat = np.diff, np.concatenate

k = 4
brks = np.array([0, 1, 1.1, 3, 5, 5.5, 7, 7.1, 7.2, 8])

t = splinetoolbox.augknt(brks, k)
tau = splinetoolbox.aveknt(t, k)

n = len(t) - k
b = (-1.)**np.arange(n-1, -1, -1)

sp = splinetoolbox.SplinePGS(tau, b, t=t)

#sp.plot()

dsp = sp.deriv()

tstar = dsp.tave()
coefs = dsp.c.ravel()

guess = tstar[:-1] - coefs[:-1] * diff(tstar) / diff(coefs)

sites = np.zeros((k, len(guess)))
sites[0] = guess
sites[1] = tau[1:-1]

values = np.zeros((k, len(guess)))
values[0] = dsp(guess).ravel()
values[1] = dsp(tau[1:-1]).ravel()

for j in xrange(1, 3):
   dvalues = values[j] - values[j-1]
   dvalues[dvalues == 0] = 1.

   sites[j+1] = sites[j] - values[j] * (sites[j] - sites[j-1]) / dvalues
   values[j+1] = dsp(sites[j+1]).ravel()

taunew = np.r_[tau[0], sites[3], tau[n-1]]
spnew = splinetoolbox.SplinePGS(taunew, b, t=t)

dsp.plot()
dsp.plot_ctrlpoly()
plot(guess, np.zeros(guess.size), 'ro', taunew, np.zeros(taunew.size), 'bx')

figure()
spnew.plot()

show()
