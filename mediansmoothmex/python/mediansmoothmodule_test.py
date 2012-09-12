import numpy as np

import matplotlib.pyplot as plt

from mediansmooth import *

x = np.linspace(0,2*np.pi,1000000);
y = np.sin(x) + 0.5*(2*np.random.random(x.size)-1)

ind2 = mediansmooth2(y, 1000000)

plt.plot(x, y, x, y[ind2])

#plt.show()

ind2 = mediansmooth2(y, 1000000)

