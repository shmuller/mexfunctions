try:
    import numpypy as np
except ImportError:
    import numpy as np

from mediansmooth import *

N = 1000000
#x = np.linspace(0, 2*np.pi, N)
x = np.arange(N) * (2*np.pi / (N - 1))
y = np.sin(x)

#y += 0.5*(2*np.random.random(x.size)-1)

ind2 = mediansmooth_old(y, 100)

print ind2[:20]

try:
    import matplotlib.pyplot as plt
except ImportError:
    pass
else:
    plt.plot(x, y, x, y[ind2])
    plt.show()

