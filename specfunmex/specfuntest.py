import numpy as np
import specfun
from matplotlib.pylab import plot, show

x = np.linspace(0.,3.,1000)
y = x.copy()

specfun.specfun('besei0', y)

plot(x, y)
show()
