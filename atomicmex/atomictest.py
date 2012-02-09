import numpy as np
import atomic
from matplotlib.pylab import plot, show

w = np.linspace(0,1e8,1000)
sigv = w.copy()

atomic.atomic(sigv)

plot(w,sigv)
show()
