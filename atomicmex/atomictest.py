import numpy as np
import atomic
from matplotlib.pylab import plot, show

w = np.linspace(0,1e8,1000)
sigv_i = w.copy()
sigv_CX = w.copy()

atomic.atomic(sigv_i , 'D', 'ion', 'BEB')
atomic.atomic(sigv_CX, 'D', 'CX')

plot(w, sigv_i, w, sigv_CX)
show()
