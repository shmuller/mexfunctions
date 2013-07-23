try:
    import numpypy as np
except ImportError:
    import numpy as np

import specfun

N = 1000
x = np.arange(N) * (3. / (N - 1))
y = x.copy()

specfun.specfun('besei0', y)

print y[:20]

try:
    from matplotlib.pylab import plot, show
except ImportError:
    pass
else:
    plot(x, y)
    show()
