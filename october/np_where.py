import numpy as np
import matplotlib.pyplot as plt

plt.figure()
a = np.arange(1, 1000)
plt.plot(a, label='before')
a[np.logical_and(200 < a, a < 500)] = -5
plt.plot(a, label='after')
plt.show()
