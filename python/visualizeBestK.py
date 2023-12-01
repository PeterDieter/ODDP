import numpy as np
import matplotlib.pyplot as plt

matrix = np.loadtxt("tuneKResults.txt", delimiter=' ')

plt.plot(matrix[:,2])
plt.yscale("log")
plt.show()