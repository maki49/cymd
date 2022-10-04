import numpy as np
import matplotlib.pyplot as plt
d=np.loadtxt("diff.txt")
plt.figure()
plt.plot(d[:,0], d[:,1])
plt.savefig("instable.png")
