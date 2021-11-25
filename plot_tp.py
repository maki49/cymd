import numpy as np
import matplotlib.pyplot as plt

temp=np.loadtxt("temperature.txt")
t=np.linspace(0, 5, 251)
dash=np.ones(251)*200
plt.plot(t,temp);
plt.plot(t,dash, 'r-.');