import numpy as np
import matplotlib.pyplot as plt

test = np.loadtxt('./out1b')
plt.plot(test.T[0],test.T[1],'.-')
plt.plot(test.T[0],test.T[2],'.-')
x = np.linspace(0,8,1000)
plt.show()