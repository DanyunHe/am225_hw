import numpy as np 
import matplotlib.pyplot as plt 

data=np.loadtxt("./out2")

fig,ax=plt.subplots()
ax.plot(data.T[1],data.T[2],'o-',label="l2 norm")
ax.plot(data.T[1],pow(data.T[1],1.7),label="power law")

ax.legend()
ax.set_xlabel("h")
ax.set_ylabel("error")
ax.set_yscale("log")
ax.set_xscale("log")
plt.show()





