import numpy as np 
import matplotlib.pyplot as plt 

time=np.loadtxt("./out1")

fig,ax=plt.subplots()
x=np.zeros(7)
n=16
for i in range(7):
	x[i]=n
	n*=2

# fit 
log_x=np.log(x)
log_y=np.log(time.T[3])
A=np.vstack([log_x,np.ones(len(log_x))]).T
m,c=np.linalg.lstsq(A,log_y,rcond=None)[0]
print(m,c)


ax.plot(x,time.T[1],'o-',label="strassen")
ax.plot(x,time.T[2],'o-',label="standard")
ax.plot(x,time.T[3],'o-',label="blas")
ax.plot(x,3e-7*x**2.806,'--',label="power law")
ax.plot(x,5e-11*x**3.644,'--',label="O($n^3$)")
ax.plot(x,5e-10*x**2.451,'--',label="O($n^3$)")
ax.legend()
ax.set_xlabel("n")
ax.set_ylabel("time(s)")
ax.set_yscale("log")
ax.set_xscale("log")
plt.show()
