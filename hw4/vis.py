import numpy as np 
import matplotlib.pyplot as plt 

data=np.loadtxt("./out2")

fig,ax=plt.subplots()
ax.plot(data.T[1],data.T[2],'o-',label="l2 norm")
ax.plot(data.T[1],pow(data.T[1],2.),label="$h^2$")

ax.legend()
ax.set_xlabel("h")
ax.set_ylabel("error")
ax.set_yscale("log")
ax.set_xscale("log")
plt.savefig("p2")
plt.close()


def basis_func(x,idx):
	if(idx==0):
		return 2*x**3-3*x**2+1
	if(idx==1):
		return x**3-2*x**2+x
	if(idx==2):
		return -2*x**3+3*x**2
	if(idx==3):
		return x**3-x**2

N=3
h=1/N

z0 = np.linspace(1,1+1/3,1000)
z1 = np.linspace(1+1/3,1+2/3,1000)
z2 = np.linspace(1+2/3,1+3/3,1000)


plt.plot(np.concatenate((z0,np.array([2]))),np.concatenate((basis_func((z0-1)/h,0),np.array([0]))),label="i=0")
plt.plot(np.concatenate((z0,np.array([2]))),0.6+np.concatenate((basis_func((z0-1)/h,1),np.array([0]))),label="i=1")
plt.plot(np.concatenate((z0,z1,np.array([2]))),1.2+np.concatenate((basis_func((z0-1)/h,2),basis_func((z1-1)/h-1,0),np.array([0]))),label="i=2")
plt.plot(np.concatenate((z0,z1,np.array([2]))),1.8+np.concatenate((basis_func((z0-1)/h,3),basis_func((z1-1)/h-1,1),np.array([0]))),label="i=3")
plt.plot(np.concatenate((np.array([1]),z1,z2)),2.4+np.concatenate((np.array([0]),basis_func((z1-1)/h-1,2),basis_func((z2-1)/h-2,0))),label="i=4")
plt.plot(np.concatenate((np.array([1]),z1,z2)),3.+np.concatenate((np.array([0]),basis_func((z1-1)/h-1,3),basis_func((z2-1)/h-2,1))),label="i=5")
plt.plot(np.concatenate((np.array([1]),z2)),3.6+np.concatenate((np.array([0]),basis_func((z2-1)/h-2,2))),label="i=6")
plt.plot(np.concatenate((np.array([1]),z2)),4.2+np.concatenate((np.array([0]),basis_func((z2-1)/h-2,3))),label="i=7")

plt.xlabel('x')
plt.ylabel('$\psi_i(x) + 0.6*i$')
plt.legend()
plt.savefig("p1a")
plt.close()

data=np.loadtxt("./out1b")
data2=np.loadtxt("./out1bb")

fig,ax=plt.subplots()
ax.plot(data.T[0],data.T[2],'o-',label="alternative")
ax.plot(data2.T[0],data2.T[2],'o-',label="original")

ax.legend()
ax.set_xlabel("N")
ax.set_ylabel("L2 error")
ax.set_yscale("log")
ax.set_xscale("log")
plt.savefig("p1b")
plt.close()

fig,ax=plt.subplots()
ax.plot(data.T[3],data.T[2],'o-',label="alternative")
ax.plot(data2.T[3],data2.T[2],'o-',label="original")

ax.legend()
ax.set_xlabel("wall clock time")
ax.set_ylabel("L2 error")
ax.set_yscale("log")
ax.set_xscale("log")
plt.savefig("p1c")
plt.close()




