import numpy as np
import matplotlib.pyplot as plt

fsal= np.loadtxt('./out1a')
euler=np.loadtxt('./am225_hw2_files/euler.conv_dat')
heun=np.loadtxt('./am225_hw2_files/heun3.conv_dat')
ralston=np.loadtxt('./am225_hw2_files/ralston.conv_dat')
rk4=np.loadtxt('./am225_hw2_files/rk4.conv_dat')
fig,ax=plt.subplots()
ax.plot(fsal.T[1],fsal.T[0],'.-',label='fsal')
ax.plot(euler.T[1],euler.T[0],'.-',label='euler')
ax.plot(heun.T[1],heun.T[0],'.-',label='heun')
ax.plot(ralston.T[1],ralston.T[0],'.-',label='ralston')
ax.plot(rk4.T[1],rk4.T[0],'.-',label='rk4')
# plt.plot(test.T[0],test.T[2],'.-')
# x = np.linspace(0,8,1000)
ax.invert_xaxis()
# ax.invert_yaxis()
ax.set_xscale("log")
ax.set_yscale("log")
ax.legend()
ax.set_xlabel("Precision")
ax.set_ylabel("Function evaluations")
plt.savefig("1a")
# plt.show()

plt.close("all")

fsal= np.loadtxt('./out1b_do')
fsal2=np.loadtxt('./out1b_int')
fig,ax=plt.subplots()
ax.plot(fsal.T[0],fsal.T[1],'-',label='y1')
ax.plot(fsal.T[0],fsal.T[2],'-',label='y2')
ax.plot(fsal2.T[0],fsal2.T[1],'o')
ax.plot(fsal2.T[0],fsal2.T[2],'o')
# x = np.linspace(0,8,1200)
ax.legend()
plt.savefig("1ba")

plt.close("all")
# x = np.linspace(0,8,1200)
y1=np.zeros(1200)
y2=np.zeros(1200)
# for i in range(1200):
# 	y1[i]=np.cos(x[i]*x[i]/2.)
# 	y2[i]=np.sin(x[i]*x[i]/2.)

fsal= np.loadtxt('./out1b_do')
for i,x in enumerate(fsal.T[0]):
	y1[i]=np.cos(x*x/2.)
	y2[i]=np.sin(x*x/2.)

fig,ax=plt.subplots()
# ax.plot(fsal.T[0],y1,'-',label='y1')
# ax.plot(fsal.T[0],y2,'-',label='y2')

ax.plot(fsal.T[0],fsal.T[1]-y1,'-',label='y1 num-y1 exact')
ax.plot(fsal.T[0],fsal.T[2]-y2,'-',label='y2 num -y2 exact')
plt.savefig("1bb")



plt.close("all")

rk=np.loadtxt('./out7b')

fig,ax=plt.subplots()
ax.plot(rk.T[1],rk.T[0],'.',label='rk')

ax.invert_xaxis()
# ax.invert_yaxis()
ax.set_xscale("log")
ax.set_yscale("log")
ax.legend()
ax.set_xlabel("Precision")
ax.set_ylabel("Function evaluations")
plt.savefig("7b")


plt.close("all")
rk=np.loadtxt('./out7b')
fsal=np.loadtxt('./out7c')
fig,ax=plt.subplots()
ax.plot(rk.T[1],rk.T[0],'.',label='rk')
ax.plot(fsal.T[1],fsal.T[0],'.',label='fsal')

ax.invert_xaxis()
# ax.invert_yaxis()
ax.set_xscale("log")
ax.set_yscale("log")
ax.legend()
ax.set_xlabel("Precision")
ax.set_ylabel("Function evaluations")
plt.savefig("7c")


plt.close("all")
fsal=np.loadtxt('./out')

fig,ax=plt.subplots()
# ax.plot(fsal.T[0],y1,'-',label='y1')
# ax.plot(fsal.T[0],y2,'-',label='y2')

ax.plot(fsal.T[1],fsal.T[2],'-',label='y1 num-y1 exact')
plt.savefig("traj")







