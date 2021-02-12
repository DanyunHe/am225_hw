import numpy as np 
import matplotlib.pyplot as plt 

# data=np.fromfile("2.out",dtype='int32',count=-1,sep="")
# data=data.reshape((80,40))

fig,ax=plt.subplots(2,4,figsize=(10,7))
for i in range(7):
	data=np.fromfile("%d.out"%(i*25),dtype='int32',count=-1,sep="")
	data=data.reshape((80,40))

	ax[i//4,i%4].imshow(data)
	ax[i//4,i%4].set_title("time %d"%(i*25))
	ax[i//4,i%4].grid(True,which='both',ls='--')


# ax.set_xticks(np.arange(0,40,1),minor=True)
# ax.set_yticks(np.arange(0,80,1),minor=True)
# ax.grid(True,which='both',ls='--')
plt.savefig("cell_snapshot")


wt=np.array([[0.000003,0.000019,0.000023,0.000028],
			[0.000004,0.000010,0.000013,0.000016],
			[0.000014,0.000021,0.000022,0.000024],
			[0.000097,0.000077,0.000072,0.000069],
			[0.000516,0.000298,0.00023,0.000184],
			[ 0.002179,0.001152,0.000847,0.000660],
			[0.008134,0.004302,0.003118,0.002415]])

plt.close("all")
x=np.arange(4,11,1)
for i in range(3):
	p=wt[:,0]/((i+2)*wt[:,i+1])
	plt.plot(x,p,'o--',label="T=%d"%(i+2))

plt.legend()
plt.xlabel("log2n")
plt.ylabel("efficiency p(n,T)")
plt.title("efficiency vs n")
plt.savefig("2c")
