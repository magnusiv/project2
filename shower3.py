import numpy as np
import matplotlib.pyplot as plt


n      = 4
omegas = [.01,.5,1.,5.]
end = 7.

ax1 = plt.subplot(1,1,1)
ax1.set_title('Eigenfunctions')
ax1.set_ylabel('|R(rho)|^2')
ax1.set_xlabel(  'rho' )
ax1.grid(True)

for i in range(0,n):
    name = 'sol%2.0f.dat'%(i+11) 
    X = np.transpose(np.loadtxt(name))
    X = X[::-1]
    N = len(X)
    RHO = np.linspace(end/(N+2),end*(1.-1./(N+2)),N)
    for j in range(0,N):
	X[j] = (X[j]/RHO[j])**2
    X = X/np.linalg.norm(X)
    RHO = np.linspace(end/(N+2),end,N+1)
    XX = []; XX.extend(X); XX.append(0)
    ax1.plot( RHO,XX, label='omega = %4.2f' % (omegas[i]) )

ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()
