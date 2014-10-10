# adal.f90 plotter. See adal.f90 for
# matrix sizes and name format. 

import numpy as np
import matplotlib.pyplot as plt


n      = 4
omegas = [.01,.5,1.,5.]

ax1 = plt.subplot(1,1,1)
ax1.set_title('Eigenfunctions')
ax1.set_ylabel('u(rho)')
ax1.set_xlabel(  'rho' )
ax1.grid(True)

for i in range(0,n):
    name = 'sol%2.0f.dat'%(i+11) 
    X = np.transpose(np.loadtxt(name))
    X = X[::-1]
    N = len(X)
    RHO = np.linspace(0.,27.0,N+2)
    XX = [0]; XX.extend(X); XX.append(0)
    ax1.plot( RHO,XX, label='omega = %4.2f' % (omega[i]) )

ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()
