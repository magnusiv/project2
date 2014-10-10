# adal.f90 plotter. See adal.f90 for
# matrix sizes and name format. 

import numpy as np
import matplotlib.pyplot as plt


n      = 4

ax1 = plt.subplot(1,2,1)
ax1.set_title('Solution curves')
ax1.set_ylabel('$u(x)$')
ax1.set_xlabel(  '$x$' )
ax1.grid(True)
"""
ax2 = plt.subplot(1,2,2)
ax2.set_title('Relative error')
ax2.set_ylabel('Sup$_i \; \epsilon_i $ ')
ax2.set_xlabel(  'log$_2 h$' )
ax2.grid(True)
"""
for i in range(0,n):
    name = 'sol%2.0f.dat'%(i+11) 
    X = np.transpose(np.loadtxt(name))
    X = X[::-1]
    N = len(X)
    RHO = np.linspace(0.,27.0,N+2)
    XX = [0]; XX.extend(X); XX.append(0)
    ax1.plot( RHO,XX, label='l = %4.0f' % (i+1) )

ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#ax2.plot(steps, errors, '-ro')
plt.show()
