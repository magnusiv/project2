import numpy as np
import matplotlib.pyplot as plt

Y = np.transpose(np.loadtxt('sol2.dat'))
X = np.transpose(np.loadtxt('omegas.dat'))

ax1 = plt.subplot(1,1,1)
ax1.set_title('Ground state energy eigenvalue')
ax1.set_ylabel('$1/\omega$')
ax1.set_xlabel(  '$\lambda`$' )
ax1.grid(True)
ax1.plot(X,Y)
plt.show()
