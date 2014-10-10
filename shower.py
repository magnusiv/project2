import numpy as np
X = np.transpose(np.loadtxt('sol.dat'))
lam = np.sort(X)
print lam[0],' and ',lam[1],' and ',lam[2], ' and ', lam[3], ' and ', lam[4], ' and ', lam[5]
print '_______________'

"""
n = 400, rho max = 7
[magnusiv@sothi proj2]$ python shower.py
2.99990429358  and  6.99952145103  and  10.9988322969  and  14.9978367961  and  18.9965349187  and  22.9949269014
_______________


n=500,rho_max = 7
[magnusiv@sothi proj2]$ rm *.o && rm executable && rm *.dat && gfortran -c jacobi.f90  && gfortran -o executable jacobi.o && time ./executable && python shower.py
 Done!
      415086

real    9m9.168s
user    9m8.923s
sys     0m0.003s
2.99993874865  and  6.99969373634  and  10.9992526988  and  14.9986156217  and  18.9977824962  and  22.9967535899
_______________

"""
