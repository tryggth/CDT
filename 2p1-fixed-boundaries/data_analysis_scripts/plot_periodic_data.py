#!/usr/bin/env python2

import numpy as np
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
import average_ensemble as ae

dataname = '/home/jonah/Dropbox/Summer_Quantum_Gravity/code/' \
    + 'volumes_aligned_table_format.dat'
print(dataname)
N3 = 10000
data = np.loadtxt(dataname)

means = np.mean(data,0)
T_slices = len(means)
stds = np.std(data,0)

popt,pcov = ae.fit_to_data(means,N3)
xdata = np.arange(0.0,T_slices,1.0)

fit = ae.de_sitter_form_testing(xdata,popt[0],popt[1],popt[2],N3,T_slices)

lines = [plt.errorbar(xdata,means,yerr=stds,fmt='ro',label='Monte Carlo'),
         plt.plot(xdata,fit,'b-',label='Fit')]

plt.xlabel('Proper Time')
plt.ylabel('Spatial Extent')
plt.legend()
#plt.title('Monte Carlo Simulation With Periodic Boundary Conditions')

plt.show()
