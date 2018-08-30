#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 18:13:03 2018

@author: thiagoyuditachibana
"""

from GrowthIID import GrowthIID
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib import cm

zgrid = np.array([1,2,3])
economia1 =  GrowthIID(2, 0.96, 0.3, 0.25, np.array([1/2,1/4,1/4]), \
                       zgrid, 17)

kgrid =  np.linspace(0.5*economia1.kss, 1.5*economia1.kss, 300)

economia1.reconstr_grid(kgrid)

#Grids with k replicated across columns and z replicated across lines

gridK, gridZ = np.meshgrid(kgrid, zgrid, indexing = 'ij')

#Initial guess for V0: consume all capital every period forever
V0 = (1/(1-economia1.beta))*economia1.utility(gridK)

V0, kprime, index = economia1.value_iteration(V0, 1e-5)

#Capital policy function
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(gridK, gridZ, kprime, cmap = cm.coolwarm)
ax.set_xlabel('$k_t$')
ax.set_ylabel('$z_t$')
ax.set_zlabel('$k^\\prime(k,h)$')
plt.savefig('capitalpolicy.pdf')
plt.show()


#Policy iteration
V0_policy, kprime_policy, index_policy = economia1.policy_iteration(V0, 1e-5, 10)

#Capital policy function
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(gridK, gridZ, kprime_policy, cmap = cm.coolwarm)
ax.set_xlabel('$k_t$')
ax.set_ylabel('$z_t$')
ax.set_zlabel('$k^\\prime(k,h)$')
plt.savefig('capitalpolicy_policyiteration.pdf')
plt.show()