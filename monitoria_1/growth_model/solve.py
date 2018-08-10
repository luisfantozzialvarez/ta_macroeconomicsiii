#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  9 20:22:44 2018

@author: luisfantozzialvarez
"""

from DiscreteGrowthCapital import DiscreteGrowthCapital
import numpy as np
import matplotlib.pyplot as plt

#Creating an economy
economia1 = DiscreteGrowthCapital(1,1/3,0.3,0.96, 2, 1)

#Resetting k0 to 0.3*kss
economia1.k0 = 0.3*economia1.kss

#Let's solve the model assuming that after T = 100, economy reaches steady-state
#Initial guess for optimal path is linear from k0 to kss
k_guess = np.linspace(economia1.k0, economia1.kss, 102)[1:-1]

#Solving model
k_optimal = economia1.solve_model(k_guess)

#Adds k0 and kss to path
k_optimal = np.concatenate(([economia1.k0], k_optimal, [economia1.kss,economia1.kss]))
#Computes consumption path
c_optimal  = economia1.compute_consumption(k_optimal[:-1],k_optimal[1:])

#Plotting results
plt.plot(np.arange(0,103,1), k_optimal,'-', color = 'red', label = 'Capital')
plt.plot(np.arange(0,102,1), c_optimal,'--', color = 'blue', label = 'Consumption')
plt.legend()
plt.grid()
plt.savefig('exemplo1.pdf')
plt.show()