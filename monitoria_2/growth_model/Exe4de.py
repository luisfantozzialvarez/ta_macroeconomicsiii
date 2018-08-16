from Growth import Growth
import matplotlib.pyplot as mp
import numpy as np
from scipy.optimize import fsolve

#We constructed a class named Growth, stored in Growth.py, in order to solve this exercise.
#This class generates an economy, given our choice of parameters, and allows us
#to compute steady state values and solve the model
#Arguments passed in the initialization are, in the following order:
#beta, delta, gamma, eta, alpha, theta, tau_c, tau_h, tau_k, k0.
#Moreover, since it was not specified, we set population growth to 0
economia1 = Growth(0.98, 0.08, 0.015, 0, 0.4, 2., 0.15, 0.25, 0.15, 4)

#Gets steady-state value of capital
k_ss = economia1.get_steady_state()[1]

#changes initial capital stock to 0.8*k_ss
economia1.__setattr__('k0', 0.8*k_ss)

#solves the model using 200 periods
solution = economia1.solve_model(200)
t_axis = np.linspace(0, 201, 202)

#Prepares to plot results
mp.figure(1, dpi=125)
mp.subplot(321)


mp.plot(t_axis, solution[0], 'b-')
mp.grid()
mp.title('Hours worked')

mp.subplot(322)
mp.plot(t_axis, solution[1], 'b-')
mp.grid()
mp.title('Capital per effective labour')

mp.subplot(323)
mp.plot(t_axis, solution[2], 'b-')
mp.grid()
mp.title('Consumption per effective labour')

mp.subplot(324)
mp.plot(t_axis, solution[3], 'b-')
mp.grid()
mp.title('Output per effective labour')

mp.subplot(325)
mp.plot(t_axis, solution[4], 'b-')
mp.grid()
mp.title('Investment per effective labour')

mp.subplot(326)
mp.plot(t_axis, solution[5], 'b-')
mp.grid()
mp.title('Gvt. spending per effective labour')

#Figure title
mp.suptitle('Initial Economy - Transition to Steady State', fontsize = 10)
#this option tries to position subplots so they do not overlap. We also leave some space for the plot title
mp.tight_layout(rect=[0,0,1,0.97])

#shows plot
mp.show()
