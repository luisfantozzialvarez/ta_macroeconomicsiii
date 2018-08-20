'''
Markov approximation of AR(1) based on Tauchen's method and simulation
'''
import numpy as np
from scipy.stats import norm

'''
Discretises AR(1) of the form:
y_{t} = \mu.(1-\rho) + \rho.y_{t-1} + u_t
Takes as arguments: rho, mu, sigma_2u (variance of u_t_, r (defaults to 3), and n (number of grids, defaults to 7)
Returns a list with the following objects (in order): grid of possible values and transition matrix
'''
def tauchen(rho, mu, sigma_2u, r = 3, n =7):
    #defines lower and upper bounds
    z1, zn = mu - r*np.sqrt(sigma_2u/(1- rho**2)), mu + r*np.sqrt(sigma_2u/(1- rho**2))
    z = np.linspace(z1, zn, n)
    #generates midpoints
    midpoints = np.empty(n+1)
    midpoints[0] = -np.Inf
    midpoints[n] = np.Inf
    for i in np.arange(1,n,1):
        midpoints[i]= 0.5*(z[i-1]+z[i])
    
    transition = np.empty(shape=(len(z), len(z)))
    
    for i in range(len(z)):
        for j in range(len(z)):
            transition[i,j] = norm.cdf(midpoints[j+1], loc = (1-rho)*mu + rho*z[i], scale = np.sqrt(sigma_2u)) - norm.cdf(midpoints[j], loc = (1-rho)*mu + rho*z[i], scale = np.sqrt(sigma_2u))
    
    return [z, transition]

"""           
simulates N repetitions of Markov process. Takes as arguments:
vector z of possible values; 
transition matrix;
starting value as an index of z.
"""
def markov_sim(z, transition, start, N):
    simul = np.empty(N+1)
    simul[0] = z[start]
    state = start    
    for t in range(N):
        probability = transition[state,:]
        #randomly draws from vector of indices, with probability given by vector probability
        state =  np.random.choice(range(len(z)), p = probability)
        simul[t+1] = z[state]
    return simul
    





            
