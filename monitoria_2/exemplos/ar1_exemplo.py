from Markov import tauchen, markov_sim, ar1
import numpy as np
import matplotlib.pyplot as mp

#We built the required functions in (A) and stored them in Markov.py

#Question B

#N=3

#We set a seed to allow for replication
np.random.seed(93)

#generates approximation
markov_approx = tauchen(0.8, 0, 0.01, r = 3, n=3)

#prints values in the grid
print('Grid when n=3:', markov_approx[0])

#index 1 equals 0, so we set it as a starting index in our simulation
simul1 = markov_sim(markov_approx[0], markov_approx[1], 1, 1000)

#resets seed to generate an AR(1) with the same RNG
np.random.seed(93)
ar_sim = ar1(0.8, 0, 0.01, 0, 1000)

#plots results
mp.plot(range(1001), ar_sim, 'y--', range(1001), simul1, 'b-')
mp.legend(['AR(1)', 'Markov approximation'], loc = 'best')
mp.title("N=3")
mp.show()


#N=7

#We set a seed to allow for replication
np.random.seed(93)

#generates approximation
markov_approx = tauchen(0.8, 0, 0.01, r = 3, n=7)
#prints values in grid
print('Grid when n=7:', markov_approx[0])

#index 3 equals 0, so we set it as a starting index in our simulation
simul1 = markov_sim(markov_approx[0], markov_approx[1], 3, 1000)

#plots results
mp.plot(range(1001), ar_sim, 'y--', range(1001), simul1, 'b-')
mp.legend(['AR(1)', 'Markov approximation'], loc = 'best')
mp.title("N=7")
mp.show()



#N=15

#We set a seed to allow for replication
np.random.seed(93)

#generates approximation
markov_approx = tauchen(0.8, 0, 0.01, r = 3, n=15)

#prints values in grid
print('Grid when n=15:', markov_approx[0])

#index 7 equals 0, so we set it as a starting index in our simulation
simul1 = markov_sim(markov_approx[0], markov_approx[1], 7, 1000)

#plots results
mp.plot(range(1001), ar_sim, 'y--', range(1001), simul1, 'b-')
mp.legend(['AR(1)', 'Markov approximation'], loc = 'best')
mp.title("N=15")
mp.show()



