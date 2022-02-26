# -*- coding: utf-8 -*-
"""
Created on 13/12/2021

@author: Yvan GARY
"""

import numpy as np
import random as rd
import matplotlib.pyplot as plt
from time import sleep
from tqdm import tqdm

from regroupement.optimizer.Init_alea_G import Init_alea_G
from regroupement.optimizer.Metropolis import Metropolis

def Recuit(nb_debris, s_min, s_max, DV, DT, Ti, Tf, alpha, n_classes, t_iter, n_iter, cost = 'only_dV'):
	''' Function computing the simulated annealing, with the corresponding dynamic of Metropolis.
		It corresponds to the succession of Markov chains computed with decreasing temperatures. At the end
		we obtain a state G_out that minimizes the energy we defined, that is to say the sum of the delta_v
		associated to each groups. G_out contains the final groups that reach this minimal "global delta_v".



	Arguments:
		nb_debris (int) : Number of debris in the given catalogue
		
		s_min (int) : Minimum number of debris contained in a group
		
		s_max (int) : Maximum number of debris contained in a group
		
		DV (Matrix) : Matrix containing the delta_v associated to each maneuver
		
		DT (Matrix): Matrix containing the elapsed time associated to each "J2 perturbation duration" between two debris
		
		Ti (float) : Initial temperature related to the dynamic of Metropolis
		
		Tf (float) : Final temperature related to the dynamic of Metropolis
		
		alpha (float) : Geometric factor to decrease Temperature (0 < alpha < 1)
		
		n_classes (array) : Number of classes for the displayed histogram (ex : range(100))
		
		t_iter (int) : Number of iterations for a Markov chain
		
		n_iter (int) : Number of Markov chains generated for each Temperature

		cost (string) : Two possible options 
			--> 'only_dV' (default parameter) : Selects cost function taking in account only the cost of dV maneuvers
			--> 'dV_and_drift' : Selectes cost function taking in account the cost of dV maneuvers and the time of the drift due to the J2 perturbation


	Returns:
		G_out (matrix) : Output state of the dynamic of Metropolis 
		
		E_out (float) : Energy associated to the new state G_out
		
		freqs (array) : Frequencies associated to each energy

	'''

	E_evol = np.zeros(t_iter) 	# Evolution of Energy along Markov chain
	freqs = np.zeros(n_classes)

	# Initializing Temperature
	T = Ti

	# Number of iterations over T
	r = np.ceil(np.log(Tf/Ti)/np.log(alpha))

	# Number of values for Energy (at each T) we keep for final display
	k = int(np.floor(t_iter/4))

	# Evolution of Energy along the whole simulated annealing process
	E_evol_global = np.zeros(int(r*k))

	# Initialization of an energy for the histogram range
	G_out,E_out = Init_alea_G(nb_debris, s_min, s_max, DV, DT, cost = cost)
	E_max = int(E_out)

	count = 1

	print('\n#########')
	print('ITERATION')
	print('#########\n')

	print('\n')

	while T > Tf:
	# for t in tqdm(range(np.int(r))):

		# sleep(3)

		print(count, '/', np.int(r))

		for i in range(n_iter):

			if T == Ti:
				# Initialization of a random state 
				G_out,E_out = Init_alea_G(nb_debris, s_min, s_max, DV, DT, cost = cost)
				E_evol[0] = E_out

			# Transitory Markov Chain to reach a minimum
			for t in range(1,t_iter):
				G_in = G_out
				E_in = E_out
				G_out, E_out = Metropolis(G_in, E_in, s_min, s_max, DV, DT, T, cost = cost)
				E_evol[t] = E_out

			# Beginning of the recorded trajectory
			E_evol[0] = E_evol[-1]

			for t in range(1,t_iter):
				G_in = G_out
				E_in = E_out
				G_out, E_out = Metropolis(G_in, E_in, s_min, s_max, DV, DT, T, cost = cost)
				E_evol[t] = E_out

			# Getting the last portion of the chain for this T (for the plot)
			if i == n_iter-1 :
				E_evol_global[k*(count-1):k*count] = E_evol[t_iter-k-1:-1]

			frequency, bins = np.histogram(E_evol, bins = n_classes, range = (0,n_classes))
			freqs += frequency


		T = alpha*T
		count += 1

	# Plotting the histogram
	x_hist = np.linspace(0, n_classes, n_classes)
	plt.plot(x_hist, freqs)
	plt.title('Frequency of Energies')
	plt.xlabel('Energy')
	plt.ylabel('Frequency')
	plt.show()

	x_E = [k for k in range(int(r*k))]

	plt.figure
	plt.plot(x_E, E_evol_global)
	plt.title('Evolution of Energy along the Simulated Annealing')
	plt.show()

	return G_out, E_out, freqs













