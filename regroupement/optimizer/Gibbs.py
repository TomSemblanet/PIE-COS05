# -*- coding: utf-8 -*-
"""
Created on 16/12/2021

@author: Yvan GARY
"""

import numpy as np
import random as rd
import matplotlib.pyplot as plt

from regroupement.optimizer.energy_computation import energy_computation
from regroupement.optimizer.Init_alea_G import Init_alea_G
from regroupement.optimizer.Metropolis import Metropolis


def Gibbs(nb_debris, card_grp, DV, T, n_classes, t_iter, n_iter):
	''' Function propagating the dynamic of Metropolis along a Markov chain for a given Temperature n_iter times.

	Arguments:
		nb_debris (int) : Nulber of debris in the given catalogue
		card_group (int): Cardinal of a group generated from this debris
		DV (Matrix) : Matrix containing the delta_v associated to each maneuver
		T (float) : Temperature related to the dynamic of Metropolis
		n_classes (array) : Number of classes for the displayed histogram (ex : range(100))
		t_iter (int) : Number of iterations for a Markov chain
		n_iter (int) : Number of Markov chains generated for each Temperature


	Returns:
		G_out (matrix) : Output state of the dynamic of Metropolis 
		E_out (float) : Energy associated to the new state G_out
		freqs (array) : Array containing the frequencies associated to each energy

	'''

	E_evol = np.zeros(t_iter) 	# Evolution of Energy along Markov chain
	freqs = np.zeros(n_classes)

	for i in range(n_iter):

		# Initialization of a random state 
		G_out,E_out = Init_alea_G(nb_debris, card_grp, DV)
		E_evol[0] = E_out

		# Transitory Markov Chain to reach a minimum
		for t in range(1,t_iter):
			G_in = G_out
			E_in = E_out
			G_out, E_out = Metropolis(G_in, E_in, DV,T)
			E_evol[t] = E_out

		# Beginning of the recorded trajectory
		E_evol[0] = E_evol[-1]

		for t in range(1,t_iter):
			G_in = G_out
			E_in = E_out
			G_out, E_out = Metropolis(G_in, E_in, DV,T)
			E_evol[t] = E_out

		# frequency, bins, patches = plt.hist(E_evol, bins = n_classes, range= (0,n_classes))
		frequency, bins = np.histogram(E_evol, bins = n_classes, range = (0,n_classes))
		freqs += frequency
	
	# Plotting the histogram	
	plt.figure
	x_hist = np.arange(n_classes)
	plt.plot(x_hist,freqs)
	plt.title('Gibbs Distribution for T = %f' %T)
	plt.xlabel('Energy')
	plt.ylabel('Frequency')
	plt.show()


	# Plotting the evolution of the energy
	# x_E = np.arange(t_iter)

	# plt.figure
	# plt.plot(x_E, E_evol)
	# plt.title('Evolution of Energy along the Chain')
	# plt.show()

	return G_out, E_out, freqs






