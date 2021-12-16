# -*- coding: utf-8 -*-
"""
Created on 16/12/2021

@author: Yvan GARY
"""

import numpy as np
import random as rd
import matplotlib.pyplot as plt

from Recuit import Recuit
from Metropolis import Metropolis
from Init_alea_G import Init_alea_G
from Gibbs import Gibbs


def Generate_test_DV(nb_debris):
	''' Function used to compute random DV matrices in order to test the Recuit() function.
	


	Arguments:
		nb_debris (int) : Nulber of debris in the given catalogue


	Returns:
		DV (Matrix) : Test matrix containing the delta_v associated to each maneuver

	'''

	DV = np.zeros([nb_debris,nb_debris])

	for i in range(nb_debris):
		for j in range(i+1,nb_debris):

			dv = rd.uniform(0,20)
			DV[i,j] = dv

	return DV

if __name__ == "__main__":


	nb_debris = 19
	card_grp = 5

	DV = Generate_test_DV(nb_debris)

	########################
	# TESTS FOR METROPOLIS #
	########################

	# n_classes = 100
	# t_iter = 100
	# n_iter = 1

	# G_in, E_in = Init_alea_G(nb_debris,card_grp,DV)

	# G, E = Metropolis(G_in, E_in, DV,T)

	###################
	# TESTS FOR GIBBS #
	###################

	# n_classes = 200
	# t_iter = 200
	# n_iter = 50

	# list_T = [1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001]

	# for T in list_T:

	# 	G, E, freqs = Gibbs(nb_debris, card_grp, DV, T, n_classes, t_iter, n_iter)
	# 	E_mean = np.argmax(freqs)
	# 	print('For T = ', T, ' we get E_min_mean ~ ', E_mean)
	# 	print()



	####################
	# TESTS FOR RECUIT #
	####################

	# Final and Initial Temperature
	Ti = 1
	Tf = 0.1

	alpha = 0.95

	n_classes = 200
	t_iter = 200
	n_iter = 10

	G, E = Recuit(nb_debris, card_grp, DV, Ti, Tf, alpha, n_classes, t_iter, n_iter)






