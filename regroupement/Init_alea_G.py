# -*- coding: utf-8 -*-
"""
Created on 08/12/2021

@author: Yvan GARY
"""

import numpy as np
import random as rd
from energy_computation import energy_computation

def Init_alea_G(nb_debris, card_grp, DV):
	''' Function used to initiate the optimization

	Arguments:
		nb_debris (int): Number of debris in the given catalogue
		card_group (int): Cardinal of a group generated from this debris
		DV (Matrix): Matrix containing the delta_v associated to each maneuver

	Returns:
		G (matrix): First state generated randomly to begin Optimization 
		E (float): Energy associated to the state G

	'''

	debris = np.array([i for i in range(nb_debris)])

	if nb_debris%card_grp == 0:
		nb_grp = np.floor(nb_debris/card_grp)
	else:
		nb_grp = np.floor(nb_debris/card_grp) + 1

	nb_grp = int(nb_grp)

	# Initializing G matrix
	G = np.zeros([nb_debris, nb_grp])


	for i in range(nb_grp-1):
		ind1 = rd.sample(range(nb_debris), card_grp)
		selection = debris[ind1]
		G[selection,i] = 1

		# Actualizing the remaining debris
		debris = np.delete(debris,ind1)
		nb_debris = np.size(debris)


	# Creating last group manually
	selection = debris
	G[selection,-1] = 1

	#################################################
	#           Computation of the Energy           #
	#################################################

	E = energy_computation(G,DV)

	#################################################

	return G,E



if __name__ == "__main__":
	nb_debris = 5
	card_grp = 2

	DV = np.array([[0,2,1,3,4],[0,0,1,2,5],[0,0,0,3,1],[0,0,0,0,3],[0,0,0,0,0]])

	G,E = Init_alea_G(nb_debris,card_grp,DV)

	print(G)
	print(E)

	# For unitary tests, check that sum(sum(G)) = nb_debris and sum(G) = [card_grp, card_grp, card_grp, ... , card_grp - nb_debris%card_grp]




