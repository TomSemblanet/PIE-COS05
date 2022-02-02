# -*- coding: utf-8 -*-
"""
Created on 08/12/2021

@author: Yvan GARY
"""

import numpy as np
import random as rd
from regroupement.optimizer.energy_computation import energy_computation

def Init_alea_G(nb_debris, s_min, s_max, DV):
	''' Function used to initiate the optimization

	Arguments:
		nb_debris (int): Number of debris in the given catalogue
		s_min (int) : Minimum number of debris contained in a group
		s_max (int) : Maximum number of debris contained in a group
		DV (Matrix): Matrix containing the delta_v associated to each maneuver

	Returns:
		G (matrix): First state generated randomly to begin Optimization 
		E (float): Energy associated to the state G
	'''

	debris = np.array([i for i in range(nb_debris)])
	nb_debris_init = nb_debris

	G = np.zeros((1,nb_debris_init));
	i = 0

	while nb_debris >= s_min:
		# Choosing randomly cardinal of the group between min and max values
		if nb_debris >= s_max:
			rand = rd.randint(0,1)
			card_grp = rand*s_min + (1-rand)*s_max
		else:
			card_grp = s_min

		ind = rd.sample(range(nb_debris), card_grp)	
		selection = debris[ind]

		if i == 0:
			G[i,selection] = 1

		else:
			new_grp = np.zeros(nb_debris_init)
			new_grp[selection] = 1

			# Actualizing the state and the remaining debris
			G = np.vstack((G,new_grp))


		debris = np.delete(debris,ind)
		nb_debris = np.size(debris)
		i+=1;

	G = np.transpose(G)

	#################################################
	#           Computation of the Energy           #
	#################################################

	E = energy_computation(G,DV)

	#################################################

	return G,E



if __name__ == "__main__":
	nb_debris = 6
	s_min = 4
	s_max = 5

	DV = np.array([[0,2,1,3,4,2],[0,0,1,2,5,1],[0,0,0,3,1,5],[0,0,0,0,3,8],[0,0,0,0,0,2],[0,0,0,0,0,0]])

	G,E = Init_alea_G(nb_debris,s_min,s_max,DV)

	print(sum(G))
	print(G)
	print(E)

	# For unitary tests, check that sum(sum(G)) = nb_debris and sum(G) = [card_grp, card_grp, card_grp, ... , card_grp - nb_debris%card_grp]




