# -*- coding: utf-8 -*-
"""
Created on 09/12/2021

@author: Yvan GARY
"""

import numpy as np

def energy_computation(G, DV, detail = False):
	''' Function used to compute the energy associated to a state

	Arguments:
		G (Matrix): Actual state for which we compute the energy
		DV (Matrix): Matrix containing the delta_v associated to each maneuver

	Returns: 
		E (float): Energy associated to the state G
		E_transfers (array) - optionnal : Energy asociated to each individual group

	'''

	non_zero = np.nonzero(np.transpose(G))
	grp_label = non_zero[0]
	debris_label = non_zero[1]

	# Getting the number of debris
	nb_debris = np.size(debris_label)
	# Getting the number of groups
	nb_grp = np.size(G,1)

	E_transfers = np.array([])

	count = 0

	for i in range(nb_grp):

		e = 0	# Energy associated to a group of debris
		actual_debris = debris_label[count]

		while grp_label[count] == i:
			target_debris = debris_label[count]
			e+=DV[actual_debris,target_debris]

			actual_debris = target_debris
			count+=1

			# Ensuring side effects
			if count == nb_debris:
				break 

		

		E_transfers = np.append(E_transfers,e)

	E = sum(E_transfers)

	if detail == False :
		return E
	else:
		return E, E_transfers
