# -*- coding: utf-8 -*-
"""
Created on 08/12/2021

@author: Yvan GARY
"""

import numpy as np
import random as rd
from regroupement.optimizer.energy_computation import energy_computation
from regroupement.optimizer.Init_alea_G import Init_alea_G

def Metropolis(G_in, E_in, DV,T):
	''' Function computing the dynamic of Metropolis. A neighbour of a state G is defined
		as a switch of two debris between two groups selected randomly. then it is kept or
		abandonned according to the Metropolis dynamic.

	Arguments:
		G_in (Matrix): Current state  i.e. current regroupments of debris
		E_in (float): Energy associated to the current state G_in
		DV (Matrix): Matrix containing the delta_v associated to each maneuver
		T (float) : Temperature related to the dynamic of Metropolis

	Returns:
		G_out (Matrix): Output state of the dynamic of Metropolis 
		E_out (float): Energy associated to the new state G_out

	'''

	nb_debris = np.size(G_in,0)
	nb_grp = np.size(G_in,1)
	G = np.copy(G_in)

	# # Max and min sizes of groups
	s_M = 5
	s_m = 4

	#########################
	# NEIGHBOUR COMPUTATION #
	#########################

	# Selecting two groups randomly
	groups = rd.sample(range(nb_grp), 2)
	grp1 = groups[0]
	grp2 = groups[1]

	# Selecting two random debris in these groups
	debris_labels_1 = np.nonzero(G_in[:,grp1])
	debris_labels_2 = np.nonzero(G_in[:,grp2])

	nb_debris_1 = np.size(debris_labels_1)
	nb_debris_2 = np.size(debris_labels_2)

	idx1 = np.random.randint(nb_debris_1)
	idx2 = np.random.randint(nb_debris_2)

	debris_1 = debris_labels_1[0][idx1]
	debris_2 = debris_labels_2[0][idx2]

	# Defining the neighbour
	case = rd.randint(0,1)

	if case == 0:
		# We move one debris from a group to another, with following constraints :
		# Min number of debris in a group : 4
		# Max number of debris in a group : 6
		if (nb_debris_1 > s_m)*(nb_debris_2 < s_M) :
			G[debris_1,grp1] = 0
			G[debris_1,grp2] = 1
		elif (nb_debris_2 > s_m)*(nb_debris_1 < s_M) :
			G[debris_2,grp2] = 0
			G[debris_2,grp1] = 1
		else :
			case = 1

	if case == 1:
		# Switching the two debris
		G[debris_1,grp1] = 0
		G[debris_2,grp1] = 1

		G[debris_2,grp2] = 0
		G[debris_1,grp2] = 1

	######################
	# ENERGY COMPUTATION #
	######################

	E = energy_computation(G,DV)

	###################
	# STATE SELECTION #
	###################

	p = rd.uniform(0,1)

	if E < E_in - T*np.log(p):
		G_out = G
		E_out = E
	else:
		G_out = G_in	
		E_out = E_in

	return G_out, E_out














