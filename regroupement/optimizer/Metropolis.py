# -*- coding: utf-8 -*-
"""
Created on 08/12/2021

@author: Yvan GARY
"""

import numpy as np
import random as rd
from regroupement.optimizer.energy_computation import energy_computation
from regroupement.optimizer.Init_alea_G import Init_alea_G

def select_random_debris(G, grps, card_grps):
	''' Function selecting randomly one debris in each group in grps (used to compute neighbours)

	Arguments:
		G (Matrix) : Current state  i.e. current regroupments of debris

		grps (1d-array) : Array containing the groups (of same size) indices

		card_grps (int) : Number of debris contained in each group (the same for every group), typically s_max

	Returns:
		selected_debris (1d-array) : Array containing the indices of the selected debris in each group (same order as groups)  

	'''

	nb_grp = len(grps)
	nb_debris = np.size(G,0)

	selected_debris = []

	for i in range(nb_grp):
		nonzero_idx = np.nonzero(G[:,grps[i]])[0]
		debris_idx = rd.sample(list(nonzero_idx),1)[0]

		selected_debris.append(debris_idx)

	return selected_debris

def split_and_fill(G, grps, grp_idx):
	''' Function selecting randomly a group to split and groups to be filled (used to compute neighbours)

	Arguments:
		G (Matrix) : Current state  i.e. current regroupements of debris

		grps (1d-array) : Array containing the groups (of same size) indices

		grp_idx (1d_array) : indices of the group that can be filled (with cardinal < s_max)

	Returns:
		G (Matrix) : Current state  i.e. current regroupements of debris after distribution

	'''

	nb_grp = np.size(G,1)

	# We select a group to split 
	if (len(grp_idx) == nb_grp):
		# Any group will do
		grp_to_split = rd.randint(0,nb_grp-1)
		grps = np.delete(grps, grp_to_split)
		grp_idx = np.setdiff1d(grp_idx, grp_to_split)

	else:
		# We take a group outside grp_idx
		other_groups = np.setdiff1d(grps, grp_idx)
		grp_to_split = rd.sample(list(other_groups), 1)[0]
		grps = np.delete(grps, grp_to_split)

	# We gather the debris inside the group to split
	debris = np.nonzero(G[:,grp_to_split])[0]
	card_grp = len(debris)

	# We put it into different groups
	grps_to_fill = rd.sample(list(grp_idx), card_grp)
	G[debris, grps_to_fill] = 1

	# We remove the split group from the state
	G = np.delete(G, grp_to_split, axis = 1)

	return G


def Metropolis(G_in, E_in, s_min, s_max, DV, DT, T, V_tol, t_tol, DV_RAAN = None):
	''' Function computing the dynamic of Metropolis. 

	Arguments:
		G_in (Matrix): Current state  i.e. current regroupments of debris

		E_in (float): Energy associated to the current state G_in

		s_min (int) : Minimum number of debris contained in a group

		s_max (int) : Maximum number of debris contained in a group

		DV (Matrix): Matrix containing the delta_v associated to each maneuver

		DT (Matrix): Matrix containing the elapsed time associated to each "J2 perturbation duration" between two debris

		T (float) : Temperature related to the dynamic of Metropolis

		V_tol (float) : Order of magnitude of tolerated dV for one mission in km/s

		t_tol (float) : Order of magnitude of tolerated duration for a mission in days

		DV_RAAN (Matrix) - optionnal : None by default. Matrix containing the delta_v associated to each maneuver, including the ones in RAAN

	Returns:
		G_out (Matrix): Output state of the dynamic of Metropolis 

		E_out (float): Energy associated to the new state G_out

	'''

	nb_debris = np.size(G_in,0)
	nb_grp = np.size(G_in,1)
	G = np.copy(G_in)

	#########################
	# NEIGHBOUR COMPUTATION #
	#########################

	# Defining the neighbour
	if s_min == s_max:
		# Some types of neigbours can not be computed in this case
		case = rd.randint(2,5)
	else:
		case = rd.randint(0,5)

	##############################################################################################
	if case == 0:
		# Rare type of neighbour : we take one debris in s_min different groups of size s_max
		# Then we build a new group of size s_min with them

		# First we get the index of the groups with a cardinal = s_max
		card_grps = sum(G)
		grp_idx = np.where(card_grps == s_max)[0]

		if (len(grp_idx) < s_min):
			# We can't create this neighbour, we create another one
			case = rd.randint(2,5)
		else:
			chosen_grps = rd.sample(list(grp_idx), s_min)	

			# Now we choose one debris in each group
			selected_debris = select_random_debris(G, chosen_grps, s_max)

			# Then we create a new group containing them
			G[selected_debris, chosen_grps] = 0

			new_group = np.zeros((nb_debris,1))
			new_group[selected_debris] = 1

			G = np.hstack((G,new_group))


	##############################################################################################
	if case == 1:
		# Rare type of neighbour : we take a random group and split it into other groups

		# First we get the index of the groups with a cardinal < s_max (groups that can be filled)
		grps = [i for i in range(nb_grp)]
		card_grps = sum(G)
		grp_idx = np.where(card_grps < s_max)[0]

		if (len(grp_idx) < s_max) or (nb_grp <= s_max):
			# We can't always create this neighbour, we create another one
			case = case = rd.randint(2,5)
		else:
			
			G = split_and_fill(G, grps, grp_idx)

	##############################################################################################
	# Preparation for the last cases
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

	##############################################################################################
	if case == 2 or case == 3:
		# We move one debris from a group to another, with following constraints :
		# Min number of debris in a group : s_min
		# Max number of debris in a group : s_max
		if (nb_debris_1 > s_min)*(nb_debris_2 < s_max) :
			G[debris_1,grp1] = 0
			G[debris_1,grp2] = 1
		elif (nb_debris_2 > s_min)*(nb_debris_1 < s_max) :
			G[debris_2,grp2] = 0
			G[debris_2,grp1] = 1
		else :
			# We switch for the always possible case
			case = 4

	##############################################################################################
	if case == 4 or case == 5:
		# Switching the two debris
		G[debris_1,grp1] = 0
		G[debris_2,grp1] = 1

		G[debris_2,grp2] = 0
		G[debris_1,grp2] = 1


	######################
	# ENERGY COMPUTATION #
	######################

	E = energy_computation(G,DV,DT,V_tol,t_tol,DV_RAAN=DV_RAAN)

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



if __name__ == '__main__':
	nb_debris = 6
	s_min = 1
	s_max = 2

	DV = np.array([[0,2,1,3,4,2],[0,0,1,2,5,1],[0,0,0,3,1,5],[0,0,0,0,3,8],[0,0,0,0,0,2],[0,0,0,0,0,0]])

	G_in,E_in = Init_alea_G(nb_debris,s_min,s_max,DV,DT)

	T = 1000 

	# G_out,E_out = Metropolis(G_in, E_in, s_min, s_max, DV, T)

	# print(G_in)
	# print(G_out)











