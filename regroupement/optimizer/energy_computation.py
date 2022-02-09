# -*- coding: utf-8 -*-
"""
Created on 09/12/2021

@author: Yvan GARY
"""

import numpy as np

def energy_computation_DV(G, DV, detail = False):
	''' Function used to compute the delta v associated to a state

	Arguments:
		G (Matrix): Actual state for which we compute the energy
		DV (Matrix): Matrix containing the delta_v associated to each maneuver
		detail (bool) : False by default - If True, gives the detail of the delta v for each group

	Returns: 
		dV (float): Global delta v associated to the state G
		dV_transfers (array) - optionnal : Global delta v associated to each individual group

	'''

	non_zero = np.nonzero(np.transpose(G))
	grp_label = non_zero[0]
	debris_label = non_zero[1]

	# Getting the number of debris
	nb_debris = np.size(debris_label)
	# Getting the number of groups
	nb_grp = np.size(G,1)

	dV_transfers = np.array([])

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

		

		dV_transfers = np.append(dV_transfers,e)

	dV = sum(dV_transfers)

	if detail == False :
		return dV
	else:
		return dV, dV_transfers

def energy_computation_DT(G, DT, detail = False):
	''' Function used to compute the global elapsed time when using J2 perturbation instead of maneuvers

	Arguments:
		G (Matrix): Actual state for which we compute the energy
		DT (Matrix): Matrix containing the elapsed time associated to each "J2 perturbation duration" between two debris
		detail (bool) : False by default - If True, gives the detail of the elapsed time dt for each group

	Returns: 
		dt (float): Global elapsed time associated to the state G
		dt_transfers (array) - optionnal : Global elapsed time associated to each individual group

	'''

	non_zero = np.nonzero(np.transpose(G))
	grp_label = non_zero[0]
	debris_label = non_zero[1]

	# Getting the number of debris
	nb_debris = np.size(debris_label)
	# Getting the number of groups
	nb_grp = np.size(G,1)

	dt_transfers = np.array([])

	count = 0

	for i in range(nb_grp):

		e = 0	# Energy associated to a group of debris
		actual_debris = debris_label[count]

		while grp_label[count] == i:
			target_debris = debris_label[count]
			e+=DT[actual_debris,target_debris]

			actual_debris = target_debris
			count+=1

			# Ensuring side effects
			if count == nb_debris:
				break 

		

		dt_transfers = np.append(dt_transfers,e)

	dt = sum(dt_transfers)

	if detail == False :
		return dt
	else:
		return dt, dt_transfers

def energy_computation(G, DV, DT, detail = False):
	''' Function used to compute the energy associated to a state

	Arguments:
		G (Matrix): Actual state for which we compute the energy
		DV (Matrix): Matrix containing the delta_v associated to each maneuver
		DT (Matrix): Matrix containing the elapsed time associated to each "J2 perturbation duration" between two debris
		detail (bool) : False by default - If True, gives the detail of the energy for each group

	Returns: 
		E (float): Energy associated to the state G
		E_transfers (array) - optionnal : Energy associated to each individual group
		E_transfers_dV (array) - optionnal : delta v associated to each individual group
		E_transfers_dt (array) - optionnal : elapsed time due to J2 perturbation associated to each individual group

	'''

	# Computation of the dV part
	dV, dV_transfers = energy_computation_DV(G, DV, detail = True)

	# Computation of the dt part
	dt, dt_transfers = energy_computation_DT(G, DT, detail = True)

	# Weighting the two parts
	max_tolerated_dV = 3.0 # [km/s]
	max_tolerated_dt = 365 # [days]

	weight_coef = max_tolerated_dV/max_tolerated_dt

	# Computation of the energy
	E = dV + weight_coef*dt
	E_transfers = dV_transfers + weight_coef*dt_transfers

	if detail:
		return E, E_transfers, dV_transfers, dt_transfers
	else:
		return E

















