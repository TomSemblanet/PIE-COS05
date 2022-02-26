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

	nb_grp = np.size(G,1)
	dV = 0
	dV_transfers = []

	for n in range(nb_grp):

		grp = G[:,n]
		debris = np.nonzero(grp)[0]

		nb_debris = len(debris)

		dV_temp = 0
		for l in range(nb_debris-1):
			dv = DV[debris[l],debris[l+1]]
			dV_temp += dv
		
		dV_transfers.append(dV_temp)
		dV += dV_temp

	dV_transfers = np.array(dV_transfers)

	if detail:
		return dV, dV_transfers
	else:
		return dV

def energy_computation_DT(G, DT, detail = False):
	''' Function used to compute the delta t (J2) associated to a state

	Arguments:
		G (Matrix): Actual state for which we compute the energy
		
		DV (Matrix): Matrix containing the delta_v associated to each maneuver
		
		detail (bool) : False by default - If True, gives the detail of the delta v for each group

	Returns: 
		dV (float): Global delta v associated to the state G
		
		dV_transfers (array) - optionnal : Global delta v associated to each individual group

	'''

	nb_grp = np.size(G,1)
	dT = 0
	dT_transfers = []

	for n in range(nb_grp):

		grp = G[:,n]
		debris = np.nonzero(grp)[0]

		nb_debris = len(debris)

		dT_temp = 0
		for l in range(nb_debris-1):
			dt = DT[debris[l],debris[l+1]]
			dT_temp += dt

		dT += dT_temp
		dT_transfers.append(dT_temp)

	dT_transfers = np.array(dT_transfers)

	if detail:
		return dT, dT_transfers
	else:
		return dT


def energy_computation(G, DV, DT, detail = False, cost = 'only_dV'):
	''' Function used to compute the energy associated to a state

	Arguments:
		G (Matrix): Actual state for which we compute the energy
		
		DV (Matrix): Matrix containing the delta_v associated to each maneuver
		
		DT (Matrix): Matrix containing the elapsed time associated to each "J2 perturbation duration" between two debris
		
		detail (bool) : False by default - If True, gives the detail of the energy for each group

		cost (string) : Two possible options 
			--> 'only_dV' (default parameter) : Selects cost function taking in account only the cost of dV maneuvers
			--> 'dV_and_drift' : Selectes cost function taking in account the cost of dV maneuvers and the time of the drift due to the J2 perturbation

	Returns: 
		E (float): Energy associated to the state G
		
		E_transfers (array) - optionnal : Energy associated to each individual group
		
		E_transfers_dV (array) - optionnal : delta v associated to each individual group
		
		E_transfers_dt (array) - optionnal : elapsed time due to J2 perturbation associated to each individual group in cost = 'dV_and_drift' case

	'''

	# Computation of the dV part
	dV, dV_transfers = energy_computation_DV(G, DV, detail = True)

	if cost == 'only_dV':

		# Computation of the energy
		E = dV
		E_transfers = dV_transfers

	elif cost == 'dV_and_drift':

		# Computation of the dt part
		dt, dt_transfers = energy_computation_DT(G, DT, detail = True)

		# Weighting the two parts
		V_tol = 1.5 # [km/s]
		t_tol = 3.0*365.0 # [days]

		weight_coef = V_tol/t_tol

		# Computation of the energy
		E = dV + weight_coef*dt
		E_transfers = dV_transfers + weight_coef*dt_transfers


	else:
		raise Exception('This cost function is not implemented, please choose between only_dV or dV_and_drift')


	if detail:
		if cost == 'only_dV':
			return E, E_transfers, dV_transfers
		elif cost == 'dV_and_drift':
			return E, E_transfers, dV_transfers, dt_transfers
	else:
		return E
