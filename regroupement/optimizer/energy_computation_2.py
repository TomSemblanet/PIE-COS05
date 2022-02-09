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

		for l in range(nb_debris-1):
			dV = DV[debris[l],debris[l+1]]
			dV_transfers.append(dV)
			dV += dV

	if detail:
		return dV, dV_transfers
	else:
		return dV

def energy_computation_DT(G, DT, detail = False):
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
	dT = 0
	dT_transfers = []

	for n in range(nb_grp):

		grp = G[:,n]
		debris = np.nonzero(grp)[0]

		nb_debris = len(debris)

		for l in range(nb_debris-1):
			dT = DT[debris[l],debris[l+1]]
			dT_transfers.append(dT)
			dT += dT

	if detail:
		return dT, dT_transfers
	else:
		return dT


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
