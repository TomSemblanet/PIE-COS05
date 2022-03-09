# -*- coding: utf-8 -*-
"""
Created on 09/12/2021

@author: Yvan GARY
"""

import numpy as np
from itertools import permutations
from regroupement.dV_computations.compute_dt_alignment import compute_dt, RAAN_evol

RAAN_evol = np.vectorize(RAAN_evol)

def single_energy_computation(debris, DV, debris_data, DV_RAAN = None, V_tol = None, t_tol = None):
	''' Function used to compute the delta t (J2) associated to a state

	Arguments:
		debris (array) : Array containing a combination of debris 
		
		DV (Matrix): Matrix containing the delta_v associated to each maneuver	

		debris_data (Dataframe): Dataframe containing the orbital parameters of all debris

		DV_RAAN (Matrix) - optionnal: None by default - Matrix containing the delta_v associated to each maneuver, including the ones in RAAN	

		V_tol (float) - optionnal, used when DV_RAAN is not None : None by default. Order of magnitude of tolerated dV for one mission in km/s

		t_tol (float) - optionnal, used when DV_RAAN is not None : None by default. Order of magnitude of tolerated duration for a mission in days

	Returns: 
		dV (float): Global delta v associated to the state G
		
		dT (float): Global elapsed time for RAAN alignement

	'''

	dV = 0
	dT = 0
	nb_debris = len(debris)

	SMAs = debris_data.values[debris,0]
	INCs = debris_data.values[debris,2]
	RAANs = debris_data.values[debris,3]
	# print(RAANs)

	if DV_RAAN is None :
		for l in range(nb_debris-1):
			dv = DV[debris[l],debris[l+1]]
			dt = compute_dt(SMAs[l], SMAs[l+1], INCs[l], INCs[l+1], RAANs[l], RAANs[l+1])
			RAANs[l+1:] = RAAN_evol(SMAs[l+1:],INCs[l+1:],dt*86400,RAANs[l+1:])

			#print()
			#print(RAANs)

			dV += dv
			dT += dt

		return dV,dT

	else:
		RAAN_man = np.zeros(nb_debris-1)

		for l in range(nb_debris-1):
			dv = DV[debris[l],debris[l+1]]
			dv_raan = DV_RAAN[debris[l],debris[l+1]]

			dt = compute_dt(SMAs[l], SMAs[l+1], INCs[l], INCs[l+1], RAANs[l], RAANs[l+1])
			RAANs[l+1:] = RAAN_evol(SMAs[l+1:],INCs[l+1:],dt*86400,RAANs[l+1:])

			if (dv_raan - dv)/V_tol < (dt/t_tol) :
				# If RAAN maneuver costs less in delta_v than the J2 drift in years
				dV+=dv_raan
				RAAN_man[l] = 1
			else:
				dV += dv
				dT += dt

		return dV,dT,RAAN_man


def minimal_energy_computation(ordered_debris, DV, debris_data, V_tol, t_tol, DV_RAAN = None):
	''' Function used to compute the energy associated to a state

	Arguments:
		ordered_debris (array) : array containing the debris in their natural order
		
		DV (Matrix): Matrix containing the delta_v associated to each maneuver
		
		debris_data (Dataframe): Dataframe containing the orbital parameters of all debris

		V_tol (float) : Order of magnitude of tolerated dV for one mission in km/s

		t_tol (float) : Order of magnitude of tolerated duration for a mission in days

		DV_RAAN (Matrix) - optionnal : None by default. Matrix containing the delta_v associated to each maneuver, including the ones in RAAN

	Returns: 
		E_min (float): Minimal energy produced by the optimal combination of debris 

		debris_opt (array) : Optimal debirs configuration for successive maneuvers

		RAAN_man_opt (array) - returned if DV_RAAN is not None: Array of booleans 
			--> RAAN_man[i] = 1 if there is a RAAN maneuver between debris i and i+1
			--> RAAN_man[i] = 0 otherwise
	'''

	E_min = 10**5
	debris_opt = ordered_debris
	RAAN_man_opt = np.zeros(len(debris_opt)-1)

	for debris in permutations(ordered_debris):
		if DV_RAAN is None:
			dV,dT = single_energy_computation(debris,DV,debris_data)
		else:
			dV, dT, RAAN_man = single_energy_computation(debris,DV,debris_data,DV_RAAN=DV_RAAN,V_tol=V_tol,t_tol=t_tol)

		# E = dV/V_tol + dT/t_tol
		E = dV + V_tol/t_tol*dT

		if E<E_min:
			E_min = E
			dV_opt = dV
			dT_opt = dT
			debris_opt = debris
			if DV_RAAN is not None:
				RAAN_man_opt = RAAN_man

	if DV_RAAN is None:
		return dV_opt, dT_opt, debris_opt
	else:
		return dV_opt, dT_opt, debris_opt, RAAN_man_opt



def energy_computation(G, DV, debris_data, V_tol, t_tol, DV_RAAN = None, show_grps = False):
	''' Function used to compute the energy associated to a state

	Arguments:
		G (Matrix): Actual state for which we compute the energy
		
		DV (Matrix): Matrix containing the delta_v associated to each maneuver
		
		debris_data (Dataframe): Dataframe containing the orbital parameters of all debris

		V_tol (float) : Order of magnitude of tolerated dV for one mission in km/s

		t_tol (float) : Order of magnitude of tolerated duration for a mission in days

		DV_RAAN (Matrix) - optionnal : None by default. Matrix containing the delta_v associated to each maneuver, including the ones in RAAN
		
		show_grps (bool) : False by default - If True, gives the detail of groups formed in the state

	Returns: 
		E (float): Energy associated to the state G
		
		grps (array) - optionnal : False by default. If true, gives the detail of the order of the maneuvers for each group

		RAAN_mans - optionnal : Displayed only if show_grps is True and DV_RAAN is not None. Indicates weather or not we use RAAN maneuver in a group

	'''

	nb_grp = np.size(G,1)
	if show_grps:
		grps = []
		RAAN_mans = []

	E = 0

	for n in range(nb_grp):

		grp = G[:,n]
		ordered_debris = np.nonzero(grp)[0]

		if DV_RAAN is None :
			dV, dT, debris_opt = minimal_energy_computation(ordered_debris,DV,debris_data,V_tol,t_tol)
		else:
			dV, dT, debris_opt, RAAN_man = minimal_energy_computation(ordered_debris,DV,debris_data,V_tol,t_tol,DV_RAAN=DV_RAAN)

		# std = np.std([dV/V_tol, dT/t_tol])
		std = np.std([dV, V_tol*dT/t_tol])
		E_n = dV + V_tol/t_tol*dT + std

		E += E_n
		if show_grps:
			grps.append(debris_opt)
			if DV_RAAN is not None:
				RAAN_mans.append(RAAN_man)

	if show_grps:
		if DV_RAAN is None:
			return E, grps
		else:
			return E, grps, RAAN_mans
	else:
		return E


